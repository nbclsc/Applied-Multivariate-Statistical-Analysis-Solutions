'''
Functions related to things defined in Chapter 6.
'''
from collections import namedtuple
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la
import pandas as pd
from IPython.display import display, Math
from matplotlib.patches import Ellipse
from scipy import stats

def plot_mean_difference_confidence_ellipse(x1: np.ndarray, x2: np.ndarray, alpha: float):
    '''
    Plot the 2D confidence ellipse as found in (6-23), Result 6.2, and Example
    6.3 found in section 6.3 on page 287. Just like the other confidence regions,
    the 1D confidence intervals are the projection of the 2D ellipse onto the
    component axis.
    Args:
        x1 (np.ndarray): The first sample. A $n_{1} \times 1$ array to create a confidence ellipse for.
        x2 (np.ndarray): The second sample. A $n_{2} \times 1$ array to create a confidence ellipse for.
        n1 (int): The number of rows from the first sample for the 2-D array.
        n2 (int): The number of rows from the second sample for the 2-D array.
        p (int): The original number of features sumultaneous confidence intervals were created for.
        alpha (float):
    Returns:
        tuple[plt, ax]: The plot graphic.
    '''
    x1 = x1.copy()
    x2 = x2.copy()
    n1 = x1.shape[0]
    n2 = x2.shape[0]
    p = 2
    # assert n1 == n2, f'Input data should have same number of columns, found #x1={x1.shape[1]} and #x2={x2.shape[1]}.'
    xbar1 = np.mean(x1, axis=0)[:, np.newaxis]
    xbar2 = np.mean(x2, axis=0)[:, np.newaxis]
    xbard = xbar1 - xbar2

    S1 = np.cov(x1, rowvar=False)
    S2 = np.cov(x2, rowvar=False)

    S_pooled = ((n1-1)*S1 + (n2-1)*S2)/(n2 + n2 - 2)

    eigenvalues, eigenvectors = la.eigh(S_pooled)
    max_idx, min_idx = np.argmax(eigenvalues), np.argmin(eigenvalues)
    lmbda1, lmbda2 = eigenvalues[max_idx], eigenvalues[min_idx]
    e1, e2 = eigenvectors[:, max_idx].copy(), eigenvectors[:, min_idx].copy()

    const = ((n1 + n2 - 2)*p)/(n1 + n2 - (p + 1))
    f_val = stats.f.ppf(1-alpha, dfn=p, dfd=n1 + n2 - (p + 1))

    ell_width = np.sqrt(lmbda1)*np.sqrt((1/n1)+(1/n2))*np.sqrt(const*f_val)
    ell_height = np.sqrt(lmbda2)*np.sqrt((1/n1)+(1/n2))*np.sqrt(const*f_val)
    ell_angle = np.degrees(np.arctan2(e1[1], e1[0]))
    # print(f'{ell_angle}, width = {ell_width:.4f}, height = {ell_height:.4f}')

    plt.figure()
    ax = plt.gca()
    ellipse = Ellipse(xy=xbard,
                    width=2*ell_width,
                    height=2*ell_height,
                    angle=ell_angle,
                    fill=False)
    ax.add_patch(ellipse)
    for i in [-1, 1]:
        plt.quiver(xbard[0],
                xbard[1],
                e1[0] * ell_width * i,
                e1[1] * ell_width * i,
                angles='xy',
                scale_units='xy',
                scale=1
                )
        plt.quiver(xbard[0],
                xbard[1],
                e2[0]* ell_height * i,
                e2[1]* ell_height * i,
                angles='xy',
                scale_units='xy',
                scale=1)
    return plt, ax

def manova_obs_breakdown(df: pd.DataFrame, trt_col: str, var_col: str) -> namedtuple:
    '''
    Breakdown observations for an input variable into mean, treatment, and residual components.
    Args:
        df (pd.DataFrame): Input data with a column for treatments and columns for each variable.
        trt_col (str): The column with the treatments (groups).
        var_col (str): The variable breakdown.
    Return:
        namedtuple: The named tuple has elements: 'Variable', 'Obs', 'Mean', 'TreatmentEffect',
        'Residual'.
    '''
    # Store the output for the observations breakdown for a given variable.
    ObsBreakdown = namedtuple('ObsBreakdown', ['Variable', 'Obs', 'Mean', 'TreatmentEffect', 'Residual'])
    df = df[[trt_col, var_col]].copy()

    # The number of groups.
    g = df[trt_col].nunique()
    # The max number of observations within a group.
    n_max = df[trt_col].value_counts().max()

    # Setup the arrays to store the data.
    obs_a = np.full([g, n_max], np.nan)
    global_mean_a = np.full([g, n_max], np.nan)
    trt_effect_a = np.full([g, n_max], np.nan)
    residual_a = np.full([g, n_max], np.nan)

    global_mean = np.mean(df[var_col])
    for i, (_, r) in enumerate(df.groupby(trt_col)):
        # Store the group data for a variable in a numpy array.
        grp_a = r[var_col].to_numpy()
        n = grp_a.shape[0]
        grp_mean = np.mean(grp_a)

        # Compute the breakdown.
        obs_a[i, :len(grp_a)] = grp_a
        global_mean_a[i, :len(grp_a)] = np.repeat(global_mean, n)
        trt_effect_a[i, :len(grp_a)] = np.repeat(grp_mean, n) - np.repeat(global_mean, n)
        residual_a[i, :len(grp_a)] = grp_a - np.repeat(grp_mean, n)
        
    return ObsBreakdown(Variable=var_col,
                        Obs=obs_a,
                        Mean=global_mean_a,
                        TreatmentEffect=trt_effect_a,
                        Residual=residual_a)

def create_array_text(a: np.ndarray) -> str:
    '''
    Create a text string with the latex code to generate an array.
    Args:
        a (np.ndarray): Contains array data to generate a latex array for.
    Returns:
        str: Latex code for the input array.
    '''
    g, n = a.shape
    array_col_num = 'r' * n
    start_array_latex = fr'\left[\begin{{array}}{{ {array_col_num} }}'
    end_array_latex = r'\end{array}\right]'
    array_str = start_array_latex
    for i in range(g):
        array_str += ' & '.join(['' if np.isnan(v) else str(int(v)) for v in a[i,:]])
        array_str += r' \\ '
    array_str += end_array_latex
    return array_str

def display_manova_obs_breakdown(data: namedtuple, spacing: list[str]) -> None:
    '''
    Take a breakdown of observations into components and display them using latex.
    Args:
        data (namedtuple): Contains breakdown of observations into components.
        Row for each group, column for each observation.
        spacing list[str]: A list with 4 string elements.
        Something like, ['0.5cm','2.0cm','2.2cm','2.5cm'].
    '''
    assert len(spacing)==4, 'Spacing must have 4 string elements.'
    obs_latex = create_array_text(data.Obs)
    mean_latex = create_array_text(data.Mean)
    trt_effect_latex = create_array_text(data.TreatmentEffect)
    residual_latex = create_array_text(data.Residual)
    display(Math(fr'\text{{Variable: {data.Variable}}}'))
    display(Math(f'{obs_latex}'
                 ' = '
                 f'{mean_latex}'
                 ' + '
                 f'{trt_effect_latex}'
                 ' + '
                 f'{residual_latex}'
                 ))
    display(Math(fr'\hspace{{ {spacing[0]} }}\text{{(observation)}}'
                 fr'\hspace{{ {spacing[1]} }}\text{{(mean)}}'
                 fr'\hspace{{ {spacing[2]} }}\text{{(treatment effect)}}'
                 fr'\hspace{{ {spacing[3]} }}\text{{(residual)}}'))
