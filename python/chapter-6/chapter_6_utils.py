'''
Functions related to things defined in Chapter 6.
'''
from collections import namedtuple
from dataclasses import dataclass, field
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

def one_way_manova_obs_breakdown(df: pd.DataFrame, trt_col: str, var_col: str) -> namedtuple:
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
    def latex_cell_value(v: float) -> str:
        if v < 0:
            return str(int(v))
        else:
            return r'\phantom{-} ' + str(int(v))
    g, n = a.shape
    array_col_num = 'r' * n
    start_array_latex = fr'\left[\begin{{array}}{{ {array_col_num} }}'
    end_array_latex = r'\end{array}\right]'
    array_str = start_array_latex
    for i in range(g):
        array_str += ' & '.join(['' if np.isnan(v) else latex_cell_value(v) for v in a[i,:]])
        array_str += r' \\ '
    array_str += end_array_latex
    return array_str

def display_1way_manova_obs_breakdown(data: namedtuple, spacing: list[str]) -> None:
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

def display_2way_manova_obs_breakdown(data: namedtuple, spacing: list[str]) -> None:
    '''
    Take a breakdown of observations into components and display them using latex.
    Args:
        data (namedtuple): Contains breakdown of observations into components.
        Row for each group, column for each observation.
        spacing list[str]: A list with 4 string elements.
        Something like, ['0.5cm','2.0cm','2.2cm','2.5cm'].
    '''
    assert len(spacing)==5, 'Spacing must have 4 string elements.'
    obs_latex = create_array_text(data.Obs)
    mean_latex = create_array_text(data.Mean)
    trt_effect1_latex = create_array_text(data.TreatmentEffect1)
    trt_effect2_latex = create_array_text(data.TreatmentEffect2)
    residual_latex = create_array_text(data.Residual)
    display(Math(fr'\text{{Variable: {data.Variable}}}'))
    display(Math(f'{obs_latex}'
                 ' = '
                 f'{mean_latex}'
                 ' + '
                 f'{trt_effect1_latex}'
                 ' + '
                 f'{trt_effect2_latex}'
                 ' + '
                 f'{residual_latex}'
                 ))
    display(Math(fr'\hspace{{ {spacing[0]} }}\text{{(observation)}}'
                 fr'\hspace{{ {spacing[1]} }}\text{{(mean)}}'
                 fr'\hspace{{ {spacing[2]} }}\text{{(treatment 1 effect)}}'
                 fr'\hspace{{ {spacing[3]} }}\text{{(treatment 2 effect)}}'
                 fr'\hspace{{ {spacing[4]} }}\text{{(residual)}}'))

def compute_manova_ss_matrices(X: np.ndarray) -> np.ndarray:
    r'''
    Use along with manova_obs_breakdown. Can use this to compute B = SS_tr or W = SS_res.
    The manova_obs_breakdown function creates the components in,
    $\textbf{X}_{v} = \textbf{M}_{v} + \textbf{T}_{v} + \textbf{E}_{v}$, where $v$ identifies
    which measurement we're looking at.
    For the sum of squares and cross-products in the MANOVA table we need a matrix result,
    that comes from some block computations. An example of what that looks like for the
    treatment effect is below. The The Hadamard product, denoted by $\circ$, performs
    elementwise multiplication of the entries of matrices of the same dimensions. The code
    for the computations implimented here is similar to the first part.

    $$
    \left[
        \begin{array}{cc}
            \text{sum}(\textbf{T}_{1} \circ \textbf{T}_{1}) & \text{sum}(\textbf{T}_{1} \circ \textbf{T}_{2}) \\
            \text{sum}(\textbf{T}_{2} \circ \textbf{T}_{1}) & \text{sum}(\textbf{T}_{2} \circ \textbf{T}_{2})
        \end{array}
    \right]
    =
    \left[
        \begin{array}{cc}
            \text{tr}(\textbf{T}_{1}^{\prime} \textbf{T}_{1}) & \text{tr}(\textbf{T}_{1}^{\prime} \textbf{T}_{2}) \\
            \text{tr}(\textbf{T}_{2}^{\prime} \textbf{T}_{1}) & \text{tr}(\textbf{T}_{2}^{\prime} \textbf{T}_{2})
        \end{array}
    \right]
    $$
    Args:
        a1 (np.ndarray): Should be either a treatment effect or residual matrix created
        by the manova_obs_breakdown function. This should be consistent with the a2 parameter.
        a2 (np.ndarray): Should be either a treatment effect or residual matrix created
        by the manova_obs_breakdown function. This should be consistent with the a1 parameter.
    Returns:
        np.ndarray: Both inputs are placed into 2 by 2 block matrices. In the result, each
        block is element-wise multiplied and summed.
    '''
    # Stack the two matrices. The first dim is the variable. The rest is the 3x5 matrix of data.
    # X = np.stack([a1, a2])
    # X = np.nan_to_num(X)
    # For [1,2], axis 1 is groups. Axis 2 is observation within group.
    Y = np.tensordot(X, X, axes=([1, 2], [1, 2]))
    return Y

def diplay_1way_manova_table(B: np.ndarray, W: np.ndarray, T: np.ndarray, nl: list[int], g: int) -> None:
    '''
    Display the MANOVA table.
    Args:
        B (np.ndarray): Between sample sum of squares and cross-product matrix.
        W (np.ndarray): Within sum of squares and cross-products matrix.
        T (np.ndarray): The total(corrected) sum of squares and cross-products matrix.
        nl (list[int]): List containing the number of observations in each group.
        g (int): The number of groups.
    '''
    n = sum(nl)
    assert len(nl) == g, f'Number of groups ({g}) and length of nl ({len(nl)}) differ.'
    display(Math(r'\begin{array}{lll}'
             r'\text{Source} & \text{Matrix of sum of squares} &  \\'
             r'\text{of variation} & \text{and cross products} & \text{Degrees of freedom} \\'
             r'\hline \\'
             r'\text{Treatment} & '
             f'{create_array_text(B)} & '
             fr'{g} - 1 = {g - 1} \\ \\'
             r'\text{Residual} & '
             f'{create_array_text(W)} &'
             fr'{" + ".join([str(ni) for ni in nl])} - {g} = {n - g} \\ \\'
             r'\hline \\'
             r'\text{Total (corrected)} & '
             f'{create_array_text(T)} & '
             f'{(n - 1)}'
             r'\end{array}'
             ))
    
def diplay_2way_manova_table(B1: np.ndarray, B2: np.ndarray, W: np.ndarray, T: np.ndarray, nl: list[int], g: int, b: int) -> None:
    '''
    Display the 2-Way MANOVA table. NO REPLICATION!
    Args:
        B (np.ndarray): Between sample sum of squares and cross-product matrix.
        W (np.ndarray): Within sum of squares and cross-products matrix.
        T (np.ndarray): The total(corrected) sum of squares and cross-products matrix.
        nl (list[int]): List containing the number of observations in each group.
        g (int): The number of groups.
    '''
    n = sum(nl)
    # assert len(nl) == g, f'Number of treatment 1 levels ({g}) and length of nl ({len(nl)}) differ.'
    display(Math(r'\begin{array}{lll}'
             r'\text{Source} & \text{Matrix of sum of squares} &  \\'
             r'\text{of variation} & \text{and cross products} & \text{Degrees of freedom} \\'
             r'\hline \\'
             r'\text{Treatment 1} & '
             f'{create_array_text(B1)} & '
             fr'{g} - 1 = {g - 1} \\ \\'
             r'\text{Treatment 2} & '
             f'{create_array_text(B2)} & '
             fr'{b} - 1 = {b - 1} \\ \\'
             r'\text{Residual} & '
             f'{create_array_text(W)} &'
             fr'({g} - {1})({b} - {1}) = {(g-1)*(b-1)} \\ \\'
             r'\hline \\'
             r'\text{Total (corrected)} & '
             f'{create_array_text(T)} & '
             f'{g}({b}) - {1} = {g*b - 1}'
             r'\end{array}'
             ))
