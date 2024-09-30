'''
Functions related to things defined in Chapter 5.
'''
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la
import pandas as pd
from matplotlib.patches import Ellipse
from scipy import stats


def plot_confidence_ellipse(x: np.ndarray | pd.DataFrame, n: int, p: int, alpha: float):
    '''
    Plot the 2D confidence ellipse as found in (5-18), (5-19), and Example
    5.3 found in section 5.4 on page 220. The 1D confidence intervals are
    the projection of the 2D ellipse onto the component axis.
    Args:
        x (np.ndarray): A 2-dimensional array to create a confidence
        ellipse for.
        n (int): The number of rows for the 2-D array.
        p (int): The original number of features sumultaneous confidence intervals were created for.
        alpha (float):
    Returns:
        tuple[plt, ax]: The plot graphic.
    '''
    if isinstance(x, pd.DataFrame):
        x = x.to_numpy()
    x = x.copy()
    assert n == x.shape[0], f'Numer of rows and parameter n do not match. Paramerter n = {n}. Rows in x = {x.shape[0]}'
    assert x.shape[1] == 2, f'Input matrix should have two columns, found {x.shape[1]}.'
    xbar = np.mean(x, axis=0)

    eigenvalues, eigenvectors = la.eigh(np.cov(x.T))
    max_idx, min_idx = np.argmax(eigenvalues), np.argmin(eigenvalues)
    lmbda1, lmbda2 = eigenvalues[max_idx], eigenvalues[min_idx]
    e1, e2 = eigenvectors[:, max_idx].copy(), eigenvectors[:, min_idx].copy()

    const = ((n - 1)*p)/(n*(n - p))
    f_val = stats.f.ppf(1-alpha,p,n-p)

    ell_width = np.sqrt(lmbda1)*np.sqrt(const*f_val)
    ell_height = np.sqrt(lmbda2)*np.sqrt(const*f_val)
    ell_angle = np.degrees(np.arctan2(e1[1], e1[0]))
    # print(f'{ell_angle}, width = {ell_width:.4f}, height = {ell_height:.4f}')

    plt.figure()
    ax = plt.gca()
    ellipse = Ellipse(xy=xbar,
                    width=2*ell_width,
                    height=2*ell_height,
                    angle=ell_angle,
                    fill=False)
    ax.add_patch(ellipse)
    for i in [-1, 1]:
        plt.quiver(xbar[0],
                xbar[1],
                e1[0] * ell_width * i,
                e1[1] * ell_width * i,
                angles='xy',
                scale_units='xy',
                scale=1
                )
        plt.quiver(xbar[0],
                xbar[1],
                e2[0]* ell_height * i,
                e2[1]* ell_height * i,
                angles='xy',
                scale_units='xy',
                scale=1)
    return plt, ax

@staticmethod
def simult_conf_int(x: np.ndarray, alpha: float, p: int=None) -> np.ndarray:
    '''
    Compute the simultaneous confidence intervals.
    '''
    n, p_data = x.shape
    if not p:
        p = p_data
    xbar = np.mean(x, axis=0).reshape(p_data, 1)
    S = np.diag(np.cov(x, rowvar=False)).reshape(p_data,1)

    const = ((n-1)*p)/(n-p)
    f_crit = const*stats.f.ppf(1-alpha, dfn=p, dfd=n-p)
    ci = xbar + np.array([-1, 1])* np.sqrt((S/n)*f_crit)
    return ci

@staticmethod
def bonferroni_conf_int(x: np.ndarray, m: int, alpha: float) -> np.ndarray:
    '''
    Compute the Bonferroni confidence intervals.
    '''
    n, p = x.shape
    xbar = np.mean(x, axis=0).reshape(p, 1)
    S = np.diag(np.cov(x, rowvar=False)).reshape(p,1)

    t_crit = stats.t.ppf(1-(alpha/(2*m)), df=n-1)
    ci = xbar + np.array([-1, 1])*t_crit*np.sqrt(S/n)
    return ci

def plot_control_ellipse(x: np.ndarray | pd.DataFrame, alpha: float):
    '''
    Plot the control ellipse.
    '''
    if isinstance(x, pd.DataFrame):
        x = x.to_numpy()
    x = x.copy()
    _, p = x.shape
    S = np.cov(x, rowvar=False)
    xbar = np.mean(x, axis=0).reshape(p, 1)

    eigenvalues, eigenvectors = la.eigh(S)
    max_idx, min_idx = np.argmax(eigenvalues), np.argmin(eigenvalues)
    lmbda1, lmbda2 = eigenvalues[max_idx], eigenvalues[min_idx]
    e1, e2 = eigenvectors[:, max_idx].copy(), eigenvectors[:, min_idx].copy()

    chi2_val = stats.chi2(df=p).ppf(1 - alpha)

    ell_width = np.sqrt(lmbda1)*np.sqrt(chi2_val)
    ell_height = np.sqrt(lmbda2)*np.sqrt(chi2_val)
    ell_angle = np.degrees(np.arctan2(e1[1], e1[0]))

    plt.figure()
    ax = plt.gca()
    ellipse = Ellipse(xy=xbar,
                    width=2*ell_width,
                    height=2*ell_height,
                    angle=ell_angle,
                    fill=False)
    ax.add_patch(ellipse)
    for i in [-1, 1]:
        plt.quiver(xbar[0],
                xbar[1],
                e1[0] * ell_width * i,
                e1[1] * ell_width * i,
                angles='xy',
                scale_units='xy',
                scale=1
                )
        plt.quiver(xbar[0],
                xbar[1],
                e2[0] * ell_height * i,
                e2[1] * ell_height * i,
                angles='xy',
                scale_units='xy',
                scale=1)
    return plt, ax

def plot_future_control_ellipse(x: np.ndarray | pd.DataFrame, n: int, p: int, alpha: float):
    '''
    Plot the 2D control ellipse for future observations as found above (5-34) and Example
    5.12 found in section 5.6 on page 248. The 1D future confidence intervals are
    the projection of the 2D ellipse onto the component axis.
    Args:
        x (np.ndarray): A 2-dimensional array to create a future observation confidence
        ellipse for.
        n (int): The number of rows for the 2-D array.
        p (int): The original number of features sumultaneous confidence intervals were created for.
        alpha (float):
    Returns:
        tuple[plt, ax]: The plot graphic.
    '''
    if isinstance(x, pd.DataFrame):
        x = x.to_numpy()
    x = x.copy()
    assert n == x.shape[0], f'Numer of rows and parameter n do not match. Paramerter n = {n}. Rows in x = {x.shape[0]}'
    assert x.shape[1] == 2, f'Input matrix should have two columns, found {x.shape[1]}.'
    xbar = np.mean(x, axis=0)

    eigenvalues, eigenvectors = la.eigh(np.cov(x.T))
    max_idx, min_idx = np.argmax(eigenvalues), np.argmin(eigenvalues)
    lmbda1, lmbda2 = eigenvalues[max_idx], eigenvalues[min_idx]
    e1, e2 = eigenvectors[:, max_idx].copy(), eigenvectors[:, min_idx].copy()

    const = ((n + 1)*(n - 1)*p)/(n*(n - p))
    f_val = stats.f.ppf(1-alpha,p,n-p)

    ell_width = np.sqrt(lmbda1)*np.sqrt(const*f_val)
    ell_height = np.sqrt(lmbda2)*np.sqrt(const*f_val)
    ell_angle = np.degrees(np.arctan2(e1[1], e1[0]))

    plt.figure()
    ax = plt.gca()
    ellipse = Ellipse(xy=xbar,
                    width=2*ell_width,
                    height=2*ell_height,
                    angle=ell_angle,
                    fill=False)
    ax.add_patch(ellipse)
    for i in [-1, 1]:
        plt.quiver(xbar[0],
                xbar[1],
                e1[0] * ell_width * i,
                e1[1] * ell_width * i,
                angles='xy',
                scale_units='xy',
                scale=1
                )
        plt.quiver(xbar[0],
                xbar[1],
                e2[0]* ell_height * i,
                e2[1]* ell_height * i,
                angles='xy',
                scale_units='xy',
                scale=1)
    return plt, ax

def my_q_q_plot(x: np.ndarray):
    '''
    Create a Q-Q plot, for an input variable.
    '''
    x = x.copy()
    assert len(x.shape) == 1, "Wrong shape."
    x.sort()
    n = x.shape[0]
    prob = (np.arange(x.size)+1 - 0.5) / n
    quant = stats.norm.ppf(prob)

    return plt.scatter(quant, x, facecolors='none', edgecolors='royalblue')

def ppcc_simulation(n: int, num_simulations: int, q_levels: list[float]):
    '''
    Function to perform the simulation for correlation coefficient test for
    normality from the paper 'Probability Plot Correlation Coefficient Test
    for Normality' by Filliben. Just like in Table 4.2 on page 181.
    '''
    ppcc_values = list()
    for i in range(num_simulations):
        # Generate a random sample from a normal distribution.
        sample_data = stats.norm.rvs(size=n)
        sample_data.sort()

        # Theoretical quantiles from a normal distribution.
        probs = [((i+1) - 0.50)/n for i in range(n)]
        theoretical_quantiles = stats.norm.ppf(probs)

        # Compute the correlation coefficient.
        correlation_coefficient = np.corrcoef(sample_data, theoretical_quantiles)[0,1]
        
        # Store the PPCC value.
        ppcc_values.append(correlation_coefficient)
    
    # Determine critical values for common significance levels.
    critical_values = np.quantile(ppcc_values, q_levels)

    # Return the distribution of PPCC values and critical values.
    return ppcc_values, critical_values
