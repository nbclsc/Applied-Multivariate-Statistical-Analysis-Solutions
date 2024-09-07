import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la
import pandas as pd
from matplotlib.patches import Ellipse
from scipy import stats


class MultivariateFunctionsCh5(object):
    '''
    Functions related to things defined in Chapter 5.
    '''

    @staticmethod
    def plot_confidence_ellipse(x: np.ndarray, n: int, p: int, alpha: float):
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
        x = x.copy()
        assert n == x.shape[0], f'Numer of rows and parameter n do not match. Paramerter n = {n}. Rows in x = {x.shape[0]}'
        assert x.shape[1] == 2, f'Input matrix should have two columns, found {x.shape[1]}.'
        xbar = np.mean(x, axis=0)
        if isinstance(x, pd.DataFrame):
            xbar = xbar.to_numpy()

        eigenvalues, eigenvectors = la.eig(np.cov(x.T))
        max_idx, min_idx = np.argmax(eigenvalues), np.argmin(eigenvalues)
        lmbda1, lmbda2 = eigenvalues[max_idx], eigenvalues[min_idx]
        e1, e2 = eigenvectors[:, max_idx].copy(), eigenvectors[:, min_idx].copy()

        const = (p*(n - 1))/(n*(n - p))
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
    def simult_conf_int(x: np.ndarray, alpha: float) -> np.ndarray:
        '''
        Compute the simultaneous confidence intervals.
        '''
        n, p = x.shape
        xbar = np.mean(x, axis=0).reshape(p, 1)
        S = np.diag(np.cov(x, rowvar=False)).reshape(p,1)

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
