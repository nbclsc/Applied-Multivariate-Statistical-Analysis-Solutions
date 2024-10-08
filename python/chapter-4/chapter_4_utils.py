import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from collections import namedtuple
from scipy import stats


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

def box_cox_power_transform(x: np.ndarray, lam_min: float, lam_max: float) -> tuple[np.ndarray|float]:
    r'''
    Univariate Box-Cox power transform to find the optimal power transformation attempting to make the input data normally distributed.
    Args:
        x (np.ndarray): Array to apply the Box-Cox transformation to.
        lam_min (float): A minimum value of lambda to apply transformation to.
        lam_max (float): A maximum value of lambda to apply transformation to.
    Returns:
        tuple[np.ndarray|float]: The first element is an np.ndarray of the computed values of \ell(\lambda), from
        the expression in (4-35). The second element is the maximum value of \lambda found in the computed values
        of \ell(\lambda).
    '''
    if (x == 0).sum() > 0:
        x = x + 1e-7
    lmbda = np.linspace(start=lam_min, stop=lam_max, num=500)

    # Reshape lambda_values to be broadcastable with x (if x is a 1D vector)
    lmbda_values = lmbda.reshape(1, -1)

    # Box-Cox power transformation (4-34).
    x_lmbda = (x[:, np.newaxis]**lmbda_values - 1) / lmbda
    # When \lambda is zero, we take the natural log.
    if (lmbda == 0).sum() > 0:
        x_lmbda[:, lmbda == 0] = np.log(x)

    # Mean of the lambda transformations, \overline{x^{\lambda}} from (4-36).
    # x_lmbda_bar = np.mean(x_lmbda, axis=0)
    # Expression from (4-35).
    n = x.shape[0]
    l_lmbda = -(n/2)*np.log(np.var(x_lmbda, axis=0, ddof=0)) + (lmbda - 1)*np.log(x).sum(axis=0)

    # Pull the value of \lambda that's the maximum from the expression in (4-35).
    max_lmbda, argmax_lmbda = np.max(l_lmbda), np.argmax(l_lmbda)
    min_lmbda = np.min(l_lmbda)
    best_lmbda = lmbda[argmax_lmbda]

    lmbda_points = np.vstack([lmbda, l_lmbda]).T
    return lmbda_points, best_lmbda

def plot_box_cox_transformation(lmbda_points: np.ndarray, max_lmbda: float, var_text: str):
    r'''
    Plot the output from the Box-Cox transformation.
    Args:
        lmbda_points (np.ndarray): The computed values of \ell(\lambda), from the expression in (4-35).
        max_lmbda (float): The maximum value of \lambda found in the computed values of \ell(\lambda).
        var_text (str): String with the variable name. Used in the title of the plot.
    '''

    plt.figure()
    ax = plt.gca()
    
    # Plot the Box-Cox curve.
    plt.plot(lmbda_points[:,0], lmbda_points[:,1])
    # Plot the dashed line.
    max_point = lmbda_points[np.where(lmbda_points[:,0] == max_lmbda),:].copy().squeeze()
    plt.plot([max_point[0], max_point[0]], [np.min(lmbda_points[:,1]), max_point[1]], linestyle='--', color='black')
    # Plot the max point.
    plt.plot(max_point[0], max_point[1], marker='o', markersize=7, color='black')
    plt.annotate(fr'$\lambda={max_lmbda:.4f}$', max_point + 0.2)
    plt.xlabel(r'$\lambda$')
    plt.ylabel(r'$\ell(\lambda)$')
    plt.title(f'Box-Cox transformation: {var_text}')
    
    return plt, ax

MultivariatePowerTransformation = namedtuple('MultivariatePowerTransformation',
                                                 ['mesh_x', 'mesh_y', 'l_lmbda', 'lmbda_x', 'lmbda_y', 'argmax_x', 'argmax_y', 'l_max'])
def power_transform_2d(X: np.ndarray, lam_min: float, lam_max: float) -> MultivariatePowerTransformation:
    r'''
    Compute the bivariate power transformation from input data. Returns lots of things used for plotting the contour plot.
    Args:
        X (np.ndarray): A 2D array of data. The data should be n x p. n= number of observations. p = number of variables (should be 2).
        lam_min (float): A minimum value of \lambda to apply transformation to.
        lam_max (float): A maximum value of \lambda to apply transformation to.
    Returns:
        MultivariatePowerTransformation: namdedtuple with 8 elements. They are:
            mesh_x: Is the meshgrid values in x-direction.
            mesh_y: Is the meshgrid values in y-direction.
            l_lmbda: Are the values computed in (4-40).
            lmbda_x: Are \lambda values in x-direction.
            lmbda_y: Are \lambda values in y-direction.
            argmax_x: Is the x-direction argmax.
            argmax_y Is the y-direction argmax.
            l_max: Is the max \ell(\lambda_{1}, \lambda_{2}) value.
    '''
    # Set up 2D lambda values.
    lmbda1 = np.linspace(lam_min, lam_max, 150)
    lmbda2 = lmbda1

    # Meshgrid 1st argument is x-axis (cols) and 2nd argument is y-axis (rows).
    l1, l2 = np.meshgrid(lmbda1, lmbda2)

    # Assuming X is already defined (e.g. X = np.random.rand(n, 2) for example purposes)
    # Transformed observations using (4-34).
    x1_lmbda = (X[:, 0].reshape(-1, 1) ** lmbda1 - 1) / lmbda1  # n x 150
    x2_lmbda = (X[:, 1].reshape(-1, 1) ** lmbda2 - 1) / lmbda2  # n x 150

    # Handle the case where lmbda1 or lmbda2 is zero from (4-34).
    x1_lmbda[:, lmbda1 == 0] = np.log(X[:, 0])[:, np.newaxis]
    x2_lmbda[:, lmbda2 == 0] = np.log(X[:, 1])[:, np.newaxis]

    # Compute values for function to maximize equation (4-40).
    l_lmbda = np.zeros((len(lmbda2), len(lmbda1)))
    n = X.shape[0]  # Number of observations.
    part1 = np.sum(np.log(X[:, 0]))  # Variable X1.
    part2 = np.sum(np.log(X[:, 1]))  # Variable X2.

    # Loop through all lambda1 and lambda2 combinations.
    for i, _ in enumerate(lmbda2):
        for j, _ in enumerate(lmbda1):
            # Compute the determinant of the covariance matrix.
            cov_matrix = np.cov(x2_lmbda[:, i], x1_lmbda[:, j])
            l_lmbda[i, j] = -(n / 2) * np.log(la.det(cov_matrix)) + (lmbda1[j] - 1) * part1 + (lmbda2[i] - 1) * part2

    # Find the maximum by stacking columns (linear).
    a = np.max(l_lmbda)
    max_idx = np.argmax(l_lmbda)

    # Convert linear index to row and column.
    l2_max, l1_max = np.unravel_index(max_idx, l_lmbda.shape)
    results = MultivariatePowerTransformation(l1,
                                    l2,
                                    l_lmbda,
                                    lmbda1,
                                    lmbda2,
                                    l1_max,
                                    l2_max,
                                    a)
    return results

def plot_power_transformation_2d(pwr: MultivariatePowerTransformation) -> plt:
    r'''
    Plot the contour plot and surface plot for a 2-D power transformation.
    Args:
        MultivariatePowerTransformation: A named tuple of output from the 2-D power transformation. The elements are...
            mesh_x: Is the meshgrid values in x-direction.
            mesh_y: Is the meshgrid values in y-direction.
            l_lmbda: Are the values computed in (4-40).
            lmbda_x: Are \lambda values in x-direction.
            lmbda_y: Are \lambda values in y-direction.
            argmax_x: Is the x-direction argmax.
            argmax_y Is the y-direction argmax.
            l_max: Is the max \ell(\lambda_{1}, \lambda_{2}) value.
    Returns:
        plt: The interface to the side-by-side plot of power transformation results. The contour plot on the left and surface plot on the right.
    '''
    plt.figure(figsize=(16, 10))
    # plt.subplots_adjust(wspace=0.2)  # Increase the space between the two plots.
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 2]) 

    # Plot the contour.
    plt.subplot(gs[0])
    plt.contour( pwr.mesh_x, pwr.mesh_y, pwr.l_lmbda, levels=20)
    plt.title(r'Contour Plot of $\ell\left( \lambda_{1}, \lambda_{2} \right)$ for Snow Data')
    plt.xlabel(r'$\lambda_{1}$')
    plt.ylabel(r'$\lambda_{2}$')
    plt.plot(pwr.lmbda_x[pwr.argmax_x], pwr.lmbda_y[pwr.argmax_y], 'o')
    plt.text(pwr.lmbda_x[pwr.argmax_x], pwr.lmbda_y[pwr.argmax_y], f'({pwr.lmbda_x[pwr.argmax_x]:.4f}, {pwr.lmbda_y[pwr.argmax_y]:.4f})', 
             verticalalignment='bottom', horizontalalignment='right')

    # Plot the surface.
    ax = plt.subplot(gs[1], projection='3d')
    ax.plot_surface(pwr.mesh_x, pwr.mesh_y, pwr.l_lmbda, cmap='viridis')
    ax.set_title(r'Surface Plot of $\ell\left( \lambda_{1}, \lambda_{2} \right)$ for Snow Data')
    ax.set_xlabel(r'$\lambda_{1}$')
    ax.set_ylabel(r'$\lambda_{2}$')
    ax.set_zlabel(r'$\ell(\lambda_1, \lambda_2)$')
    ax.plot([pwr.lmbda_x[pwr.argmax_x]], [pwr.lmbda_y[pwr.argmax_y]], [pwr.l_max], 'bo', markersize=10, markerfacecolor='b')
    ax.text(pwr.lmbda_x[pwr.argmax_x], pwr.lmbda_y[pwr.argmax_y], pwr.l_max + 2, f'({pwr.lmbda_x[pwr.argmax_x]:.3f}, {pwr.lmbda_y[pwr.argmax_y]:.3f})', 
            verticalalignment='top', horizontalalignment='right')
    
    plt.tight_layout()

    return plt