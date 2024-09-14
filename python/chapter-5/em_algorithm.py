"""
Functions and class necessary to compute EM algorithm estimates for
the sample mean vector mu and the sample covariance matrix S.
"""
import numpy as np
from collections import namedtuple


class Partition(object):
    """
    Partition the mean vector and covariance matrix based on an input, "missing", vector.
    The missing vector should be p-dimensions.
    Also takes care of reordering variables in the mean vector and covariance matrix so
    they're in block mu1 and S_11.
    Don't really need this for example 5.13, but I did it anyway.
    """

    XBarBlock = namedtuple('XBarBlock', ['xbar1' , 'xbar2'])
    SBlock = namedtuple('SBlock', ['S11', 'S12', 'S21', 'S22'])
    def __init__(self, xbar: np.ndarray, S: np.ndarray, missings: np.array):
        self.xbar = xbar
        self.S = S
        self.missings = missings
        self.missing_vars = np.where(self.missings)[0]  # Row idx values with missing.
        self.n_missing = self.missing_vars.shape[0]  # How many variables with missing.
    
    def swap_value(self, x: np.ndarray) -> np.ndarray:
        """
        Swap columns in the input array so values true in the missing vector are at the
        first column(s) (starting at column index 0).
        Args:
            x (np.ndarray): An input array with missing data that needs its columns swapped.
        Returns:
            np.ndarray: An array where columns with any missing values are moved to the
            columns on the most-left positions.
        """
        x, missings = x.copy(), self.missings.copy()
        order = list(range(x.shape[1]))
        # Loop through index list of rows with missing values.
        for v in self.missing_vars:
            if v > 0:
                # Any columns before this one that are not missing.
                possible_places = [c for c in range(v) if ~missings[c]]
                if len(possible_places) == 0:
                    # Nowhere to swap.
                    continue
                #  The closest place to column 0 to swap with.
                free_spot = min(possible_places)
                print(f'\tthen swap column {v} with column {free_spot}')
                # Keep track of the swap.
                order[free_spot], order[v] = order[v], order[free_spot]
                # Swap elements in target array.
                x[:, [v, free_spot]] = x[:, [free_spot, v]]
                missings[[v, free_spot]] = missings[[free_spot, v]]
        return x, missings, order
    
    def swap_row_columns(self, x: np.ndarray) -> tuple[np.ndarray, dict[str, list[int]]]:
        """
        Swap columns in the input array so values true in the missing vector are at the
        first row(s)/column(s) (starting at row index [0, 0]). For a matrix with missing
        values in rows and columns, rearrange the matrix so missing data appears in the
        upper left of the matrix. This is so we can partition the matrix with those
        variables together.
        Args:
            x (np.ndarray): A matrix with missing data for values in rows and columns.
        Returns:
            np.ndarray: A matrix where the missing data positions are moved to the
            top-left.
        """
        if self.n_missing == 0:
            print(f'Observation has no missing data.')
            return x
        else:
            print(f'Observation has {self.n_missing} missing variables.')

        x = x.copy()
        # Swap rows.
        x, rmissing_arrange, rorder = self.swap_value(x.T)
        x = x.copy().T
        # Swap columns.
        x, cmissing_arrange, corder = self.swap_value(x)
        for i in range(self.n_missing):
            assert cmissing_arrange[i], f'Variable column {i} not moved!'
        return x, {'col-order': corder, 'row-order': rorder}
    
    def partition_mean(self) -> XBarBlock:
        """
        Partition the mean vector into 2 blocks.
        Args:
            None
        Returns:
            XBarBlock: The named tuple containing the partitioned mean vector.
            The first partition contains the variables with missing data.
        """
        xbar = self.xbar.copy()
        # Swap any rows that need it.
        xbar, _, rorder = self.swap_value(xbar.T)
        xbar = xbar.T
        stop_idx = self.n_missing

        xbar1 = xbar[0:stop_idx]
        xbar2 = xbar[stop_idx:]
        return self.XBarBlock(xbar1, xbar2), {'row-order': rorder}

    def partition_covar(self) -> SBlock:
        """
        Partition the covariance matrix into 4 blocks.
        Ags:
            None
        Returns:
            SBlock: The named tuple containing the partitioned covariance matrix.
            The S11 block contains the variables with missing data.
        """
        S = self.S.copy()
        S, order = self.swap_row_columns(S)
        stop_idx = self.n_missing
        S11 = S[0:stop_idx, 0:stop_idx]
        S12 = S[0:stop_idx, stop_idx:]
        S21 = S[stop_idx:, 0:stop_idx]
        S22 = S[stop_idx:, stop_idx:]
        return self.SBlock(S11, S12, S21, S22), order
    
def extract_notna_submatrix(x: np.ndarray) -> np.ndarray:
    """
    Used in the EM algorithm when computing values for T2_list, where we have a matrix
    with certain rows and columns missing and only a block matrix
    contains data. Extract the block of data without missing values.
    For example, the input 4 x 4 matrix:
    [NA NA NA NA]
    [NA NA NA NA]
    [NA NA 4  5 ]
    [NA NA 6  7 ]
    would return the 2 x 2 output matrix:
    [4 5]
    [6 7]
    Args:
        x (np.ndarray): An input matrix. It's going to be square and p x p.
    Returns:
        np.ndarray: The block matrix without missing data.
    """
    notna_cols = ~np.all(np.isnan(x), axis=0)
    notna_rows = ~np.all(np.isnan(x), axis=1)
    return x[notna_rows][:, notna_cols]

def compute_xxt(x: np.ndarray, predictions: dict[str, np.ndarray]) -> np.ndarray:
    """
    Compute X * X^{\prime} for an input p x 1 vector x.
    Replace NA values with values predicted using the EM algorithm.
    Args:
        x (np.ndarray): An observation of data. Basically, a p x 1 array.
        predictions (dict[str, np.ndarray]): An dictionary containing the estimates from EM.
        The key values are 'a', 'b', and 'c'. Here, we only use 'b' and 'c'. The key 'b'
        corresponds to the predictions for crossterms (like x_{1}x_{2}) when both values are missing. The key 'c' corresponds to the predictions when only one of the crossterm values (like x_{1}x_{2}) are missing. 
    Returns:
        np.ndarray: A p x p matrix that has contains the results of X * X^{\prime},
        where the NA values in the results are filled in using the EM predictions.
    """
    dotted = x @ x.T
    s = np.hstack((predictions['b'], predictions['c']))
    t = np.hstack((predictions['c'].T, extract_notna_submatrix(dotted)))
    return np.vstack((s,t))

def compute_em_estimates(X: np.ndarray) -> namedtuple:
    """
    Use the EM algorithm to estimate the mean and covariance matrix when we have missing at random (MAR) data.
    Args:
        X (np.ndarray): n x p matrix of data. n observations with p measurements for each variable.
    Returns:
        namedtuple: The estimates for the vector mu and covariance vector sigma.
    """
    EMEstimates = namedtuple('EMEstimates', ['mu' , 'sigma'])
    X = X.copy()
    n, p = X.shape
    X_missing = np.isnan(X)
    xbar_til = np.nanmean(X, axis=0, ).reshape(p, 1)
    X_til = np.nan_to_num(X, nan=xbar_til.squeeze())
    S_til = np.cov(X_til, ddof=0, rowvar=False)

    T1_pre = X.copy()
    T2_list = list()
    for i, row in enumerate(X_missing):
        X_row = X[i,:].reshape(p, 1).copy()
        if row.sum() > 0:

            # Split the mean vector and covariance matrix into blocks.
            i_missing = X_missing[i]
            part = Partition(xbar_til, S_til, i_missing)
            xbar_til_p, _ = part.partition_mean()
            S_til_p, S_order = part.partition_covar()
            predictions = dict()

            # Estimate single terms for missing data found in (5-38) on page 252.
            mu1 = xbar_til_p.xbar1
            s12 = S_til_p.S12
            s22_inv = np.linalg.inv(S_til_p.S22)
            x2 = X[i,~i_missing].reshape(~i_missing.sum(), 1)
            mu2 = xbar_til[~i_missing]

            x1 = mu1 + s12 @ s22_inv @ (x2 - mu2)
            predictions['a'] = x1

            # Estimate cross(squared) terms, where both are missing found in (5-39) on page 252.
            s11 = S_til_p.S11
            s21 = S_til_p.S21

            x1_sq = s11 - s12 @ s22_inv @ s21 + x1 @ x1.T
            predictions['b'] = x1_sq
            
            # Estimate cross (squared) terms, where one is missing and one is not
            #  found in (5-39) on the top of page 253.
            # cross_term_est = x1 * X[i,~i_missing]
            cross_term_est = x2.squeeze() * x1
            # print(cross_term_est)
            predictions['c'] = cross_term_est
            T1_pre[i, row] = x1.squeeze()

            # Fill missings with predictions.
            xxt = compute_xxt(X_row, predictions)
            restored_arr = xxt[np.argsort(S_order['row-order']), :]  # Restore original row order.
            restored_arr = restored_arr[:, np.argsort(S_order['col-order'])]  # Restore original column order
            T2_list.append(restored_arr)
        else:
            T2_list.append(X_row @ X_row.T)

    T1_til = T1_pre.sum(axis=0).reshape(p, 1)
    T2_til = sum(T2_list)

    # Estimation of mu and sigma found in (5-40) on page 253.
    mu_til = T1_til / n
    sigma_til = T2_til / n - mu_til @ mu_til.T

    return EMEstimates(mu=mu_til, sigma=sigma_til)
