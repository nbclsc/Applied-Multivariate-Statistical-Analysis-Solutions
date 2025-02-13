'''
Compute simultaneous Bonferroni confidence intervals for two-way MANOVA.
'''
from collections import namedtuple
from itertools import combinations
from scipy import stats
from IPython.display import display, Math
import numpy as np
from TwoWayManova import TwoWayManova


class BonferroniSimultaneousConfidenceIntervals(object):
    '''
    Compute and display Bonferroni simultaneous confidence intervals for two-way MANOVA.
    '''

    def __init__(self, manova: TwoWayManova, alpha: float):
        self.manova = manova
        self._design = manova._design
        self.alpha = alpha
        self.bonferroni_output: namedtuple = self._bonferroni_two_way()

    @staticmethod
    def _create_contrast_matrix(n: int) -> tuple[list[tuple], np.ndarray]:
        '''
        Args:
            n (int): The number of factor levels.
        Returns:
            tuple[list[tuple], np.ndarray]:
        '''
        pair_indices = list(combinations(range(n), 2))

        # Build the contrast matrix:
        # For each pair (i, j), create a contrast vector with 0 at i, -1 at j, and 0 elsewhere.
        contrast_list = []
        for i, j in pair_indices:
            contrast = np.zeros(n)
            contrast[i] = 1
            contrast[j] = -1
            contrast_list.append(contrast)
        return pair_indices, np.vstack(contrast_list)

    def _bonferroni_two_way(self) -> namedtuple:
        '''
        The formulas used here is from (6-70) on pages 317 and 318.
        Slice a row out of the tau tensor. Make sure they're columns. There's a column for every variable. The variable is the first arg in shape.
        Returns:
            namedtuple:
        '''
        g = self._design.metadata.g
        b = self._design.metadata.b
        p = self._design.metadata.p
        n = self._design.metadata.n

        def _compute_bonf_ci(target_levels: int, other_levels: int, factor: np.ndarray) -> list[np.ndarray]:
            ci_dict = {}
            v = target_levels * other_levels * (n - 1)
            t_val = stats.t.ppf(1-(self.alpha/(p*target_levels*(target_levels-1))), df=v)
            pairs, C = self._create_contrast_matrix(target_levels)
            ci_dict.update({'pairs': pairs})
            diff_matrix = C @ factor.T
            E = np.diag(self.manova.manova.W)*np.ones((target_levels-1,p))

            # Loop through each of the p variables and compute the CI.
            for i in range(diff_matrix.shape[1]):
                height = diff_matrix[:, i].shape[0]
                tau_diff = diff_matrix[:, i].reshape(height,1)
                E_ii = E[0,i] * np.ones((height,1))
                ci = tau_diff + np.array([-1, 1]) * t_val * np.sqrt(E_ii * 1/v * 2/(other_levels*n))
                ci_dict.update({f'x{i+1}': ci})
            return ci_dict
    
        BonferroniCI = namedtuple('BonferroniCI', ['Factor1', 'Factor2'])
        return BonferroniCI(
            # All columns are the same for factor1, so pick first column.
            Factor1=_compute_bonf_ci(g, b, self.manova.x_breakdown.anova_values.Factor1Effect[:,:,0]),
            # All rows are the same for factor2, so pick first row.
            Factor2=_compute_bonf_ci(b, g, self.manova.x_breakdown.anova_values.Factor2Effect[:,0,:])
            )

    def display_two_way_manova_bonf(self) -> None:
        '''
        Create a nice display for the confidence intervals.
        '''
        for factor_name, variables in self.bonferroni_output._asdict().items():
            variables = variables.copy()
            print(factor_name)
            pairs = variables.pop('pairs')
            # print(pairs)
            for i, cis in enumerate(variables.values()):
                factor_count_label = {1: 'g', 2: 'b'}
                for j in range(cis.shape[0]):
                    ci = cis[j,:]
                    ell = pairs[j][0]+1
                    m = pairs[j][1]+1
                    display(Math(fr'x_{i+1}:'))
                    display(Math(fr'\qquad\qquad\ell = {ell}, m = {m}, i = {i+1}'))
                    display(Math(fr'\qquad\qquad\tau_{{ {ell}, {i+1} }} - \tau_{{ {m}, {i+1} }}'
                                 r'\pm '
                                 r't_{v}'
                                 fr'\left( \frac{{ \alpha }}{{ p{factor_count_label[i+1]}({factor_count_label[i+1]}-1) }} \right)'
                                 fr'\sqrt{{ \frac{{ E_{{ {i+1}, {i+1} }} }}{{ v }} \frac{{ 2 }}{{ {factor_count_label[2-i]}n }} }}'
                                '='
                                f'[{ci[0]:.4f}, {ci[1]:.4f}]'
                                ))
