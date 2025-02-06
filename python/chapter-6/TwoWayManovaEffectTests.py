'''
Test for effects in two-way MANOVA.
'''
from collections import namedtuple
from IPython.display import display, Math
import numpy as np
from scipy import stats
import scipy.linalg as la
from TwoWayManova import TwoWayManova


class TwoWayManovaEffectTests(object):
    '''
    A class to compute the tests for interaction and possibly factor effects in two-way MANOVA.
    '''
    def __init__(self, manova: TwoWayManova, alpha: float):
        self._manova = manova
        self.alpha = alpha

        self.g = manova._design.metadata.g
        self.b = manova._design.metadata.b
        self.p = manova._design.metadata.p
        self.n = manova._design.metadata.n

        self.df_values = self._compute_df()
        self.lambda_star = self._compute_lambda_star()
        self.f_values = self._compute_f_test_values()
        self.f_critical_values = self._compute_f_critical()

    def _compute_lambda_star(self) -> namedtuple:
        r'''
        Compute the $\Lambda^{\star}$ values used in the F-tests.
        Returns:
            namedtuple:
        '''
        lambda_star_interaction = la.det(self._manova.manova.W)/(la.det(self._manova.manova.I + self._manova.manova.W))
        lambda_star_factor1 = la.det(self._manova.manova.W) / la.det(self._manova.manova.B1 + self._manova.manova.W)
        lambda_star_factor2 = la.det(self._manova.manova.W) / la.det(self._manova.manova.B2 + self._manova.manova.W)
        LambdaStarResults = namedtuple('LambdaStarResults', ['SumOfSquaresMatrix', 'LambdaStar'])
        LambdaStar = namedtuple('LambdaStar', ['Interaction', 'Factor1', 'Factor2'])
        return LambdaStar(Interaction=LambdaStarResults(SumOfSquaresMatrix=self._manova.manova.I, LambdaStar=lambda_star_interaction),
                          Factor1=LambdaStarResults(SumOfSquaresMatrix=self._manova.manova.B1, LambdaStar=lambda_star_factor1),
                          Factor2=LambdaStarResults(SumOfSquaresMatrix=self._manova.manova.B2, LambdaStar=lambda_star_factor2))

    def _compute_f_test_values(self) -> namedtuple:
        '''
        Compute the test values for the F-tests based on the data.
        Returns:
            namedtuple:
        '''
        g = self.g
        b = self.b
        p = self.p
        n = self.n

        f_interaction = ((1-self.lambda_star.Interaction.LambdaStar)/self.lambda_star.Interaction.LambdaStar)*((g*b*(n-1)-p+1)/2) / ((np.abs((g-1)*(b-1) - p) + 1)/2)
        f_factor1 = ((1-self.lambda_star.Factor1.LambdaStar)/self.lambda_star.Factor1.LambdaStar)*((g*b*(n-1)-p+1)/2) / ((np.abs((g-1) - p) + 1)/2)
        f_factor2 = ((1-self.lambda_star.Factor2.LambdaStar)/self.lambda_star.Factor2.LambdaStar)*((g*b*(n-1)-p+1)/2) / ((np.abs((b-1) - p) + 1)/2)
        FValues = namedtuple('FValues', ['Interaction', 'Factor1', 'Factor2'])
        return FValues(Interaction=f_interaction,
                       Factor1=f_factor1,
                       Factor2=f_factor2)
    
    def _compute_df(self) -> namedtuple:
        '''
        Compute the numerator and denominator degrees of freedom for the F-tests.
        Returns:
            namedtuple:
        '''
        g = self.g
        b = self.b
        p = self.p
        n = self.n

        v1_int = np.abs((g-1)*(b-1) - p) + 1
        v1_int_text = fr'\nu_{{1}} = \left| (g-1)(b-1) - p \right| + 1 = \left| {g-1}({b-1}) - {p} \right| + 1 = {v1_int}'

        v1_fac1 = np.abs((g-1) - p) + 1
        v1_fac1_text = fr'\nu_{1} = |(g-1) - p| + 1 = |{g-1} - {p}| + 1 = {v1_fac1}'
        
        v1_fac2 = np.abs((b-1) - p) + 1
        v2_fac2_text = fr'\nu_{1} = |(b-1) - p| + 1 = |{b-1} - {p}| + 1 = {v1_fac2}'

        v2 = g*b*(n-1) - p + 1
        v2_text = fr'\nu_{{2}} = gb(n - 1) - p + 1 = {g}({b})({n-1}) - {p} + 1 = {v2}'

        DegreesOfFreedom = namedtuple('DegreesOfFreedom', ['v1', 'v1Text', 'v2', 'v2Text'])
        TwoWayManovaDF = namedtuple('TwoWayManovaDF', ['Interaction', 'Factor1', 'Factor2'])
        return TwoWayManovaDF(Interaction=DegreesOfFreedom(v1=v1_int, v1Text=v1_int_text, v2=v2, v2Text=v2_text),
                              Factor1=DegreesOfFreedom(v1=v1_fac1, v1Text=v1_fac1_text, v2=v2, v2Text=v2_text),
                              Factor2=DegreesOfFreedom(v1=v1_fac2, v1Text=v2_fac2_text, v2=v2, v2Text=v2_text))

    def _compute_f_critical(self) -> namedtuple:
        '''
        Compute the critical values using the F-distribution to compare the test value against.
        Returns:
            namedtuple:
        '''
        f_crit_interaction = stats.f.ppf(1-self.alpha, self.df_values.Interaction.v1, self.df_values.Interaction.v2)
        f_crit_factor1 = stats.f.ppf(1-self.alpha, self.df_values.Factor1.v1, self.df_values.Factor1.v2)
        f_crit_factor2 = stats.f.ppf(1-self.alpha, self.df_values.Factor2.v1, self.df_values.Factor2.v2)
        FCritical = namedtuple('FCritical', ['Interaction', 'Factor1', 'Factor2'])
        return FCritical(Interaction=f_crit_interaction,
                         Factor1=f_crit_factor1,
                         Factor2=f_crit_factor2)
    
    def build_latex_parameters(self) -> namedtuple:
        '''
        Build the string that describes the null hypothesis for the different tests.
        Returns:
            namedtuple:
        '''
        interaction_params = []
        added_text = r'\textbf{0}'
        for i in range(self.g):
            for j in range(self.b):
                interaction_params.append(fr' \bm{{ \gamma }}_{{ {i},{j} }} ')
        interaction_params.append(added_text)

        # interaction_params = '='.join(interaction_params_list)
        factor1_params = [fr' \bm{{ \tau }}_{{ {i+1} }} ' for i in range(self.g)]
        factor1_params.append(added_text)

        factor2_params = [fr' \bm{{ \beta }}_{{ {i+1} }} ' for i in range(self.b)]
        factor2_params.append(added_text)

        ParameterStatement = namedtuple('ParameterStatement', ['Interaction', 'Factor1', 'Factor2'])
        return ParameterStatement(Interaction='='.join(interaction_params),
                                  Factor1='='.join(factor1_params),
                                  Factor2='='.join(factor2_params))

    def display_test_results(self) -> None:
        '''
        Display the computations and results of tests for effects.
        '''
        g = self.g
        b = self.b
        p = self.p
        n = self.n
        parameter_text = self.build_latex_parameters()
        groups = {'Interaction': 'int', 'Factor1': 'fac1', 'Factor2': 'fac2'}
        test_index = {'Interaction': '', 'Factor1': '1', 'Factor2': '2'}
        for group, data in self.lambda_star._asdict().items():
            label = groups.get(group)
            index = test_index.get(group)
            df = self.df_values._asdict().get(group)
            f_test = self.f_values._asdict().get(group)
            f_crit = self.f_critical_values._asdict().get(group)
            params = parameter_text._asdict().get(group)
            
            display(Math(fr'\textrm{{Test for {group}: }}'))
            display(Math(r'\qquad\qquad\Lambda_{1}^{\star}'
                 '='
                 fr'\frac{{\left| \text{{SSP}}_{{\text{{res}}}} \right|}}{{\left| \text{{SSP}}_{{\text{{ {label} }}}} + \text{{SSP}}_{{\text{{res}}}} \right|}}'
                 '='
                 fr'\frac{{ {la.det(self._manova.manova.W):.4f} }}{{ {la.det(data.SumOfSquaresMatrix +self._manova.manova.W):.4f} }}'
                 '='
                 f'{data.LambdaStar:.4f}'
                 ))
            display(Math(r'\qquad\qquad' + df.v1Text))
            display(Math(fr'\qquad\qquad\nu_{2} = gb(n - 1) - p + 1 = {g}({b})({n-1}) - {p} + 1 = {df.v2}'))
            display(Math(fr'\qquad\qquad F_{{ {index} }}'
                         '='
                         r'\left( \frac{1 - \Lambda^{\star}}{\Lambda^{\star}} \right)'
                         fr'\frac{{ v_{{2}} /2 }}{{ v_{{1}} /2 }}'
                         '='
                         fr'\left( \frac{{ 1 - {data.LambdaStar:.4f} }}{{ {data.LambdaStar:.4f} }} \right)'
                         fr'\frac{{ {df.v2/2} }}{{ {df.v1/2} }}'
                         '='
                         f'{f_test:.2f}'
                         ))
            display(Math(fr'\qquad\qquad F_{{\nu_{{1}}, \nu_{{2}}}}(\alpha) = F_{{ {df.v1}, {df.v2} }}({self.alpha}) = {f_crit:.2f}'))
            
            if f_test > f_crit:
                display(Math(fr'\text{{We have that }} F_{{ {index} }} = {f_test:.2f} > F_{{ \text{{crit}} }} = F_{{ {df.v1}, {df.v2} }} \left( {self.alpha} \right) = '
                             fr'{f_crit:.2f} \text{{, so we would reject the null hypothesis that }} '
                             f'{params}'
                             ))
            else:
                display(Math(fr'\text{{We have that }} F_{{ {index} }} = {f_test:.2f} < F_{{ \text{{crit}} }} = F_{{ {df.v1}, {df.v2} }} \left( {self.alpha} \right) = '
                            fr'{f_crit:.2f} \text{{, so we would fail to reject the null hypothesis that }} '
                            f'{params}'
                            ))
            if group == 'Interaction':
                if f_test > f_crit:
                    display(Math(r'\textrm{The test for interaction was significant, so testing for individual factor effects makes no sense.}'))
                    break
