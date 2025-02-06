'''
Tools to compute two-way MANOVA and display the MANOVA table.
'''
from collections import namedtuple
from dataclasses import dataclass, field
from IPython.display import display, Math
from scipy import stats
import numpy as np
import pandas as pd
from chapter_6_utils import create_array_text, compute_manova_ss_matrices

# Used to store metadata for the two-way MANOVA.
Metadata = namedtuple('Metadata', ['g', 'b', 'p', 'n'])

@dataclass
class TwoWayDesignMatrix:
    df: pd.DataFrame
    factor1_name: str
    factor2_name: str
    variable_names: list[str]
    metadata: Metadata = field(init=False)

    def __post_init__(self):
        self.metadata = self._compute_metadata()
        self.df = self._subset_df()

    def _subset_df(self) -> pd.DataFrame:
        keep_cols =  [self.factor1_name, self.factor2_name] + self.variable_names
        df = self.df[keep_cols].copy()
        df['RepCol'] = df.groupby(['Factor1', 'Factor2']).cumcount()
        return df

    def _compute_metadata(self) -> Metadata:

        return Metadata(g = self.df[self.factor1_name].nunique(),
                        b = self.df[self.factor2_name].nunique(),
                        p=len(self.variable_names),
                        n = self._determine_reps()
                        )
    
    def _determine_reps(self) -> int:
        rep_counts = self.df.value_counts([self.factor1_name, self.factor2_name])
        if rep_counts.eq(1).all():
            return 1
        else:
            assert rep_counts.all() & rep_counts.gt(1).all(), 'Data is unbalanced!'
            return rep_counts.unique().item()

class ObsBreakdownTwoWayManova(object):
    def __init__(self, design: Metadata) -> namedtuple:
        self._design = design

        self.means = self.comp_means()
        self.anova_values = self.comp_two_way_anova()
        self.obs_breakdown = self.comp_obs_breakdown()

    def comp_means(self) -> namedtuple:
        MeansTwoWayANOVA = namedtuple('MeansTwoWayANOVA',
                                      ['GlobalMean',
                                       'Trt1Mean',
                                       'Trt2Mean',
                                       'Trt1Trt2Mean'])

        xbars = self._design.df[self._design.variable_names].mean()
        global_xbar_a = np.stack([ xbar * np.ones([self._design.metadata.g, self._design.metadata.b])
                                  for xbar in xbars])

        xbar_ell = self._design.df.groupby(self._design.factor1_name)[self._design.variable_names].mean().to_numpy()
        xbar_ell_a = np.stack([np.repeat(xbar_ell[:, i:i+1], self._design.metadata.b, axis=1)
                               for i in range(self._design.metadata.p)])

        xbar_k = self._design.df.groupby(self._design.factor2_name)[self._design.variable_names].mean().to_numpy()
        xbar_k_a = np.stack([np.repeat(xbar_k[:, i:i+1].T, self._design.metadata.g, axis=0)
                             for i in range(self._design.metadata.p)])

        xbar_ell_k_a = self._design.df.groupby([self._design.factor1_name, self._design.factor2_name])[self._design.variable_names]\
            .mean().to_numpy().T\
            .reshape((self._design.metadata.p, self._design.metadata.g, self._design.metadata.b))
        
        return MeansTwoWayANOVA(GlobalMean=global_xbar_a,
                                Trt1Mean=xbar_ell_a,
                                Trt2Mean=xbar_k_a,
                                Trt1Trt2Mean=xbar_ell_k_a
                                )

    def comp_two_way_anova(self) -> namedtuple:
        ANOVAPieces = namedtuple('ANOVAPieces', ['Mean',
                                                 'Trt1Effect',
                                                 'Trt2Effect',
                                                 'Interaction'])
        a = self.means
        trt1_effect_a = a.Trt1Mean - a.GlobalMean
        trt2_effect_a = a.Trt2Mean - a.GlobalMean
        interaction_a = a.Trt1Trt2Mean - a.Trt1Mean - a.Trt2Mean + a.GlobalMean

        return ANOVAPieces(Mean=a.GlobalMean,
                           Trt1Effect=trt1_effect_a,
                           Trt2Effect=trt2_effect_a,
                           Interaction=interaction_a)

    def obs_single_rep(self) -> namedtuple:
        '''
        Breakdown observations for an input variable into mean, treatment, maybe interaction (if there are replications), and residual components.
        Args:
            df (pd.DataFrame): Input data with a column for treatments and columns for each variable.
            trt1_col (str): The column with the treatment one (groups).
            trt2_col (str): The column with the treatment two (groups).
            var_col (str): The variable breakdown. Only one!
        Return:
            namedtuple: The named tuple has elements: 'Variable', 'Obs', 'Mean', 'TreatmentEffect1', 'TreatmentEffect2', 'Residual'.
        '''
        # Store the output for the observations breakdown for a given variable.
        tuple_values = ['Variable',
                        'Obs',
                        'Mean',
                        'TreatmentEffect1',
                        'TreatmentEffect2',
                        'Residual']
        ObsBreakdown = namedtuple('ObsBreakdown', tuple_values)

        obs_a = pd.pivot(self._design.df,
                         index=self._design.factor1_name,
                         columns=self._design.factor2_name,
                         values=self._design.variable_names).to_numpy()
        return ObsBreakdown(Variable=self._design.variable_names,
                            Obs=obs_a,
                            Mean=self.anova_values.Mean,
                            TreatmentEffect1=self.anova_values.Trt1Effect,
                            TreatmentEffect2=self.anova_values.Trt2Effect,
                            Residual=self.anova_values.Interaction)

    def obs_multi_rep(self) ->list[namedtuple]:
        '''
        Breakdown observations for an input variable into mean, treatment, interaction, and residual components.
        Args:
            df (pd.DataFrame): Input data with a column for treatments and columns for each variable.
            trt1_col (str): The column with the treatment one (groups).
            trt2_col (str): The column with the treatment two (groups).
            var_col (str): The variable breakdown. Only one!
            rep_col (str): Optional parameter for the replication column. Only include if there are replications for treatment combinations.
        Return:
            namedtuple: The named tuple has elements: 'Variable', 'Obs', 'Mean', 'TreatmentEffect1', 'TreatmentEffect2', 'Interaction',
            'Residual'.
        '''
        # Store the output for the observations breakdown for a given variable.
        tuple_values = ['Replication',
                        'Obs',
                        'Mean',
                        'TreatmentEffect1',
                        'TreatmentEffect2',
                        'Interaction',
                        'Residual']

        ObsBreakdown = namedtuple('ObsBreakdown', tuple_values)
        
        obs_list = list()

        for i in self._design.df.RepCol.unique():
            subset_df = self._design.df[self._design.df.RepCol.eq(i)].copy()
            # Need to create an obs ndarray that's of dimension (# of vars, g, b).
            # The depth is the number of variables (measurements). Example, in Exercise 6.14 we have x1 and x2.
            obs_a = np.stack([pd.pivot(subset_df,
                                       index=self._design.factor1_name,
                                       columns=self._design.factor2_name,
                                       values=v).to_numpy() for v in self._design.variable_names])
            
            # Compute the breakdown.
            residual_a = obs_a - self.means.Trt1Trt2Mean
            obs = ObsBreakdown(Replication=i,
                                Obs=obs_a,
                                Mean=self.anova_values.Mean,
                                TreatmentEffect1=self.anova_values.Trt1Effect,
                                TreatmentEffect2=self.anova_values.Trt2Effect,
                                Interaction=self.anova_values.Interaction,
                                Residual=residual_a)
            obs_list.append(obs)
        return obs_list

    def comp_obs_breakdown(self) ->list[namedtuple]:
        '''
        Breakdown observations for an input variable into mean, treatment, and residual components. In this case there are no replications, so the interaction term from the replication case becomes the residual.
        Args:
            df (pd.DataFrame): Input data with a column for treatments and columns for each variable.
            trt1_col (str): The column with the treatment one (groups).
            trt2_col (str): The column with the treatment two (groups).
            var_col (str): The variable breakdown. Only one!
            rep_col (str): Optional parameter for the replication column. Only include if there are replications for treatment combinations.
        Return:
            namedtuple: The named tuple has elements: 'Variable', 'Obs', 'Mean', 'TreatmentEffect1', 'TreatmentEffect2',
            'Residual'.
        '''
        # Store the output for the observations breakdown for a given variable.
        if self._design.metadata.n == 1:
            return self.obs_single_rep()
        else:
            return self.obs_multi_rep()
    
    def display_obs_breakdown(self, spacing: list[str]) -> None:
        '''
        Take a breakdown of observations into components and display them using latex.
        Args:
            data (namedtuple): Contains breakdown of observations into components.
            Row for each group, column for each observation.
            spacing list[str]: A list with 4 string elements.
            Something like, ['0.5cm','2.0cm','2.2cm','2.5cm'].
        '''
        if self._design.metadata.n == 1:
            self.display_obs_breakdown_one_rep(spacing)
        else:
            self.display_obs_breakdown_multi_rep(spacing)
    
    def display_obs_breakdown_multi_rep(self, spacing: list[str]) -> None:
        for obs in self.obs_breakdown:
            # Loop through all variables.
            for i, var in enumerate(self._design.variable_names):
                data = obs
                obs_latex = create_array_text(data.Obs[i])
                mean_latex = create_array_text(data.Mean[i])
                trt_effect1_latex = create_array_text(data.TreatmentEffect1[i])
                trt_effect2_latex = create_array_text(data.TreatmentEffect2[i])
                interaction_latex = create_array_text(data.Interaction[i])
                residual_latex = create_array_text(data.Residual[i])

                assert len(spacing)==6, 'Spacing must have 4 string elements.'
                text_list = [mean_latex, trt_effect1_latex, trt_effect2_latex, interaction_latex, residual_latex]
                latex_str = obs_latex + ' = '
                latex_str += ' + '.join(text_list)

                display(Math(fr'\text{{Variable: {var}}}'))
                display(Math(latex_str))
                display(Math(fr'\hspace{{ {spacing[0]} }}\text{{(observation)}}'
                        fr'\hspace{{ {spacing[1]} }}\text{{(mean)}}'
                        fr'\hspace{{ {spacing[2]} }}\text{{(treatment 1 effect)}}'
                        fr'\hspace{{ {spacing[3]} }}\text{{(treatment 2 effect)}}'
                        fr'\hspace{{ {spacing[4]} }}\text{{(interaction)}}'
                        fr'\hspace{{ {spacing[5]} }}\text{{(residual)}}'))
            
    def display_obs_breakdown_one_rep(self, spacing: list[str]) -> None:
        data = self.obs_breakdown
        obs_latex = create_array_text(data.Obs)
        mean_latex = create_array_text(data.Mean)
        trt_effect1_latex = create_array_text(data.TreatmentEffect1)
        trt_effect2_latex = create_array_text(data.TreatmentEffect2)
        residual_latex = create_array_text(data.Residual)
        
        assert len(spacing)==5, 'Spacing must have 4 string elements.'
        text_list = [mean_latex, trt_effect1_latex, trt_effect2_latex, residual_latex]
        latex_str = obs_latex + ' = '
        latex_str += ' + '.join(text_list)

        display(Math(fr'\text{{Variable: {data.Variable}}}'))
        display(Math(latex_str))
        display(Math(fr'\hspace{{ {spacing[0]} }}\text{{(observation)}}'
                    fr'\hspace{{ {spacing[1]} }}\text{{(mean)}}'
                    fr'\hspace{{ {spacing[2]} }}\text{{(treatment 1 effect)}}'
                    fr'\hspace{{ {spacing[3]} }}\text{{(treatment 2 effect)}}'
                    fr'\hspace{{ {spacing[4]} }}\text{{(residual)}}'))
            
class TwoWayManova(object):
    def __init__(self, df: pd.DataFrame, factor1_col: str, factor2_col: str, var_cols: list[str]):
        self._design = TwoWayDesignMatrix(df=df,
                                          factor1_name=factor1_col,
                                          factor2_name=factor2_col,
                                          variable_names=var_cols)
        self.x_breakdown = ObsBreakdownTwoWayManova(self._design)
        self.manova = self._compute_manova_table()

    def _compute_manova_table(self):
        MANOVA = namedtuple('MANOVA', ['B1','B2','I','W','T'])
        B1, B2, I, W, T = 0, 0, 0, 0, 0
        for obs in self.x_breakdown.obs_breakdown:
            B1 += compute_manova_ss_matrices(obs.TreatmentEffect1)
            B2 += compute_manova_ss_matrices(obs.TreatmentEffect2)
            I  += compute_manova_ss_matrices(obs.Interaction)
            W  += compute_manova_ss_matrices(obs.Residual)
            T  += compute_manova_ss_matrices(obs.Obs) - compute_manova_ss_matrices(obs.Mean)

        return MANOVA(B1=B1,
                      B2=B2,
                      I=I,
                      W=W,
                      T=T)
    
    def display_manova_table(self) -> None:
        '''
        Display the 2-Way MANOVA table.
        '''
        g = self._design.metadata.g
        b = self._design.metadata.b
        n = self._design.metadata.n
        display(Math(r'\begin{array}{lll}'
                r'\text{Source} & \text{Matrix of sum of squares} &  \\'
                r'\text{of variation} & \text{and cross products} & \text{Degrees of freedom} \\'
                r'\hline \\'
                r'\text{Treatment 1} & '
                f'{create_array_text(self.manova.B1)} & '
                fr'{g} - 1 = {g - 1} \\ \\'
                r'\text{Treatment 2} & '
                f'{create_array_text(self.manova.B2)} & '
                fr'{b} - 1 = {b - 1} \\ \\'
                r'\text{Interaction} & '
                f'{create_array_text(self.manova.I)} & '
                fr'({g} - {1})({b} - {1}) = {(g - 1)*(b - 1)} \\ \\'
                r'\text{Residual} & '
                f'{create_array_text(self.manova.W)} &'
                fr'({g})({b})({n} - 1) = {g*b*(n - 1)} \\ \\'
                r'\hline \\'
                r'\text{Total (corrected)} & '
                f'{create_array_text(self.manova.T)} & '
                f'{g}({b})({n}) - {1} = {g*b*n - 1}'
                r'\end{array}'
                ))

