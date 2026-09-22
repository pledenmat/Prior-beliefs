# Manipulating Prior Beliefs Causally Induces Under- and Overconfidence

Replication materials for Van Marcke, Le Denmat, Verguts, & Desender (2024),
*Psychological Science*, 35(4), 358-375. https://doi.org/10.1177/09567976241231572

## Two experiments

- **Exp. 1** ("fake feedback"): prior beliefs manipulated via fake
  comparative feedback during training.
- **Exp. 2** ("training difficulty"): prior beliefs manipulated via the
  objective difficulty of training trials.

## Pipeline / run order

There are two largely independent analysis paths that both start from the
same raw data (`Data/Data_1/` = Exp. 1, `Data/Data_2/` = Exp. 2, one CSV per
participant).

**1. Behavioral / mixed-model statistics**:
- `BehaviouralDataAnalysis_Exp1.R`
- `BehaviouralDataAnalysis_Exp2.R`

Each script is self-contained: it loads its own raw data, cleans it, and
runs the linear mixed-effects models reported in the paper.

**2. Drift-diffusion / quantitative modeling** (heat-map based confidence
model, Figures 3-4, Table 2):
1. `1_preprocessing.R` - loads the raw per-participant CSVs, cleans them,
   writes `Data/Aggregated/*.csv`.
2. `2_model_fit_and_analysis.R` - sources `1_preprocessing.R`, then
   `quantile_fit_DDM.R`, fits the DDM models. This is the pipeline used for
   the paper's qualitative DDM results (Fig. 4-style parameters).
3. `quantitative_fit.R` - sources `1_preprocessing.R` and
   `quantile_optim_DDM_Vs_bias.R`; combines precomputed per-participant fits
   (under `Fits/Exp{1,2}/quantitative_fit/`) into simulated datasets and
   writes `quantitative_prediction_exp{1,2}.csv`, consumed by
   `Fig3and4.Rmd` to produce Figures 3-4.

Both paths need the raw data under `Data/` and the precomputed fits under
`Fits/`. Model fitting itself (`2_model_fit_and_analysis.R`'s DEoptim calls)
is not re-run by `quantitative_fit.R` - it reads the `.Rdata` results that
are already checked in.

## Dependencies

Standard CRAN packages: reshape, effects, lmerTest, scales, DEoptim, car,
MALDIquant, Rcpp, RcppZiggurat (+ RcppGSL, needs system `libgsl-dev`),
readxl, multcomp, lme4, afex, ggplot2, ggthemes, ggridges, matrixStats,
phia, quickpsy, BayesFactor, fields, mnormt, zoo, emmeans, timeSeries.

## What was fixed (Sept 2026)

An independent reproducibility report flagged several issues with this
project as it was archived on OSF (https://osf.io/8bf3r/). Fixed here:

1. **This README** - none existed before.
2. **`postchecks` referencing non-existent columns in
   `BehaviouralDataAnalysis_Exp1.R`** - `SubjectiveInfluenceFb`,
   `FeedbackCredibility`, `ManipulationAwareness` are column names from a previously missing file 
   "RecodedPostchecks_Exp1.csv". This file has now been added to the repository and is referenced in the analyses.
3. **`quantile_optim_DDM_Vs_bias.R` hardcoded a `coh` column** that has
   never existed in this project's data; `1_preprocessing.R` names it
   `trialdifflevel`. All internal `$coh` references now say `$trialdifflevel`.
4. **`quantitative_fit.R` called `quantile_optim_DDM_Vs_biasfixed(...,
   condition_name = "selfconf")`** for Exp. 1, but `1_preprocessing.R` names
   that column `fbcond`. Fixed to `"fbcond"`.
5. **`library(myPackage)`** - a personal package never included in the
   replication materials. The only function actually needed from it,
   `fast_hm()`, is now sourced directly from `fast_hm.R` to remove the dependency on the package.
6. **`library(prob)`** - archived from CRAN on 2022-04-29. Turns out this package was never actually used anywhere in this pipeline, so we removed the call to this package.
7. **`BehaviouralDataAnalysis_Exp1.R`**  specified a  3x27 post-hoc contrast matrix C which didn't
    match `fit3` as it was leftover code from a follow-up analysis.The analysis now specifies the correct code for
    post hoc tests of fit3, and the C-matrix is cut.
8. **Undefined `fit11` / `fit9`** - not flagged by the original report, but
   found while fixing the above: both behavioral scripts have a whole
   post-hoc section (summary/Anova/plots/glht/testInteractions) written
   against a model variable that's never defined (`fit11` in Exp1, `fit9` in
   Exp2). Repointed at the actual winning model (`fit3` in Exp1, per the
   comment directly above it; `fit1` in Exp2, the model being summarized
   immediately before).
9. **`BehaviouralDataAnalysis_Exp{1,2}.R`** specified Anova(type=3) for the analysis of confidence, 
    which produces chi square statistics. Upon request of the editor, we reported F-values instead, 
    which are obtained by anova() as now updated in the code.
9. **`quantitative_fit.R` wrongly loaded estimated parameters** - Drift rates were assigned per condition instead of per difficulty level due to pre-publication variable naming changes, leading to reproduced plots in Figure 3 to be way off. Now fixed.


