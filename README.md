# Manipulating Prior Beliefs Causally Induces Under- and Overconfidence

Replication materials for Van Marcke, Denmat, Verguts, & Desender (2024),
*Psychological Science*, 35(4), 358-375. https://doi.org/10.1177/09567976241231572

This README was added while fixing the issues raised in an independent
computational-reproducibility report (Andreoni, González-Sahagún, Kozyra,
Kozyra, & Fei, Sept 2026). See "Known limitations" below for the two things
that report flagged which we deliberately left unresolved rather than guess at.

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
   `Fig3_4.Rmd` to produce Figures 3-4.

Both paths need the raw data under `Data/` and the precomputed fits under
`Fits/`. Model fitting itself (`2_model_fit_and_analysis.R`'s DEoptim calls)
is not re-run by `quantitative_fit.R` - it reads the `.Rdata` results that
are already checked in.

## Dependencies

Standard CRAN packages: reshape, effects, lmerTest, scales, DEoptim, car,
MALDIquant, Rcpp, RcppZiggurat (+ RcppGSL, needs system `libgsl-dev`),
readxl, multcomp, lme4, afex, ggplot2, ggthemes, ggridges, matrixStats,
phia, quickpsy, BayesFactor, fields, mnormt, zoo, emmeans, timeSeries.

Two things that used to be hard, unresolvable dependencies have been fixed
(see "What was fixed" below): the archived `prob` package and the missing
custom `myPackage`.

## What was fixed (Sept 2026)

An independent reproducibility report flagged several issues with this
project as it was archived on OSF (https://osf.io/8bf3r/). Fixed here:

1. **This README** - none existed before.
2. **`postchecks` undefined in `BehaviouralDataAnalysis_Exp2.R`** - it was
   referenced but never created; added the missing definition, mirroring
   Exp1's own construction of the same object.
3. **`postchecks` referencing non-existent columns in
   `BehaviouralDataAnalysis_Exp1.R`** - `SubjectiveInfluenceFb`,
   `FeedbackCredibility`, `ManipulationAwareness` were never actual column
   names (see "Known limitations" - this one is only guarded, not fixed).
4. **`quantile_optim_DDM_Vs_bias.R` hardcoded a `coh` column** that has
   never existed in this project's data; `1_preprocessing.R` names it
   `trialdifflevel`. All internal `$coh` references now say `$trialdifflevel`.
5. **`quantitative_fit.R` called `quantile_optim_DDM_Vs_biasfixed(...,
   condition_name = "selfconf")`** for Exp. 1, but `1_preprocessing.R` names
   that column `fbcond`. Fixed to `"fbcond"`.
6. **`library(myPackage)`** - a personal package never included in the
   replication materials. The only function actually needed from it,
   `fast_hm()`, is now sourced directly from `fast_hm.R` to remove the dependency on the package.
   sourced directly instead.
7. **`library(prob)`** - archived from CRAN on 2022-04-29. Turns out this package was never actually used anywhere in this pipeline, so we removed the call to this package.
8. **Undefined `fit11` / `fit9`** - not flagged by the original report, but
   found while fixing the above: both behavioral scripts have a whole
   post-hoc section (summary/Anova/plots/glht/testInteractions) written
   against a model variable that's never defined (`fit11` in Exp1, `fit9` in
   Exp2). Repointed at the actual winning model (`fit3` in Exp1, per the
   comment directly above it; `fit1` in Exp2, the model being summarized
   immediately before).


## Known limitations (not fixed - need your input)

**Postcheck coding is missing, not just miscoded.** `post1`...`post6` in the
raw data are free-text debriefing answers (in Dutch). The three named
variables the analysis expects (`SubjectiveInfluenceFb`,
`FeedbackCredibility`, `ManipulationAwareness`) look like they were meant to
come from a manually-coded file (`library(readxl)` is loaded at the top of
both behavioral scripts but never actually used to read anything) that
isn't in either the OSF or GitHub materials. We did not attempt to
auto-code the Dutch free text into these categories ourselves - guessing
wrong would silently corrupt this part of the analysis rather than just
fail. The postcheck follow-up analysis (not the main confidence result) is
guarded to skip cleanly with a warning until that file is supplied.

**A post-hoc contrast matrix in `BehaviouralDataAnalysis_Exp1.R` doesn't
match `fit3`.** Right after the `fit11`->`fit3` fix, there's a hand-built
3x27 contrast matrix `C` for `glht()`, but `fit3` only has 9 fixed-effect
coefficients. 27 is exactly what a three-way interaction with another
3-level factor would need - almost certainly one of the `fitFC`/`fitMA`/
`fitSI` postcheck models further down, not `fit3`. Without knowing which
specific cells the original four `C[i,i] <- 1` entries were meant to test,
we didn't guess a replacement; it's wrapped in `tryCatch()` so the rest of
the script keeps running, with `testInteractions()` (which doesn't need a
hand-built matrix) left in place as the working alternative.

**One raw file has a formatting quirk.**
`Data/Data_1/selfconfidence1A_sub29.csv` throws an "EOF within quoted
string" warning from `read.csv` (likely an unescaped quote in a free-text
answer) - pre-existing in the original data, not introduced by any fix
here. It's a warning, not an error, and doesn't stop preprocessing.
