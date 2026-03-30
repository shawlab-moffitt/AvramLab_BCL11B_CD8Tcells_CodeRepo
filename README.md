
# 📊 Tumor Growth Analysis -- BCL11B PMEL/B16 Model

This repository contains code for analyzing longitudinal tumor growth in
a **B16 melanoma model with PMEL adoptive T cell transfer**, designed to
evaluate the impact of **BCL11B perturbation** on tumor progression and
immune-mediated tumor control.

------------------------------------------------------------------------

## 🔬 Overview

BCL11B is a key transcription factor that regulates **T cell lineage
commitment, differentiation, and functional state**. In this study, we
assess how perturbation of BCL11B alters **tumor growth kinetics in
vivo**, using a **longitudinal modeling framework**.

Rather than relying on endpoint tumor size, this analysis quantifies
**tumor growth trajectories over time**, enabling detection of subtle
but biologically meaningful differences in tumor control.

------------------------------------------------------------------------

## 📁 Repository Structure

    .
    ├── TumorGrowth_20260324.R
    ├── data/
    │   └── PMEL_TumorGrowthData.tsv
    ├── output/
    └── README.md

------------------------------------------------------------------------

## ⚙️ Requirements

``` r
install.packages(c("lme4", "lmerTest", "emmeans", "ggplot2"))
```

------------------------------------------------------------------------

## 📊 Data Format

  Column      Description
  ----------- ----------------------------------
  ID          Unique identifier for each mouse
  Days        Time point
  TumorSize   Tumor measurement
  Group       WT, NT, KO

------------------------------------------------------------------------

## 🚀 Analysis Workflow

### Model

``` r
log2(TumorSize + 1) ~ Days + Group + Days:Group + (1 | ID)
```

-   Fixed: Days, Group, Interaction\
-   Random: Mouse (ID)

------------------------------------------------------------------------

## 📈 Interpretation

-   Interaction term → difference in tumor growth rates\
-   Mixed model → handles repeated measures\
-   Output → slope differences + P-values

------------------------------------------------------------------------

## 👨‍🔬 Authors

Shaw Lab, Moffitt Cancer Center
