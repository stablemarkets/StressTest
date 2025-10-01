# Bayesian Sensitivity Analyses in Causal Inference

<p align="center">
  <img src="/mnar_outcome/bnp_example.png" width="500" />
</p>

---

This is the companion repository for the paper "Stress-Testing Assumptions: A Practical Guide to Bayesian Sensitivity Analyses in Causal Inference""

### Directory

---
This repository contains the following sub-directories:
- `gcomp_example`: contains code implementing model run in Section 2.
- `misclassification_example`: contains code implementing model discussed in Section 3.1. Generates Figure 2 (left).
- `unmeasured_confounding`: contains code implementing model discussed in Section 3.2. Generates Figure 2 (right).
- `mnar_outcome`: contains code estimating model discussed in Section 3.3/4. Generates Figures 3 and 4.

Each sub-directory contains a .R file and a .stan file, as per the workflow described in the paper.