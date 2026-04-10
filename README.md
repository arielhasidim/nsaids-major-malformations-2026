# First Trimester Nonsteroidal Anti-inflammatory Drugs Exposure and Risk of Major Congenital Malformations: a retrospective register-based cohort study

[![DOI](https://zenodo.org/badge/1182896951.svg)](https://doi.org/10.5281/zenodo.19043979)

## Abstract

**Background:** Pain and fever are common in early pregnancy, yet their management poses a major clinical dilemma. Although not confirmed, recent studies have raised safety concerns regarding acetaminophen. Evidence on the use of nonsteroidal anti-inflammatory drugs (NSAID) in the first trimester remains inconclusive. This uncertainty has left clinicians with limited evidence to guide treatment decisions. This study evaluated the association between first-trimester NSAID exposure and the risk of major congenital malformations (MCMs) in a large, population-based cohort of pregnancies.

**Methods and Findings:** We conducted a population-based retrospective cohort study within the Southern Israeli Pregnancy Registry (siPREG) project, including all singleton pregnancies of women aged 15–45 years resulting in live births, stillbirths, or elective terminations for fetal malformations at a Soroka University Medical Center between 1998 and 2018. Pregnancies exposed to established teratogens, multiple gestations, and those with documented genetic or chromosomal anomalies were excluded. First-trimester NSAID exposure was defined by pharmacy dispensations (overall and by specific agents). MCMs were identified from linked clinical, hospitalization, and termination records through the first postnatal year.

Propensity scores were estimated using covariates selected via a directed acyclic graph, including maternal age, ethnicity, diabetes, medical indication for NSAID use, exposure to other antipyretics, obesity, smoking, folic-acid use, gravidity, perinatal care, and year of pregnancy. Generalized full matching was used to balance covariates. Adjusted risk ratios were derived using weighted Poisson-regression with G-computation, and two-way cluster-robust standard errors, jointly clustering by maternal identifier and matching subclass. Sensitivity analyses included a dose–response assessment across defined-daily-dose (DDD) categories and a tipping-point analysis evaluating the impact of potential misclassification from unrecorded over-the-counter NSAID use.

A total of 264,858 singleton pregnancies were included in the final cohort; 20,202 (7.6%) were exposed to NSAID, most commonly ibuprofen (5.1%), diclofenac (1.6%), and naproxen (1.2%). NSAID exposure, in total and as individual agents, was not associated with MCMs overall (8.2% vs. 7.0%; matched-adjusted-Relative-Risk (aRR) = 0.99 (95% CI [0.90,1.10])) or with organ-system–specific MCMs, including cardiovascular (matched-aRR=1.05 (95%CI [0.92,1.20]), musculoskeletal (matched-aRR =1.03 (95%CI [0.77,1.39])), central nervous system (matched-aRR=0.77 (95%CI [0.53,1.11])), cleft palate (matched-aRR=0.95 (95%CI [0.47–1.91])), gastrointestinal (matched-aRR=1.03 (95%CI [0.64–1.63])), and genitourinary (matched-aRR=0.99 (95%CI [0.72,1.35])) malformations. Dose–response analyses showed no significant association with MCMs across cumulative NSAID exposure: short-term (1–7 DDD, matched-aRR = 1.06 (95% CI [0.97,1.15]), medium-term (8–21 DDD, matched-aRR = 1.10 (95% CI [0.99,1.22]), and long-term (>21 DDD, matched-aRR = 1.24 (95% CI [0.94,1.63])). The main limitation was the potential for minor exposure misclassification due to over-the-counter availability of ibuprofen, although sensitivity analyses simulating such misclassification suggested minimal impact on the risk estimates.


**Conclusion:** In this large, population-based cohort, we found no evidence supporting an association between first-trimester exposure to NSAID and MCMs, providing reassuring evidence regarding their fetal safety in early pregnancy.

## Additional data files

This repository also includes two CSV files containing the numerical values underlying selected study figures, provided to support data transparency.

- **`annual_exposure_rates.csv`**: numerical data underlying **Figure 3**. Each row represents a calendar year (`BirthTerminationYear`), and each drug-specific column contains the **absolute number of pregnancies** with first-trimester dispensation exposure to that NSAID in that year. Units: **count of exposed pregnancies per year**.

- **`UpSet.csv`**: numerical data underlying **Supplementary Figure S1.1** (UpSet plot of overlapping antipyretic/NSAID exposures). Each row represents one exposure combination (`combination`), `freq` gives the **absolute number of pregnancies** in that combination, and the binary drug columns indicate whether each medication/exposure category is included in the combination (`1` = included, `0` = not included). `n_sets` indicates the **number of exposure categories** included in that combination. Units: **count of pregnancies** for `freq`; **binary indicator (0/1)** for exposure-membership columns.

`Unnamed: 0` is a row index exported from R and does not represent a study variable.

## Citation

If you use this codebase or build upon it in academic work, please cite both the associated article and the archived software record.

### Associated article

Hasidim AA, Ben Shitrit I, Idan D, Michael T, Levy A, Pariente G, Lunenfeld E, Daniel S. *First Trimester Nonsteroidal Anti-inflammatory Drugs Exposure and Risk of Major Congenital Malformations: a retrospective register-based cohort study*. *PLOS Medicine*. 2026. https://doi.org/10.1371/journal.pmed.1005063

### Software archive

Hasidim AA. *nsaids-major-malformations-2026*. Zenodo. 2026. https://doi.org/10.5281/zenodo.19043979

```bibtex
@article{hasidim2026nsaids,
  author  = {Hasidim, Ariel Avraham and Ben Shitrit, Itamar and Idan, Daphna and Michael, Tal and Levy, Amalia and Pariente, Gali and Lunenfeld, Eitan and Daniel, Sharon},
  title   = {First Trimester Nonsteroidal Anti-inflammatory Drugs Exposure and Risk of Major Congenital Malformations: a retrospective register-based cohort study},
  journal = {PLOS Medicine},
  year    = {2026},
  doi     = {10.1371/journal.pmed.1005063}
}

@software{hasidim2026software,
  author    = {Hasidim, Ariel Avraham},
  title     = {nsaids-major-malformations-2026},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.19043979},
  url       = {https://doi.org/10.5281/zenodo.19043979}
}
