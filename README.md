# Bayesian-vs-Synthetic

## Description

This contains R code to test a Bayesian dynamic borrowing (BDB) and synthetic control method (SCM) on a case study of pediatric atopic dermatitis. Six historical randomised control trials looking at treatments for pediatric atopic dermatitis were selected to use as the base data.  The same studies were used in both methods. The synthetic controls were generated using the CART method from the synthpop library. The RBest library was used to test the BDB method. Control response rate, type 1 error rate, and power are used as comparison metrics.


## Installation

The script can be downloaded and run in R Studio.

## Libraries

The R Libraries used include the following:

- synthpop
- dplyr
- ggplot2
- RBest

## Usage

This code is for comparing the quality of BDB and SCM in a case study. The study information can be changed and fit to other diseases of interest.

## Citations

### More information on Synthpop Synthesis Methods

- Nowok, B., Raab, G.M. and Dibben, C. (2016) ‘synthpop: Bespoke Creation of Synthetic Data in R’, Journal of Statistical Software, 74, pp. 1–26. Available at: https://doi.org/10.18637/jss.v074.i11.

### More information on RBest Bayesian Methods

- Weber, S. et al. (2021) ‘Applying Meta-Analytic-Predictive Priors with the R Bayesian Evidence Synthesis Tools’, Journal of Statistical Software, 100, pp. 1–32. Available at: https://doi.org/10.18637/jss.v100.i19.


### Historical Randomised Control Trials

- Simpson, E.L. et al. (2020) ‘Efficacy and Safety of Dupilumab in Adolescents With Uncontrolled Moderate to Severe Atopic Dermatitis: A Phase 3 Randomized Clinical Trial’, JAMA Dermatology, 156(1), pp. 44–56. Available at: https://doi.org/10.1001/jamadermatol.2019.3336.
- Paller, A.S. et al. (2022) ‘Dupilumab in children aged 6 months to younger than 6 years with uncontrolled atopic dermatitis: a randomised, double-blind, placebo-controlled, phase 3 trial’, The Lancet, 400(10356), pp. 908–919. Available at: https://doi.org/10.1016/S0140-6736(22)01539-2.
- Paller, A.S. et al. (2020) ‘Efficacy and safety of dupilumab with concomitant topical corticosteroids in children 6 to 11 years old with severe atopic dermatitis: A randomized, double-blinded, placebo-controlled phase 3 trial’, Journal of the American Academy of Dermatology, 83(5), pp. 1282–1293. Available at: https://doi.org/10.1016/j.jaad.2020.06.054.
- Ebisawa, M. et al. (2024) ‘Efficacy and safety of dupilumab with concomitant topical corticosteroids in Japanese pediatric patients with moderate-to-severe atopic dermatitis: A randomized, double-blind, placebo-controlled phase 3 study’, Allergology International, 73(4), pp. 532–542. Available at: https://doi.org/10.1016/j.alit.2024.04.006.
- Torrelo, A. et al. (2023) ‘Efficacy and safety of baricitinib in combination with topical corticosteroids in paediatric patients with moderate-to-severe atopic dermatitis with an inadequate response to topical corticosteroids: results from a phase III, randomized, double-blind, placebo-controlled study (BREEZE-AD PEDS)’, British Journal of Dermatology, 189(1), pp. 23–32. Available at: https://doi.org/10.1093/bjd/ljad096.
- Paller, A.S. et al. (2023) ‘Efficacy and Safety of Tralokinumab in Adolescents With Moderate to Severe Atopic Dermatitis: The Phase 3 ECZTRA 6 Randomized Clinical Trial’, JAMA Dermatology, 159(6), pp. 596–605. Available at: https://doi.org/10.1001/jamadermatol.2023.0627.


