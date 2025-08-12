library(RBesT)


# citations for previous studies:

# Simpson, E.L. et al. (2020) 
# ‘Efficacy and Safety of Dupilumab in Adolescents With Uncontrolled Moderate to Severe Atopic Dermatitis:
# A Phase 3 Randomized Clinical Trial’, JAMA Dermatology, 156(1), pp. 44–56. 
# Available at: https://doi.org/10.1001/jamadermatol.2019.3336.


# Paller, A.S. et al. (2022) 
# ‘Dupilumab in children aged 6 months to younger than 6 years with uncontrolled atopic dermatitis: 
# a randomised, double-blind, placebo-controlled, phase 3 trial’, 
# The Lancet, 400(10356), pp. 908–919.
# Available at: https://doi.org/10.1016/S0140-6736(22)01539-2.


# Paller, A.S. et al. (2020) 
# ‘Efficacy and safety of dupilumab with concomitant topical corticosteroids in children 6 to 11 years 
# old with severe atopic dermatitis: A randomized, double-blinded, placebo-controlled 
# phase 3 trial’, Journal of the American Academy of Dermatology, 83(5), pp. 1282–1293. 
# Available at: https://doi.org/10.1016/j.jaad.2020.06.054.


# Ebisawa, M. et al. (2024) 
# ‘Efficacy and safety of dupilumab with concomitant topical corticosteroids in Japanese pediatric patients 
# with moderate-to-severe atopic dermatitis:
# A randomized, double-blind, placebo-controlled phase 3 study’, 
# Allergology International, 73(4), pp. 532–542. 
# Available at: https://doi.org/10.1016/j.alit.2024.04.006.


# Torrelo, A. et al. (2023) 
# ‘Efficacy and safety of baricitinib in combination with topical corticosteroids in paediatric patients with 
# moderate-to-severe atopic dermatitis with an inadequate response to topical corticosteroids: 
# results from a phase III, randomized, double-blind, placebo-controlled study (BREEZE-AD PEDS)’,
# British Journal of Dermatology, 189(1), pp. 23–32. 
# Available at: https://doi.org/10.1093/bjd/ljad096.


# Paller, A.S. et al. (2023)
# ‘Efficacy and Safety of Tralokinumab in Adolescents With Moderate to Severe Atopic Dermatitis: 
# The Phase 3 ECZTRA 6 Randomized Clinical Trial’, JAMA Dermatology, 159(6),
# pp. 596–605. Available at: https://doi.org/10.1001/jamadermatol.2023.0627.
















# summary level data
data <- data.frame(
  study = c("Paller2022", "Simpson2020", "Paller2020", "Ebisawa2024", "Torrelo2023", "Paller2023"), # multiple studies would go here
  r     = c(28, 7, 33, 6, 39, 6),   # number of responders 
  n     = c(79, 85, 123, 32, 122, 94)    # sample size OF CONTROL GROUP ! ONLY 
  
)

# fit the prior (meta-analysis)
fit <- gMAP(
  cbind(r, n - r) ~ 1 | study, # ~ 1 | study means each study gets its own random effect
  data = data,
  family = binomial(),
  tau.dist = "HalfNormal",
  tau.prior = 1,   # scale of half-normal
  chains = 2 # change to 4 for default
)



# automix the prior
map_prior <- automixfit(fit, Nc = 1:3) # number of normal distributions


# output
print(map_prior)


plot(map_prior)



set.seed(1)
sim_effects <- rmix(map_prior, n = 1000)


hist(sim_effects, main = "Simulated effects from MAP prior", xlab = "Effect size")

jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures_Bayes/bayesian_forest.jpg", width = 600, height = 400)

forest_plot(fit)
dev.off()

print(map_prior)

n_active <- 30
n_placebo <- 30

robust_prior <- robustify(map_prior,weight=0.5) # 50% to comp 3 / vague component (new), other 50% dist among other comps

treat_prior <- mixbeta(c(1, 1, 1)) # prior for treatment used in trial # vague prior # first value is weight, second is shape a,
# next is shape b 

decision <- decision2S(0.975, 0, lower.tail = FALSE) # second number is delta
design_robust<- oc2S(treat_prior, robust_prior, n_active, n_placebo, decision) # numbers are sample size of each arm

true_placebo <- summary(map_prior)["mean"]
true_active <- true_placebo + 0.24 # value of interest / target (normally) # that 24 comes from clinical discussion

type1_robust <- design_robust(true_placebo, true_placebo) # type 1 error
power_robust <- design_robust(true_active, true_placebo) # power

ess_prior <- ess(robust_prior)  # how much weight prior carries in sample size


ess_prior_inital <- ess(map_prior)

plot(fit)
