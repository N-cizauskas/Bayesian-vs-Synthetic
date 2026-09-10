library(RBesT)
library(dplyr)
library(ggplot2)

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
sim_effects <- rmix(map_prior, n = 10000)


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

# bias and coverage metrics:

# posterior mean from robust prior
posterior_mean <- summary(robust_prior)["mean"]

# bias
bias_bdb <- posterior_mean - true_placebo

# draw samples from robust posterior (approximation)
samples <- rmix(robust_prior, n = 10000)

# 95% credible interval
ci <- quantile(samples, c(0.025, 0.975))

# check coverage
coverage_bdb <- (true_placebo >= ci[1] & true_placebo <= ci[2])







# power curve plot:
# plausible range of placebo response rates
placebo_rates <- seq(0.06, 0.36, by = 0.01)

# calculate type I error and power at each placebo response rate
type1_curve <- numeric(length(placebo_rates))
power_curve <- numeric(length(placebo_rates))

for (i in seq_along(placebo_rates)) {

  p_placebo <- placebo_rates[i]

  # type 1 error:
  # true treatment effect = 0
  p_active_null <- p_placebo

  type1_curve[i] <- design_robust(
    p_active_null,
    p_placebo
  )

  # power:
  # true treatment effect = 0.24
  p_active_alt <- p_placebo + 0.24

  # only calculate if active response rate is <= 1
  if (p_active_alt <= 1) {
    power_curve[i] <- design_robust(
      p_active_alt,
      p_placebo
    )
  } else {
    power_curve[i] <- NA
  }
}



oc_curve <- data.frame(
  placebo_rate = placebo_rates,
  type1_error = type1_curve,
  power = power_curve
)

print(oc_curve)

jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures_Bayes/oc_curve_power.jpg", width = 600, height = 400)




plot(
  oc_curve$placebo_rate,
  oc_curve$power,
  type = "l",
  lwd = 2,
  xlab = "Placebo response rate",
  ylab = "Power",
  ylim = c(0, 1),
  main = "BDB Power Across Placebo Response Rates"
)

# mark  original analysis point
points(
  true_placebo,
  power_robust,
  pch = 19
)

abline(
  v = true_placebo,
  lty = 2
)

dev.off()


jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures_Bayes/oc_curve_t1e.jpg", width = 600, height = 400)



plot(
  oc_curve$placebo_rate,
  oc_curve$type1_error,
  type = "l",
  lwd = 2,
  xlab = "Placebo response rate",
  ylab = "Type I error",
  ylim = c(0, 1),
  main = "BDB Type I Error Across Placebo Response Rates"
)

# mark  original analysis point
points(
  true_placebo,
  type1_robust,
  pch = 19
)

abline(
  v = true_placebo,
  lty = 2
)
dev.off()


# reporting results

map_samples <- rmix(map_prior, n = 10000)
robust_samples <- rmix(robust_prior, n = 10000)

map_ci <- quantile(map_samples, c(.025, .975))
robust_ci <- quantile(robust_samples, c(.025, .975))

design_map <- oc2S(treat_prior, map_prior, n_active, n_placebo, decision)

data.frame(
  Method = c("MAP", "Robustified MAP"),
  Mean = c(summary(map_prior)["mean"], summary(robust_prior)["mean"]),
  CI_lower = c(map_ci[1], robust_ci[1]),
  CI_upper = c(map_ci[2], robust_ci[2]),
  Power = c(design_map(true_active, true_placebo), power_robust),
  Type1_Error = c(design_map(true_placebo, true_placebo), type1_robust),
  Bias = c(summary(map_prior)["mean"] - true_placebo, bias_bdb),
  Coverage = c(
    true_placebo >= map_ci[1] & true_placebo <= map_ci[2],
    true_placebo >= robust_ci[1] & true_placebo <= robust_ci[2]
  ),
  ESS = c(ess(map_prior), ess(robust_prior))
)




# bdb coverage:
set.seed(123)

nSim <- 1000

true_control <- sum(data$r) / sum(data$n)

coverage_map <- logical(nSim)
coverage_robust <- logical(nSim)

bias_map <- numeric(nSim)
bias_robust <- numeric(nSim)

posterior_mean_map <- numeric(nSim)
posterior_mean_robust <- numeric(nSim)

map_lower <- numeric(nSim)
map_upper <- numeric(nSim)

robust_lower <- numeric(nSim)
robust_upper <- numeric(nSim)

for (i in 1:nSim) {

  # Simulate the new placebo arm
  sim_r <- rbinom(
    n = 1,
    size = n_placebo,
    prob = true_control
  )

  # Combine the prior with the new trial data
  post_map <- postmix(
    map_prior,
    r = sim_r,
    n = n_placebo
  )

  post_robust <- postmix(
    robust_prior,
    r = sim_r,
    n = n_placebo
  )

  # Draw samples from the posterior distributions
  map_samples <- rmix(post_map, n = 10000)
  robust_samples <- rmix(post_robust, n = 10000)

  # 95% posterior credible intervals
  map_ci <- quantile(
    map_samples,
    c(0.025, 0.975)
  )

  robust_ci <- quantile(
    robust_samples,
    c(0.025, 0.975)
  )

  # Posterior means
  map_mean <- summary(post_map)["mean"]
  robust_mean <- summary(post_robust)["mean"]

  posterior_mean_map[i] <- map_mean
  posterior_mean_robust[i] <- robust_mean

  # Bias
  bias_map[i] <- map_mean - true_control
  bias_robust[i] <- robust_mean - true_control

  # Store intervals
  map_lower[i] <- map_ci[1]
  map_upper[i] <- map_ci[2]

  robust_lower[i] <- robust_ci[1]
  robust_upper[i] <- robust_ci[2]

  # Coverage
  coverage_map[i] <-
    true_control >= map_ci[1] &&
    true_control <= map_ci[2]

  coverage_robust[i] <-
    true_control >= robust_ci[1] &&
    true_control <= robust_ci[2]
}

print(data.frame(
  Method = c("MAP", "Robustified MAP"),
  Posterior_Mean = c(
    mean(posterior_mean_map),
    mean(posterior_mean_robust)
  ),
  Bias = c(
    mean(bias_map),
    mean(bias_robust)
  ),
  Coverage = c(
    mean(coverage_map),
    mean(coverage_robust)
  ),
  Mean_CI_Width = c(
    mean(map_upper - map_lower),
    mean(robust_upper - robust_lower)
  )
))

print(mean(map_upper))
print(mean(map_lower))
print(mean(robust_upper))
print(mean(robust_lower))































# forest graph including posterior mean and cI



ci_study1 <- binom.test(
  28,
  79
)$conf.int
ci_study2 <- binom.test(
  7,
  85
)$conf.int
ci_study3 <- binom.test(
  33,
  123
)$conf.int
ci_study4 <- binom.test(
  6,
  32
)$conf.int
ci_study5 <- binom.test(
  39,
  122
)$conf.int
ci_study6 <- binom.test(
  6,
  94
)$conf.int

ci_lower_list<-c(ci_study1[1], ci_study2[1], ci_study3[1], ci_study4[1], ci_study5[1], ci_study6[1])
ci_upper_list<-c(ci_study1[2], ci_study2[2], ci_study3[2], ci_study4[2], ci_study5[2], ci_study6[2])






study_data<- data.frame(
  study=data$study,
  response_rate = data$r / data$n,
  ci_lower=ci_lower_list,
  ci_upper=ci_upper_list,
  n=data$n
)




study_data <- rbind(
  study_data,
  data.frame(
    study="Posterior",
    response_rate = 0.22,
    ci_lower=0.11,
    ci_upper=0.36,
    n=30
  )
)

study_data <- rbind(
  study_data,
  data.frame(
    study="Prior",
    response_rate =0.24 ,
    ci_lower=0.02,
    ci_upper=0.69,
    n=30
  )
)

study_data <- rbind(
  study_data,
  data.frame(
    study="Synthetic Control",
    response_rate =0.26 ,
    ci_lower=0.17,
    ci_upper=0.36,
    n=89
  )
)

study_data <- rbind(
  study_data,
  data.frame(
    study="Synthetic Control (Size Matched)",
    response_rate =0.27 ,
    ci_lower=0.13,
    ci_upper=0.46,
    n=33
  )
)




study_data <- rbind(
  study_data,
  data.frame(
    study="Posterior (Robustified)",
    response_rate =0.22 ,
    ci_lower=0.11,
    ci_upper=0.37,
    n=30
  )
)

study_data$study <- factor(
  study_data$study,
  levels = c(
    "Synthetic Control (Size Matched)",
    "Synthetic Control",
    "Posterior (Robustified)",
    "Posterior",
    "Prior",
    "Paller2023",
    "Torrelo2023",
    "Ebisawa2024",
    "Paller2020",
    "Simpson2020",
    "Paller2022"
  ))

jpeg(
  "C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures_Bayes/everything_forest_post.jpg",
  width = 5.5, height = 6.5, units = "in", res = 300
)


ggplot(
  study_data,
  aes(
    x = response_rate,
    y = study
  )
) +
  geom_point(
    size = 3,
    color = "chartreuse3"
  ) +
  geom_errorbarh(
    aes(
      xmin = ci_lower,
      xmax = ci_upper
    ),
    height = 0.2,
    color = "chartreuse3"
  ) +
  geom_vline(
    xintercept = 0.222, # using this number because it is the MEAN response rate of the pooled historical studies!!
    linetype = "dashed",
    color = "darkseagreen3"
  ) +
  labs(
    title = "Forest Plot of Response Rates",
    x = "Response Rate (95% CI)",
    y = "Study"
  ) +
  xlim(
    c(0.0, 0.8)
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black")
  )

dev.off()



print(data.frame(
  Method = c("MAP", "Robustified MAP"),
  Prior_Mean = c(
    summary(map_prior)["mean"],
    summary(robust_prior)["mean"]
  ),
  Posterior_Mean = c(
    summary(post_map)["mean"],
    summary(post_robust)["mean"]
  )
))


