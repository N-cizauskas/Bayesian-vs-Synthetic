library(synthpop)
library(dplyr)
library(ggplot2)

# for paller 2022:
n <- 79
responders <- 28

response <- c(rep(1, responders), rep(0, n - responders))
response <- sample(response)  # randomise across rows

Paller2022 <- data.frame(
  response = response
)


# for simpson 2020:
n <- 85
responders <- 7

response <- c(rep(1, responders), rep(0, n - responders))
response <- sample(response)  # randomise across rows

Simpson2020 <- data.frame(
  response = response
)


# for paller 2020:
n <- 123
responders <- 33

response <- c(rep(1, responders), rep(0, n - responders))
response <- sample(response)  # randomise across rows

Paller2020 <- data.frame(
  response = response
)


# for ebisawa 2024:
n <- 32
responders <- 6

response <- c(rep(1, responders), rep(0, n - responders))
response <- sample(response)  # randomise across rows

Ebisawa2024 <- data.frame(
  response = response
)


# for torrelo 2023:
n <- 122
responders <- 39

response <- c(rep(1, responders), rep(0, n - responders))
response <- sample(response)  # randomise across rows

Torrelo2023 <- data.frame(
  response = response
)


# for paller 2023:
n <- 94
responders <- 6

response <- c(rep(1, responders), rep(0, n - responders))
response <- sample(response)  # randomise across rows

Paller2023 <- data.frame(
  response = response
)


# pooled historical response rate

total_r <- sum(c(28, 7, 33, 6, 39, 6))
total_n <- sum(c(79, 85, 123, 32, 122, 94))

true_control <- total_r / total_n


mydata <- bind_rows(
  Ebisawa2024,
  Paller2020,
  Paller2022,
  Paller2023,
  Simpson2020,
  Torrelo2023
)

mydata <- cbind(
  group = rep(0, nrow(mydata)),
  mydata
)


# synthetic control sample size
# run first with n = 89 for the main SCM analysis
# then rerun with n = 33 for comparison with BDB

n <- 89

#n <- 33


nSim <- 10000

syn_results <- vector("list", length = nSim)

# p-values comparing original historical data with synthetic data
fisher_results <- vector("list", length = nSim)

# standardised mean difference
smd_response <- vector("list", length = nSim)

# store type 1 error and power decisions
type1_flags <- logical(nSim)
power_flags <- logical(nSim)

# store bias and coverage
bias_scm <- numeric(nSim)
coverage_scm <- logical(nSim)

# store synthetic control response rates
synthetic_control_rates <- numeric(nSim)

successful_flags <- logical(nSim)

alpha <- 0.05

# assumed treatment effect
effect_size <- 0.24

# treatment arm sample size
n_treatment <- 30


for (i in 1:nSim) {

  set.seed(i)

  # synthesise new data

  codebook.syn(mydata)


  # minimum number of observations is set to 5
  # k is the number of synthetic observations
  # default syn method is CART
  # synthesise new data

  mysyn <- tryCatch(
    syn(
      mydata,
      k = n,
      cont.na = NULL,
      minnumlevels = 2,
      maxfaclevels = 50
    ),
    error = function(e) {
      message(
        "ERROR AT SIMULATION: ",
        i,
        "\n",
        e$message
      )
      return(NULL)
    }
  )

  # skip failed synthesis

  if (is.null(mysyn)) {
    next
  }

  summary(mysyn)


  # make mysyn into readable df

  mysyn <- mysyn$syn

  synthetic_response <- as.numeric(
    as.character(mysyn$response)
  )
  # store results

  syn_results[[i]] <- mysyn


  # synthetic control estimate

  p_syn <- mean(mysyn$response, na.rm = TRUE)

  synthetic_control_rates[i] <- p_syn


  # bias

  bias_scm[i] <- p_syn - true_control


  # 95% CI for synthetic control

  ci <- binom.test(
    sum(mysyn$response),
    nrow(mysyn)
  )$conf.int


  # coverage

  coverage_scm[i] <- (
    true_control >= ci[1] &
      true_control <= ci[2]
  )


  # empirical type 1 error
  # under the null hypothesis theta = 0
  # treatment response rate is equal to synthetic control response rate

  p_treatment_null <- p_syn

  treatment_null <- rbinom(
    n = n_treatment,
    size = 1,
    prob = p_treatment_null
  )


  # treatment versus synthetic control under the null

  null_table <- table(
    group = c(
      rep("treatment", n_treatment),
      rep("synthetic", nrow(mysyn))
    ),
    response = factor(
      c(
        treatment_null,
        synthetic_response
      ),
      levels = c(0, 1)
    )
  )


  # one-sided Fisher's exact test
  # alternative hypothesis is that treatment response is greater than synthetic control

  null_test <- fisher.test(
    null_table,
    alternative = "greater"
  )


  # flag false positive

  type1_flags[i] <- null_test$p.value < alpha


  # empirical power
  # under the alternative hypothesis theta = 0.24
  # treatment response rate is 0.24 higher than synthetic control

  p_treatment_alt <- min(
    p_syn + effect_size,
    0.99
  )


  treatment_alt <- rbinom(
    n = n_treatment,
    size = 1,
    prob = p_treatment_alt
  )


  # treatment versus synthetic control under the alternative

  alternative_table <- table(
    group = c(
      rep("treatment", n_treatment),
      rep("synthetic", nrow(mysyn))
    ),
    response = factor(
      c(
        treatment_alt,
        synthetic_response
      ),
      levels = c(0, 1)
    )
  )


  # one-sided Fisher's exact test

  alternative_test <- fisher.test(
    alternative_table,
    alternative = "greater"
  )


  # flag successful treatment conclusion

  power_flags[i] <- alternative_test$p.value < alpha


  # compare original pooled data with synthetic data
  # this is a fidelity check rather than a type 1 error calculation

  tbl <- table(
    group = c(
      rep("original", nrow(mydata)),
      rep("synthetic", nrow(mysyn))
    ),
    response = c(
      mydata$response,
      mysyn$response
    )
  )


  # Fisher's exact test for original versus synthetic data

  fisher_test <- fisher.test(tbl)

  fisher_results[[i]] <- fisher_test$p.value


  # standard mean difference - smd

  sd1 <- sd(mydata$response, na.rm = TRUE)
  sd2 <- sd(mysyn$response, na.rm = TRUE)

  res_sd <- sqrt((sd1^2 + sd2^2) / 2)

  res_md <- mean(mydata$response, na.rm = TRUE) -
    mean(mysyn$response, na.rm = TRUE)

  smd_response[[i]] <- res_md / res_sd
  successful_flags[i] <- TRUE
}

successful_flags <- vapply(
  syn_results,
  function(x) {
    !is.null(x)
  },
  logical(1)
)
type1_error_rate <- mean(
  type1_flags[successful_flags]
)

power_rate <- mean(
  power_flags[successful_flags]
)

bias_mean <- mean(
  bias_scm[successful_flags],
  na.rm = TRUE
)

coverage_rate <- mean(
  coverage_scm[successful_flags]
)

# calculate proportion of simulations where
# original and synthetic response rates differed significantly

fidelity_difference_rate <- mean(
  unlist(fisher_results) < alpha,
  na.rm = TRUE
)


# print results

cat("\n")
cat("SCM simulation results\n")
cat("----------------------\n")
cat("synthetic control sample size:", n, "\n")
cat("number of simulations:", nSim, "\n")
cat("treatment sample size:", n_treatment, "\n")
cat("assumed treatment effect:", effect_size, "\n")
cat("\n")
cat("type 1 error:", round(type1_error_rate, 3), "\n")
cat("power:", round(power_rate, 3), "\n")
cat("mean bias:", round(bias_mean, 4), "\n")
cat("coverage:", round(coverage_rate, 3), "\n")
cat(
  "original versus synthetic significant difference:",
  round(fidelity_difference_rate, 3),
  "\n"
)

# use the final synthetic dataset generated

valid_syn <- which(
  vapply(
    syn_results,
    function(x) {
      !is.null(x) &&
        is.data.frame(x) &&
        nrow(x) > 0
    },
    logical(1)
  )
)

if (length(valid_syn) == 0) {
  stop("No valid synthetic datasets were generated.")
}

last_valid <- max(valid_syn)

final_syn <- syn_results[[last_valid]]


# graphing:


# distribution of synthetic control response rates

hist(
  synthetic_control_rates,
  main = paste(
    "Synthetic Control Response Rates (n =",
    n,
    ")"
  ),
  xlab = "Response rate",
  col = "steelblue",
  border = "white",
  breaks = 20
)

abline(
  v = true_control,
  col = "red",
  lty = 2,
  lwd = 2
)


# standardised mean difference

df <- data.frame(
  smd = unlist(smd_response)
)

ggplot(df, aes(x = smd)) +
  geom_histogram(
    fill = "steelblue",
    color = "white",
    bins = 20
  ) +
  labs(
    title = "Simulated Treatment Effects (SMD)",
    x = "SMD",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank()
  )


# distribution of p-values for original versus synthetic data

hist(
  unlist(fisher_results),
  main = "Distribution of P-values (Synthetic Control)",
  xlab = "P-value",
  col = "salmon",
  border = "white"
)

abline(
  v = 0.05,
  col = "red",
  lty = 2
)


# forest plot:


study_list <- list(
  Paller2022 = Paller2022,
  Simpson2020 = Simpson2020,
  Paller2020 = Paller2020,
  Ebisawa2024 = Ebisawa2024,
  Torrelo2023 = Torrelo2023,
  Paller2023 = Paller2023
)


study_data <- lapply(
  names(study_list),
  function(name) {

    df <- study_list[[name]]

    r <- sum(
      df$response == 1,
      na.rm = TRUE
    )

    n <- nrow(df)

    prop <- r / n

    ci <- binom.test(
      r,
      n
    )$conf.int

    data.frame(
      study = name,
      response_rate = prop,
      ci_lower = ci[1],
      ci_upper = ci[2],
      n = n
    )
  }
) %>%
  bind_rows()


syn_response <- sum(
  mysyn$response == 1,
  na.rm = TRUE
)

syn_n <- nrow(mysyn)

syn_prop <- syn_response / syn_n

syn_ci <- binom.test(
  syn_response,
  syn_n
)$conf.int


study_data <- rbind(
  study_data,
  data.frame(
    study = "Synthetic Control",
    response_rate = syn_prop,
    ci_lower = syn_ci[1],
    ci_upper = syn_ci[2],
    n = syn_n
  )
)


# add overall pooled estimate
# original studies only

total_r <- sum(
  study_data$response_rate[
    study_data$study != "Synthetic Control"
  ] *
    study_data$n[
      study_data$study != "Synthetic Control"
    ]
)

total_n <- sum(
  study_data$n[
    study_data$study != "Synthetic Control"
  ]
)

overall_rate <- total_r / total_n

overall_ci <- binom.test(
  round(total_r),
  total_n
)$conf.int


study_data <- rbind(
  study_data,
  data.frame(
    study = "Pooled Mean",
    response_rate = overall_rate,
    ci_lower = overall_ci[1],
    ci_upper = overall_ci[2],
    n = total_n
  )
)


study_data$study <- factor(
  study_data$study,
  levels = c(
    "Synthetic Control",
    "Pooled Mean",
    "Paller2023",
    "Torrelo2023",
    "Ebisawa2024",
    "Paller2020",
    "Simpson2020",
    "Paller2022"
  )
)


synthetic_rate <- study_data$response_rate[
  study_data$study == "Synthetic Control"
]


jpeg(
  "C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures_Bayes/synthetic_forest_new.jpg",
  width = 600,
  height = 400
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
    color = "steelblue"
  ) +
  geom_errorbarh(
    aes(
      xmin = ci_lower,
      xmax = ci_upper
    ),
    height = 0.2,
    color = "steelblue"
  ) +
  geom_vline(
    xintercept = synthetic_rate,
    linetype = "dashed",
    color = "blue"
  ) +
  labs(
    title = "Forest Plot of Control Group Response Rates",
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
# results to report

# final synthetic control results

final_syn_response <- as.numeric(
  as.character(final_syn$response)
)

final_syn_rate <- mean(
  final_syn_response,
  na.rm = TRUE
)

final_syn_ci <- binom.test(
  sum(final_syn_response, na.rm = TRUE),
  sum(!is.na(final_syn_response))
)$conf.int


final_results <- data.frame(
  synthetic_response_rate = final_syn_rate,
  ci_lower = final_syn_ci[1],
  ci_upper = final_syn_ci[2],
  power = power_rate,
  bias = bias_mean,
  type1_error = type1_error_rate,
  coverage = coverage_rate
)

print(final_results)

cat(
  "successful simulations:",
  sum(successful_flags),
  "\n"
)
