library(synthpop)
library(dplyr)
library(ggplot2)



# set the relationship between covariates and response

# female has 1.5 times the odds of response compared with male
or_female <- 1.5

# under 2 has 2 times the odds of response compared with over 2
or_under_2 <- 2.0

# convert odds ratios to log odds coefficients
beta_female <- log(or_female)
beta_under_2 <- log(or_under_2)



# function to create ideal data

make_ideal_data <- function(
  n,
  responders,
  sex_female_n,
  age_over_2_n
) {

  # create sex variable

  sex_female <- c(
    rep(1, sex_female_n),
    rep(0, n - sex_female_n)
  )

  # create age variable

  age_over_2 <- c(
    rep(1, age_over_2_n),
    rep(0, n - age_over_2_n)
  )

  # randomise covariates

  sex_female <- sample(sex_female)
  age_over_2 <- sample(age_over_2)

  # create under 2 indicator

  age_under_2 <- 1 - age_over_2

  # calculate response weights

  # male + over 2 = baseline weight of 1
  # female + over 2 = 1.5
  # male + under 2 = 2
  # female + under 2 = 3

  log_weight <-
    beta_female * sex_female +
    beta_under_2 * age_under_2

  response_weight <- exp(log_weight)

  # select exactly the required number of responders

  # people with higher response weights are more likely
  # to be selected as responders

  responder_ids <- sample(
    1:n,
    size = responders,
    replace = FALSE,
    prob = response_weight
  )

  # create response variable

  response <- rep(0, n)

  response[responder_ids] <- 1

  # return dataset

  data.frame(
    response = factor(response),
    sex_female = factor(sex_female),
    age_over_2 = factor(age_over_2)
  )
}



# Paller 2022

Paller2022 <- make_ideal_data(
  n = 79,
  responders = 28,
  sex_female_n = 24,
  age_over_2_n = 74
)



# Simpson 2020

Simpson2020 <- make_ideal_data(
  n = 85,
  responders = 7,
  sex_female_n = 32,
  age_over_2_n = 85
)



# Paller 2020

Paller2020 <- make_ideal_data(
  n = 123,
  responders = 33,
  sex_female_n = 62,
  age_over_2_n = 123
)



# Ebisawa 2024

Ebisawa2024 <- make_ideal_data(
  n = 32,
  responders = 6,
  sex_female_n = 13,
  age_over_2_n = 32
)



# Torrelo 2023

Torrelo2023 <- make_ideal_data(
  n = 122,
  responders = 39,
  sex_female_n = 64,
  age_over_2_n = 122
)



# Paller 2023

Paller2023 <- make_ideal_data(
  n = 94,
  responders = 6,
  sex_female_n = 43,
  age_over_2_n = 94
)



# pooled response rate

total_r <- sum(c(
  28, 7, 33, 6, 39, 6
))

total_n <- sum(c(
  79, 85, 123, 32, 122, 94
))

true_control <- total_r / total_n

print(true_control)



# combine studies

mydata <- bind_rows(
  Ebisawa2024,
  Paller2020,
  Paller2022,
  Paller2023,
  Simpson2020,
  Torrelo2023
)



# check structure

str(mydata)



# check distributions

table(mydata$response)
table(mydata$sex_female)
table(mydata$age_over_2)



# check the relationship between sex and response

cat("\n")
cat("SEX BY RESPONSE - IDEAL DATA\n")

print(
  table(
    mydata$sex_female,
    mydata$response
  )
)

cat("\n")
cat("SEX BY RESPONSE - IDEAL DATA, ROW PROPORTIONS\n")

print(
  prop.table(
    table(
      mydata$sex_female,
      mydata$response
    ),
    margin = 1
  )
)



# check the relationship between age and response

cat("\n")
cat("AGE BY RESPONSE - IDEAL DATA\n")

print(
  table(
    mydata$age_over_2,
    mydata$response
  )
)

cat("\n")
cat("AGE BY RESPONSE - IDEAL DATA, ROW PROPORTIONS\n")

print(
  prop.table(
    table(
      mydata$age_over_2,
      mydata$response
    ),
    margin = 1
  )
)



# check the three-way relationship

cat("\n")
cat("SEX BY AGE BY RESPONSE - IDEAL DATA\n")

print(
  table(
    mydata$sex_female,
    mydata$age_over_2,
    mydata$response
  )
)



# sample size

# mean study size

n <- mean(c(
  nrow(Ebisawa2024),
  nrow(Paller2020),
  nrow(Paller2022),
  nrow(Paller2023),
  nrow(Simpson2020),
  nrow(Torrelo2023)
))

# for normal results use 89
# for ESS comparison to BDB use 33

#n <- 33
n <-89


# sim settings

nSim <- 10000

syn_results <- vector("list", length = nSim)

fisher_results <- vector("list", length = nSim)

smd_response <- vector("list", length = nSim)

smd_sex <- vector("list", length = nSim)

smd_age <- vector("list", length = nSim)

power <- vector("list", length = nSim)

ess <- vector("list", length = nSim)

alpha <- 0.05

type1_flags <- logical(nSim)

power_flags <- logical(nSim)

bias_scm <- numeric(nSim)

coverage_scm <- logical(nSim)

synthetic_control_rates <- numeric(nSim)

effect_size <- 0.24

n_treatment <- 30



# sim loop

for (i in 1:nSim) {

  set.seed(i)

  # synth new data

  codebook.syn(mydata)

  mysyn <- syn(
    mydata,
    k = n,
    visit.sequence = c(
      "sex_female",
      "age_over_2",
      "response"
    ),
    cont.na = NULL,
    minnumlevels = 2,
    maxfaclevels = 50
  )

  # save synth data

  mysyn <- mysyn$syn

  # store synth data

  syn_results[[i]] <- mysyn


  # synth control response estimate

  p_syn <- mean(
    as.numeric(as.character(mysyn$response)),
    na.rm = TRUE
  )

  synthetic_control_rates[i] <- p_syn


  # bias

  bias_scm[i] <- p_syn - true_control


  # 95% CI

  synthetic_response <- as.numeric(
    as.character(mysyn$response)
  )

  ci <- binom.test(
    sum(synthetic_response),
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


  # fischers exact

  original_response <- as.numeric(
    as.character(mydata$response)
  )

  synthetic_response <- as.numeric(
    as.character(mysyn$response)
  )

  tbl <- table(
    group = c(
      rep("original", nrow(mydata)),
      rep("synthetic", nrow(mysyn))
    ),
    response = c(
      original_response,
      synthetic_response
    )
  )

  fisher_test <- fisher.test(tbl)

  fisher_results[[i]] <- fisher_test$p.value


  # response smd

  sd1 <- sd(
    original_response,
    na.rm = TRUE
  )

  sd2 <- sd(
    synthetic_response,
    na.rm = TRUE
  )

  res_sd <- sqrt(
    (sd1^2 + sd2^2) / 2
  )

  res_md <- mean(
    original_response,
    na.rm = TRUE
  ) -
    mean(
      synthetic_response,
      na.rm = TRUE
    )

  smd_response[[i]] <- res_md / res_sd


  # sex smd

  original_sex <- as.numeric(
    as.character(mydata$sex_female)
  )

  synthetic_sex <- as.numeric(
    as.character(mysyn$sex_female)
  )

  sd1 <- sd(
    original_sex,
    na.rm = TRUE
  )

  sd2 <- sd(
    synthetic_sex,
    na.rm = TRUE
  )

  sex_sd <- sqrt(
    (sd1^2 + sd2^2) / 2
  )

  sex_md <- mean(
    original_sex,
    na.rm = TRUE
  ) -
    mean(
      synthetic_sex,
      na.rm = TRUE
    )

  smd_sex[[i]] <- sex_md / sex_sd


  # age smd

  original_age <- as.numeric(
    as.character(mydata$age_over_2)
  )

  synthetic_age <- as.numeric(
    as.character(mysyn$age_over_2)
  )

  sd1 <- sd(
    original_age,
    na.rm = TRUE
  )

  sd2 <- sd(
    synthetic_age,
    na.rm = TRUE
  )

  age_sd <- sqrt(
    (sd1^2 + sd2^2) / 2
  )

  age_md <- mean(
    original_age,
    na.rm = TRUE
  ) -
    mean(
      synthetic_age,
      na.rm = TRUE
    )

  smd_age[[i]] <- age_md / age_sd


  # power / ess

  n_control <- 30

  n_treatment <- 30

  # sample 30 from original data

  control_sample <- mydata[
    sample(nrow(mydata), n_control),
  ]

  control_response <- as.numeric(
    as.character(control_sample$response)
  )

  # control response rate

  num_responders <- sum(
    control_response == 1,
    na.rm = TRUE
  )

  p1 <- num_responders / n_control

  # hypothetical treatment response rate
  # adding a 0.24 effect size

  p2 <- min(
    p1 + 0.24,
    0.99
  )

  # power

  power_result <- power.prop.test(
    n = n_control,
    p1 = p1,
    p2 = p2,
    sig.level = alpha,
    alternative = "one.sided"
  )

  power[i] <- power_result$power


  # ESS

  ess_result <- power.prop.test(
    power = 0.8,
    p1 = p1,
    p2 = p2,
    sig.level = alpha,
    alternative = "one.sided"
  )

  ess[i] <- ess_result$n
}



# results

print(type1_flags)

type1_error_rate <- mean(
  type1_flags
)

power_rate <- mean(
  power_flags
)

power_mean <- mean(
  unlist(power),
  na.rm = TRUE
)

bias_mean <- mean(
  bias_scm,
  na.rm = TRUE
)

coverage_rate <- mean(
  coverage_scm
)

fidelity_difference_rate <- mean(
  unlist(fisher_results) < alpha,
  na.rm = TRUE
)

mean_smd_response <- mean(
  unlist(smd_response),
  na.rm = TRUE
)

mean_smd_sex <- mean(
  unlist(smd_sex),
  na.rm = TRUE
)

mean_smd_age <- mean(
  unlist(smd_age),
  na.rm = TRUE
)

mean_ess <- mean(
  unlist(ess),
  na.rm = TRUE
)



# print main results

cat("\n")
cat("SIMULATION RESULTS\n")
cat("----------------------------------------\n")

cat(
  "True pooled response rate:",
  round(true_control, 4),
  "\n"
)

cat(
  "Mean synthetic response rate:",
  round(true_control + bias_mean, 4),
  "\n"
)

cat(
  "Mean bias:",
  round(bias_mean, 4),
  "\n"
)

cat(
  "95% CI coverage:",
  round(coverage_rate, 4),
  "\n"
)

cat(
  "Type I error rate:",
  round(type1_error_rate, 4),
  "\n"
)

cat(
  "Empirical power:",
  round(power_rate, 4),
  "\n"
)

cat(
  "Mean power:",
  round(power_mean, 4),
  "\n"
)

cat(
  "Mean ESS:",
  round(mean_ess, 2),
  "\n"
)

cat(
  "Original versus synthetic significant difference:",
  round(fidelity_difference_rate, 4),
  "\n"
)

cat("\n")
cat("STANDARDISED MEAN DIFFERENCES\n")
cat("----------------------------------------\n")

cat(
  "Response SMD:",
  round(mean_smd_response, 4),
  "\n"
)

cat(
  "Sex SMD:",
  round(mean_smd_sex, 4),
  "\n"
)

cat(
  "Age SMD:",
  round(mean_smd_age, 4),
  "\n"
)



# smd hist

# response SMD

df_response_smd <- data.frame(
  smd = unlist(smd_response)
)

ggplot(
  df_response_smd,
  aes(x = smd)
) +
  geom_histogram(
    bins = 20
  ) +
  labs(
    title = "Synthetic Control: Response SMD",
    x = "SMD",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank()
  )



# sex SMD

df_sex_smd <- data.frame(
  smd = unlist(smd_sex)
)

ggplot(
  df_sex_smd,
  aes(x = smd)
) +
  geom_histogram(
    bins = 20
  ) +
  labs(
    title = "Synthetic Control: Sex SMD",
    x = "SMD",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank()
  )



# age SMD

df_age_smd <- data.frame(
  smd = unlist(smd_age)
)

ggplot(
  df_age_smd,
  aes(x = smd)
) +
  geom_histogram(
    bins = 20
  ) +
  labs(
    title = "Synthetic Control: Age SMD",
    x = "SMD",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank()
  )



# fisher p value hist

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



# check og vs synth covariates

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

final_syn <- syn_results[[nSim]]



# response

cat("\n")
cat("RESPONSE DISTRIBUTION\n")
cat("----------------------------------------\n")

print(
  prop.table(
    table(mydata$response)
  )
)

print(
  prop.table(
    table(final_syn$response)
  )
)



# sex

cat("\n")
cat("SEX DISTRIBUTION\n")
cat("----------------------------------------\n")

print(
  prop.table(
    table(mydata$sex_female)
  )
)

print(
  prop.table(
    table(final_syn$sex_female)
  )
)



# age

cat("\n")
cat("AGE DISTRIBUTION\n")
cat("----------------------------------------\n")

print(
  prop.table(
    table(mydata$age_over_2)
  )
)

print(
  prop.table(
    table(final_syn$age_over_2)
  )
)



# joint dist

cat("\n")
cat("SEX BY RESPONSE - IDEAL DATA\n")

print(
  prop.table(
    table(
      mydata$sex_female,
      mydata$response
    ),
    margin = 1
  )
)

cat("\n")
cat("SEX BY RESPONSE - SYNTHETIC DATA\n")

print(
  prop.table(
    table(
      final_syn$sex_female,
      final_syn$response
    ),
    margin = 1
  )
)


cat("\n")
cat("AGE BY RESPONSE - IDEAL DATA\n")

print(
  prop.table(
    table(
      mydata$age_over_2,
      mydata$response
    ),
    margin = 1
  )
)

cat("\n")
cat("AGE BY RESPONSE - SYNTHETIC DATA\n")

print(
  prop.table(
    table(
      final_syn$age_over_2,
      final_syn$response
    ),
    margin = 1
  )
)



# three way relationship

cat("\n")
cat("SEX BY AGE BY RESPONSE - IDEAL DATA\n")

print(
  table(
    mydata$sex_female,
    mydata$age_over_2,
    mydata$response
  )
)

cat("\n")
cat("SEX BY AGE BY RESPONSE - SYNTHETIC DATA\n")

print(
  table(
    final_syn$sex_female,
    final_syn$age_over_2,
    final_syn$response
  )
)



# forest plot

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
      as.numeric(
        as.character(df$response)
      ) == 1,
      na.rm = TRUE
    )

    n_study <- nrow(df)

    prop <- r / n_study

    ci <- binom.test(
      r,
      n_study
    )$conf.int

    data.frame(
      study = name,
      response_rate = prop,
      ci_lower = ci[1],
      ci_upper = ci[2],
      n = n_study
    )
  }
) %>%
  bind_rows()



# synth control

syn_response <- sum(
  as.numeric(
    as.character(final_syn$response)
  ) == 1,
  na.rm = TRUE
)

syn_n <- nrow(final_syn)

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



# overall pooled estimate

original_studies <- study_data[
  study_data$study != "Synthetic Control",
]


total_r <- sum(
  original_studies$response_rate *
    original_studies$n
)

total_n <- sum(
  original_studies$n
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



# forest

jpeg(
  "C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures_Bayes/synthetic_forest_cov_fake.jpg",
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
    title = "Forest Plot of Control Group Response Rates - Ideal Covariates",
    x = "Response Rate (95% CI)",
    y = "Study"
  ) +

  xlim(
    c(0.0, 0.8)
  ) +

  theme_minimal() +

  theme(
    panel.grid = element_blank(),
    axis.line = element_line(
      color = "black"
    )
  )


dev.off()




# results to report
print(
  table(coverage_scm)
)

final_syn_response <- as.numeric(
  as.character(final_syn$response)
)

final_syn_rate <- mean(
  final_syn_response,
  na.rm = TRUE
)

final_syn_ci <- binom.test(
  sum(final_syn_response),
  length(final_syn_response)
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
