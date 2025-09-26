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
n <-  123
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








mydata <- bind_rows(Ebisawa2024, Paller2020, Paller2022, Paller2023, Simpson2020, Torrelo2023)
mydata <- cbind(group = rep(0, nrow(mydata)), mydata)
n <- mean(c(nrow(Ebisawa2024), nrow(Paller2020), nrow(Paller2022), nrow(Paller2023), nrow(Simpson2020), nrow(Torrelo2023)))
nSim <- 10000
syn_results <- vector("list", length = nSim) # create empty df to store results
fisher_results <- vector("list", length = nSim) 
smd_response <- vector("list", length = nSim) 
power <- vector("list", length = nSim)
ess <- vector("list", length = nSim)

alpha <- 0.05  # significance threshold
type1_flags <- logical(nSim)  # store whether p < alpha in each iteration

for (i in 1:nSim) {
  set.seed(i) 
  # set n to be equal to length of original data (adjustable)
  
  
  
  # synthesise new data
  
  codebook.syn(mydata)
  
  
  # minimum number of observations is set to 5
  # k is the number of synthetic observations (n)
  # default syn method is CART
  mysyn <- syn(mydata, k = n, cont.na = NULL, minnumlevels = 2, maxfaclevels = 50)
  
  summary(mysyn)
  
  # store results
  syn_results[[i]] <- mysyn$syn
  
  # make mysyn into readable df
  mysyn <- mysyn$syn
  
  # test power, type 1 error
  # 
  tbl <- table(
    group = c(rep("original", nrow(mydata)), rep("synthetic", nrow(mysyn))),
    response = c(mydata$response, mysyn$response)
  )
  
  # Fisher's exact test 
  fisher_test <- fisher.test(tbl)
  fisher_results[[i]] <- fisher_test$p.value
  
  # Flag if p < alpha
  type1_flags[i] <- fisher_test$p.value < alpha
  
  # standard mean difference - smd
  
  sd1 <- sd(mydata$response, na.rm = T) 
  sd2 <- sd(mysyn$response, na.rm = T)
  res_sd <- sqrt((sd1^2 + sd2^2)/2)
  res_md <- mean(mydata$response, na.rm = T)-mean(mysyn$response, na.rm = T)
  smd_response[[i]] <- res_md/res_sd
  
  n_control <- 30
  n_treatment <- 30
  
  # sample 30 from original data
  control_sample <- mydata[sample(nrow(mydata), n_control), ]
  
  #  control response rate
  num_responders <- sum(control_sample$response == 1, na.rm = TRUE)
  p1 <- num_responders / n_control
  
  # hypothetical treatment response rate (adding a 0.24 effect size)
  p2 <- min(p1 + 0.24, 0.99)  # ensure it doesn't exceed 1
  
  

  power_result <-power.prop.test(n = n_control, p1 = p1, p2 = p2, sig.level = alpha, alternative = "one.sided")
  power[i] <- power_result$power
  ess_result <- power.prop.test(power = 0.8, p1 = p1, p2 = p2, sig.level = alpha) # change needed power
  ess[i] <- ess_result$n #pointless measurement but included anyway
}


print(type1_flags)
type1_error_rate <- mean(type1_flags)
power_mean <- mean(unlist(power), na.rm = TRUE)


# graphing:


df <- data.frame(smd = unlist(smd_response))
ggplot(df, aes(x = smd)) +
  geom_histogram(fill = "steelblue", color = "white", bins = 20) +
  labs(title = "Simulated Treatment Effects (SMD)", x = "SMD", y = "Count") +
  theme_minimal() +
  theme(panel.grid = element_blank())



hist(unlist(fisher_results), 
     main = "Distribution of P-values (Synthetic Control)", 
     xlab = "P-value", 
     col = "salmon", border = "white")
abline(v = 0.05, col = "red", lty = 2)










# forest plot:


study_list <- list(
  Paller2022 = Paller2022,
  Simpson2020 = Simpson2020,
  Paller2020 = Paller2020,
  Ebisawa2024 = Ebisawa2024,
  Torrelo2023 = Torrelo2023,
  Paller2023 = Paller2023
)



study_data <- lapply(names(study_list), function(name) {
  df <- study_list[[name]]
  r <- sum(df$response == 1, na.rm = TRUE)
  n <- nrow(df)
  prop <- r / n
  ci <- binom.test(r, n)$conf.int
  data.frame(
    study = name,
    response_rate = prop,
    ci_lower = ci[1],
    ci_upper = ci[2],
    n = n
  )
}) %>% bind_rows()

syn_response <- sum(mysyn$response == 1, na.rm = TRUE)
syn_n <- nrow(mysyn)
syn_prop <- syn_response / syn_n
syn_ci <- binom.test(syn_response, syn_n)$conf.int


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

# add overall pooled estimate (original studies only)
total_r <- sum(study_data$response_rate[study_data$study != "Synthetic Control"] * study_data$n[study_data$study != "Synthetic Control"])
total_n <- sum(study_data$n[study_data$study != "Synthetic Control"])
overall_rate <- total_r / total_n
overall_ci <- binom.test(round(total_r), total_n)$conf.int

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

study_data$study <- factor(study_data$study,
                           levels = c("Synthetic Control",
                                      "Pooled Mean","Paller2023", "Torrelo2023", "Ebisawa2024",
                                      "Paller2020", "Simpson2020", "Paller2022" ))


synthetic_rate <- study_data$response_rate[study_data$study == "Synthetic Control"]

jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures_Bayes/synthetic_forest.jpg", width = 600, height = 400)

ggplot(study_data, aes(x = response_rate, y = study)) +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2, color = "steelblue") +
  geom_vline(xintercept = synthetic_rate, linetype = "dashed", color = "blue") + 
  labs(
    title = "Forest Plot of Control Group Response Rates",
    x = "Response Rate (95% CI)",
    y = "Study"
  ) +
  xlim(c(0.0, 0.8))+
  theme_minimal() +
  
          theme(
            panel.grid = element_blank(),
            axis.line = element_line(color = "black") 
          )
dev.off()

