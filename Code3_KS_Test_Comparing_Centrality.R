

library(gtools)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)


##' datasets: central_data and noncentral data with columns:
##'       Transition1 e.g. "KH"
##'       Transition2 e.g. "RC"
##'       Transition e.g. "x: KH y: RC"
##'       Frequency e.g. "1"
##'       Sequence_length e.g. "5"
##'       x_position e.g. "1"
##'       y_position e.g. "2"
##'       Model e.g. "Tiger"
##'       Distance e.g. "1"

### generate transition probabilities using higher order markovian analysis

# read the datasets
central_data <- read.csv("/xxx/xxx/xxx/xxx/xxx/xxx/xxx.csv")
noncentral_data <- read.csv("/xxx/xxx/xxx/xxx/xxx/xxx/xxx.csv")

# aggregate frequencies
aggregate_central <- aggregate(Frequency ~ Transition1 + Transition2 + Distance, 
                               data = central_data, FUN = sum)
aggregate_noncentral <- aggregate(Frequency ~ Transition1 + Transition2 + Distance, 
                                  data = noncentral_data, FUN = sum)

# compute transition probabilities
central_probs <- aggregate_central
central_probs$Probability <- 
  central_probs$Frequency / sum(central_probs$Frequency)

noncentral_probs <- aggregate_noncentral
noncentral_probs$Probability <- 
  noncentral_probs$Frequency / sum(noncentral_probs$Frequency)

# step 4: comparative analysis
# compare transition probabilities between central and noncentral datasets
comparison <- merge(central_probs, noncentral_probs, by = 
                      c("Transition1", "Transition2", "Distance"), 
                    suffixes = c("_central", "_noncentral"))






### compare central & non-central unique frequencies using a X^2 test 

# create contingency table
central_table <- table(central_data$Transition1, central_data$Transition2)
noncentral_table <- table(noncentral_data$Transition1, noncentral_data$Transition2)

# perform chi-square test
chi_sq_result <- chisq.test(central_table, noncentral_table)
print(chi_sq_result)

# assuming you already ran this:
chi_sq_result <- chisq.test(central_table, noncentral_table)

# extract needed components
chi_sq <- chi_sq_result$statistic
n <- sum(central_table) + sum(noncentral_table)
k <- min(nrow(central_table), ncol(central_table))  # Use the smaller dimension

# calculate Cramér’s V (effect size)
cramers_v <- sqrt(chi_sq / (n * (k - 1)))
cat("Cramér's V:", cramers_v, "\n")

df <- (nrow(central_table) - 1) * (ncol(central_table) - 1)
cat("Degrees of Freedom:", df, "\n")

# lower and upper chi-squared bounds (95% CI)
chi_sq_lower <- qchisq(0.025, df)
chi_sq_upper <- qchisq(0.975, df)

# calculate lower and upper bounds for Cramér's V
cramers_v_lower <- sqrt(chi_sq_lower / (n * (k - 1)))
cramers_v_upper <- sqrt(chi_sq_upper / (n * (k - 1)))

cat("95% CI for Cramér's V: [", round(cramers_v_lower, 3), ",", round(cramers_v_upper, 3), "]\n")




### compare central and non-central transitions using KS test 

# perform Kolmogorov-Smirnov test
ks_test_result <- ks.test(central_probs$Probability, noncentral_probs$Probability)
print(ks_test_result)

# print minimum and maximum probabilities for central dataset
min_prob_central <- min(central_probs$Probability)
max_prob_central <- max(central_probs$Probability)
cat("Central Dataset:\n")
cat("Minimum Probability:", min_prob_central, "\n")
cat("Maximum Probability:", max_prob_central, "\n")

# print minimum and maximum probabilities for non-central dataset
min_prob_noncentral <- min(noncentral_probs$Probability)
max_prob_noncentral <- max(noncentral_probs$Probability)
cat("\nNon-Central Dataset:\n")
cat("Minimum Probability:", min_prob_noncentral, "\n")
cat("Maximum Probability:", max_prob_noncentral, "\n")

# find the row with the maximum frequency in the central dataset
most_frequent_central <- central_probs[which.max(central_probs$Frequency), ]
most_frequent_central
# find the row with the maximum frequency in the non-central dataset
most_frequent_noncentral <- noncentral_probs[which.max(noncentral_probs$Frequency), ]
most_frequent_noncentral


# Kolmogorov–Smirnov test
ks_test_result <- ks.test(central_probs$Probability, noncentral_probs$Probability)

# bootstrap to estimate Confidence Interval around D
library(boot)

ks_stat <- function(data, indices) {
  d <- data[indices]
  d1 <- d[1:(length(d)/2)]
  d2 <- d[(length(d)/2 + 1):length(d)]
  ks.test(d1, d2)$statistic
}

combined <- c(central_probs$Probability, noncentral_probs$Probability)
boot_obj <- boot(data = combined, statistic = ks_stat, R = 1000)
ci_ks <- boot.ci(boot_obj, type = "perc")

cat("KS D:", ks_test_result$statistic, "\n")
cat("KS p-value:", ks_test_result$p.value, "\n")
cat("95% CI for D:", ci_ks$percent[4], "-", ci_ks$percent[5], "\n")





