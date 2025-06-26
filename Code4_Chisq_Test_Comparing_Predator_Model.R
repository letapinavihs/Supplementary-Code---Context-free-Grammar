
###' datasets: (1) nonnatural, contains columns:
#'                    Individual = subject e.g. "Kelly"
#'                    Model = predator model e.g. "Tiger"
#'                    Recording = recording name e.g. "Kelly_tiger"
#'                    ngram = ngram length e.g. "1" 
#'                    RecordingDuration = time (s) of recording e.g. "100"
#'                    RecMaxNgram = maximum ngram length of recording e.g. "10"
#'                    RecTotalN = total number of ngrams in recording e.g. "8"
#'                    RecCallRate = RecTotalN / RecordingDuration e.g. "10"
#'                    KSKH_PREpalindorome = does the ngram have a "KS"/"KH" call pre-fix? e.g "Yes"/"No"
#'                    Palindrome = is the ngram a palindrome? "Yes"/"No"
#'             (2) lower_pal
#'                    Individual = subject e.g. "Kelly"
#'                    Model = predator model e.g. "Tiger"
#'                    Call type = ngram or single unit e.g. "KS_RC_GR"
#'                    Palindrome = is the ngram a palindrome? "Yes"/"No"

library(ggplot2)
library(dplyr)
library(ggtext)
library(tidyr)


### chi-squared tests

##chi-squared test comparing palindrome frequency across predator models

# Create a contingency table
contingency_pal <- table(nonnatural$Model, nonnatural$Palindrome)
# Perform the chi-squared test
chi_squared_pal <- chisq.test(contingency_pal)
# Print the results
print(chi_squared_pal)

##chi-squared test comparing kiss-squeak/kiss-squeak-hand prefix frequency across predator models

# Create a contingency table
contingency_prepal <- table(nonnatural$Model, nonnatural$KSKH_PREpalindrome)
# Perform the chi-squared test
chi_squared_prepal <- chisq.test(contingency_prepal)
# Print the results
print(chi_squared_prepal)



### residuals

#chi-squared test is saved as chi_squared_prepal:

# expected counts under null hypothesis
expected_counts <- chi_squared_prepal$expected
print(expected_counts)

# pearson residuals
residuals <- chi_squared_prepal$residuals
print(residuals)

# standardised residuals
std_residuals <- residuals / sqrt(expected_counts)
print(std_residuals)

# measures of association for 2x2 tables:

# phi coefficient
phi <- sqrt(chi_squared_prepal$statistic / sum(contingency_prepal))
cat("Phi coefficient:", phi, "\n")

# cramér's V (same as Phi for 2x2)
cramers_v <- sqrt(chi_squared_prepal$statistic / (sum(contingency_prepal) * (min(dim(contingency_prepal)) - 1)))
cat("Cramér's V:", cramers_v, "\n")

# from the test object
cramers_v_pal <- sqrt(chi_squared_pal$statistic / sum(contingency_pal))
cat("Cramér's V (Phi):", round(cramers_v_pal, 3), "\n")






### confidence intervals

## confidence intervals for palindrome/predator test

# convert the contingency table to a data frame
df_pal <- as.data.frame(contingency_pal)

# expand it to row-level observations
expanded_pal <- df_pal[rep(1:nrow(df_pal), df_pal$Freq), c("Var1", "Var2")]
colnames(expanded_pal) <- c("Model", "Palindrome")

library(boot)

# define the function to compute Cramér's V
cramers_v_stat <- function(data, indices) {
  d <- data[indices, ]
  tbl <- table(d$Model, d$Palindrome)
  chi <- suppressWarnings(chisq.test(tbl))
  v <- sqrt(chi$statistic / (sum(tbl) * (min(dim(tbl)) - 1)))
  return(as.numeric(v))
}

# run bootstrap
set.seed(123)
boot_obj <- boot(data = expanded_pal, statistic = cramers_v_stat, R = 1000)

# extract 95% confidence interval
ci <- boot.ci(boot_obj, type = "perc")
cat("Bootstrap 95% CI for Cramér’s V:", round(ci$percent[4], 3), "-", round(ci$percent[5], 3), "\n")



## confidence intervals for prefix/predator test

# convert the contingency table to a data frame
df_prepal <- as.data.frame(contingency_prepal)

# expand it to row-level observations
expanded_prepal <- df_prepal[rep(1:nrow(df_prepal), df_prepal$Freq), c("Var1", "Var2")]
colnames(expanded_prepal) <- c("Model", "KSKH_PREpalindrome")

# define the function to compute Cramér's V
cramers_v_stat_prepal <- function(data, indices) {
  d <- data[indices, ]
  tbl <- table(d$Model, d$KSKH_PREpalindrome)
  chi <- suppressWarnings(chisq.test(tbl))
  v <- sqrt(chi$statistic / (sum(tbl) * (min(dim(tbl)) - 1)))
  return(as.numeric(v))
}

# run bootstrap
library(boot)
set.seed(123)  # for reproducibility
boot_obj_prepal <- boot(data = expanded_prepal, statistic = cramers_v_stat_prepal, R = 1000)

# extract 95% confidence interval
ci_prepal <- boot.ci(boot_obj_prepal, type = "perc")
cat("Bootstrap 95% CI for Cramér’s V:", round(ci_prepal$percent[4], 3), "-", round(ci_prepal$percent[5], 3), "\n")

