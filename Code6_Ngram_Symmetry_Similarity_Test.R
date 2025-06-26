

## statistical similarity test
empirical_likelihood <- function(set_lengths,repeats,number_of_matches) {
  result <- c()
  for (i in 1:repeats) {
    sum <- 0
    for (j in 1:length(set_lengths)) {
      length <- set_lengths[j]
      permutation <- sample(1:length,length,replace=FALSE)
      matches <- sum(1:length==permutation)
      sum <- sum+matches
    }
    result <- c(result,sum)
  }
  expected <- sum(number_of_matches<=result)
  p_value <- expected/repeats
  statement <- paste("The approximate p-value for seeing at least",number_of_matches,"matches is:",p_value)
  print(statement)
  hist(result,xlim=range(0,sum(set_lengths)))
  abline(v=number_of_matches,col="red")
}
# Run the function=> "The approximate p-value for seeing at least 23 matches is: 0.000746"
empirical_likelihood(c(5,2,3,5,2,3,8,4,9,10,2),1000000
                     ,23)



## ^statistical similarity test with bootstrapping

library(boot)

empirical_likelihood <- function(set_lengths, repeats, number_of_matches) {
  # Vector to store permutation sums (null distribution)
  result <- numeric(repeats)
  
  for (i in 1:repeats) {
    sum_matches <- 0
    for (j in seq_along(set_lengths)) {
      length <- set_lengths[j]
      permutation <- sample(1:length, length, replace = FALSE)
      matches <- sum(1:length == permutation)
      sum_matches <- sum_matches + matches
    }
    result[i] <- sum_matches
  }
  
  # calculate mean and SD of null distribution
  mean_null <- mean(result)
  sd_null <- sd(result)
  
  # calculate effect size (Cohen's d)
  effect_size <- (number_of_matches - mean_null) / sd_null
  
  # bootstrap function to calculate effect size for bootstrap samples
  boot_fun <- function(data, indices) {
    d <- data[indices]
    (number_of_matches - mean(d)) / sd(d)
  }
  
  # bootstrap 1000 samples to estimate 95% confidence interval for effect size
  boot_results <- boot(result, boot_fun, R = 1000)
  ci <- boot.ci(boot_results, type = "perc")$percent[4:5]  # percentile CI lower and upper
  
  # calculate exact p-value from permutation distribution
  p_value <- mean(result >= number_of_matches)
  
  # view results
  statement <- sprintf(
    "Observed matches = %d; mean under null = %.2f (SD = %.2f); effect size d = %.2f, 95%% CI [%.2f, %.2f]; p = %.6f (permutation test with %d repeats).",
    number_of_matches, mean_null, sd_null, effect_size, ci[1], ci[2], p_value, repeats
  )
  print(statement)
  
  # plot histogram with observed value marked
  hist(result, xlim = range(0, sum(set_lengths)), main = "Permutation distribution of matches", xlab = "Number of matches")
  abline(v = number_of_matches, col = "red")
  
  # return results as a list (optional, for further use)
  return(list(
    observed = number_of_matches,
    null_distribution = result,
    mean_null = mean_null,
    sd_null = sd_null,
    effect_size = effect_size,
    effect_size_ci = ci,
    p_value = p_value
  ))
}

res <- empirical_likelihood(c(5,2,3,5,2,3,8,4,9,10,2), 1000000, 23)


