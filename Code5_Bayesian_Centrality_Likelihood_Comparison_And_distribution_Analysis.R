

library(gtools)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

##' datasets required: central_data and noncentral data with columns:
##'       Transition1 e.g. "KH"
##'       Transition2 e.g. "RC"
##'       Transition e.g. "x: KH y: RC"
##'       Frequency e.g. "1"
##'       Sequence_length e.g. "5"
##'       x_position e.g. "1"
##'       y_position e.g. "2"
##'       Model e.g. "Tiger"
##'       Distance e.g. "1"

### generate probabilities using higher order markovian analysis 

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





### dirichlet Modeling - log likelihood: does one group's data generalise to the other?

# ensure the vectors are aligned and in the same order
comparison <- merge(central_probs, noncentral_probs,
                    by = c("Transition1", "Transition2", "Distance"),
                    suffixes = c("_central", "_noncentral"))

central_vec <- comparison$Probability_central
noncentral_vec <- comparison$Probability_noncentral

epsilon <- 1e-6
central_vec <- central_vec + epsilon
noncentral_vec <- noncentral_vec + epsilon

# dirichlet log-likelihood function
dirichlet_loglik <- function(x, alpha) {
  lgamma(sum(alpha)) - sum(lgamma(alpha)) + sum((alpha - 1) * log(x))
}

# use one group's probs as alpha (scaled)
alpha_central <- central_vec * 100
alpha_noncentral <- noncentral_vec * 100

ll_central <- dirichlet_loglik(central_vec, alpha_noncentral)
ll_noncentral <- dirichlet_loglik(noncentral_vec, alpha_central)

cat("Log-likelihood (central under non-central model):", ll_central, "\n")
cat("Log-likelihood (non-central under central model):", ll_noncentral, "\n")






### dirichlet modeling - posterior distributions: how exactly do the distributions differ?

# aggregate frequencies
aggregate_central <- aggregate(Frequency ~ Transition1 + Transition2 + Distance,
                               data = central_data, FUN = sum)

aggregate_noncentral <- aggregate(Frequency ~ Transition1 + Transition2 + Distance,
                                  data = noncentral_data, FUN = sum)

# merge the two data frames
comparison <- merge(aggregate_central, aggregate_noncentral,
                    by = c("Transition1", "Transition2", "Distance"),
                    suffixes = c("_central", "_noncentral"))

# create count vectors and Dirichlet priors
central_counts <- comparison$Frequency_central
noncentral_counts <- comparison$Frequency_noncentral
alpha_prior <- rep(1, length(central_counts))

posterior_central <- central_counts + alpha_prior
posterior_noncentral <- noncentral_counts + alpha_prior

# sample from posterior Dirichlet distributions
samples_central <- rdirichlet(5000, posterior_central)
samples_noncentral <- rdirichlet(5000, posterior_noncentral)
posterior_diffs <- samples_central - samples_noncentral

# summarise posterior differences
posterior_summary <- data.frame(
  Transition = paste(comparison$Transition1, comparison$Transition2, comparison$Distance, sep = "_"),
  Mean_Diff = colMeans(posterior_diffs),
  Lower_CI = apply(posterior_diffs, 2, quantile, probs = 0.025),
  Upper_CI = apply(posterior_diffs, 2, quantile, probs = 0.975)
)

# add significance flag
posterior_summary$Significance <- ifelse(
  posterior_summary$Lower_CI > 0 | posterior_summary$Upper_CI < 0,
  "Significant", "Non-significant"
)

# save to a full transitions data frame and view
full_transitions <- posterior_summary
View(full_transitions)

write.csv(full_transitions, "full_transitions.csv", row.names = FALSE)







