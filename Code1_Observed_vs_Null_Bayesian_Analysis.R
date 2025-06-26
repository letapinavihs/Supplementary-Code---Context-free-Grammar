

#' generate mock data 
#' control for:
#'    frequency of each call
#'    frequency of each n-gram types e.g. 1-gram, 2-gram, etc.
#'    frequency of calls per n-gram type
#'    number of n-grams produced for each individual in response to each predator model
#'  Bayesian inference approach (McElreath, 2018)
mock_dt <- expand_grid(
  Individual = c("Kelly", "Elisa", "Puji", "Sina", "Yet", "Chris"),
  Model = c("Tiger", "Pattern", "Spots", "White")
) |> 
  group_by(Individual, Model) |> 
  reframe(
    Call.type = replicate(
      n = 20, # number of n-grams per response
      paste(sample(c("KS", "KH", "RC", "RO", "CC", "GR", "SI", "SH"), sample(1:5, 1)), collapse = "_")) 
  )
mock_dt



### function 2.0 to generate mock data 
#' control for:
#'  total frequency of 2-letter-codes
#'  frequencies of 2-letter-codes per n-gram type
#'  total frequency of call types
#'  total unique call types
#'  performs all the above operations^ per individual before ungrouping
#'  visualises bar graph for top 20 most (1) expected and (2) unexpected n-grams
actual_data <- read.csv("/xxx/xxx/xxx/xxx/xxx/xxx.csv") #!insert file path


# load libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# load the actual data
actual_data <- read.csv("_____.csv")

# count the frequency of each two-letter code across the whole dataset
code_frequencies_actual <- actual_data %>%
  separate_rows(Call.type, sep = "_") %>%
  count(Call.type, name = "Actual_Count")

# get  distribution of the number of two-letter codes per "Call type" (length)
actual_data <- actual_data %>%
  mutate(Length = lengths(strsplit(as.character(Call.type), "_")))

length_distribution <- actual_data %>%
  count(Length)

# count number of total two-letter codes in the actual data
total_codes_actual <- sum(actual_data$Length)

# calculate total number of Call types and unique call types
num_call_types_actual <- nrow(actual_data)  # number of rows in the actual data
num_unique_call_types_actual <- n_distinct(actual_data$Call.type)  # number of unique call types

# create placeholder for the mock Call.type data
mock_call_types <- vector("list", num_call_types_actual)

# track total number of two-letter codes generated in the mock data
total_codes_mock <- 0

# track unique call types
unique_call_types_mock <- character(0)

i <- 1
while (total_codes_mock < total_codes_actual || length(unique(unique_call_types_mock)) < num_unique_call_types_actual) {
  # sample a call type length based on the actual length distribution
  call_length <- sample(length_distribution$Length, 1, prob = length_distribution$n)
  
  # ensure we don't exceed the total two-letter codes allowed
  if (total_codes_mock + call_length > total_codes_actual) {
    call_length <- total_codes_actual - total_codes_mock  # adjust the last call type to match total codes
  }
  
  # sample the two-letter codes based on the actual frequency of codes
  codes <- sample(
    code_frequencies_actual$Call.type,
    call_length,
    replace = TRUE,
    prob = code_frequencies_actual$Actual_Count
  )
  
  # collapse them into the Call type format
  new_call_type <- paste(codes, collapse = "_")
  
  # store the generated call type only if it maintains the number of unique call types
  mock_call_types[[i]] <- new_call_type
  
  # add the generated call type to the list of unique call types
  unique_call_types_mock <- c(unique_call_types_mock, new_call_type)
  
  # update the total number of two-letter codes generated so far
  total_codes_mock <- total_codes_mock + call_length
  
  # increment the call type index
  i <- i + 1
  
  # break early if both conditions are met
  if (total_codes_mock >= total_codes_actual && length(unique(unique_call_types_mock)) >= num_unique_call_types_actual) {
    break
  }
}

# combine the mock call types into a data frame
mock_dt <- tibble(Call.type = unlist(mock_call_types))

# count the frequency of each full "Call type" in the mock dataset
call_type_frequencies_mock <- mock_dt %>%
  count(Call.type, name = "Mock_Count")

# count the frequency of each full "Call type" in the actual dataset
call_type_frequencies_actual <- actual_data %>%
  count(Call.type, name = "Actual_Count")

# merge the actual and mock frequencies into one dataset
combined_frequencies <- call_type_frequencies_actual %>%
  full_join(call_type_frequencies_mock, by = "Call.type") %>%
  mutate(
    Actual_Count = ifelse(is.na(Actual_Count), 0, Actual_Count),
    Mock_Count = ifelse(is.na(Mock_Count), 0, Mock_Count)
  )

# individual-specific call frequencies
individual_frequencies <- actual_data %>%
  separate_rows(Call.type, sep = "_") %>%
  group_by(Individual) %>%
  count(Call.type) %>%
  group_by(Individual) %>%
  mutate(Frequency = n / sum(n)) %>%
  ungroup()

# updated function definition (without the prior_alpha and prior_beta arguments)
perform_bayesian_test_individual <- function(actual_count, mock_count, prior_alpha = 1, prior_beta = 1) {
  
  total_count <- actual_count + mock_count
  
  # posterior for actual data
  posterior_alpha_actual <- prior_alpha + actual_count
  posterior_beta_actual <- prior_beta + total_count - actual_count
  
  # posterior for mock data
  posterior_alpha_mock <- prior_alpha + mock_count
  posterior_beta_mock <- prior_beta + total_count - mock_count
  
  # simulate draws from posterior distributions
  posterior_samples_actual <- rbeta(10000, posterior_alpha_actual, posterior_beta_actual)
  posterior_samples_mock <- rbeta(10000, posterior_alpha_mock, posterior_beta_mock)
  
  # probability that the call type is more frequent in actual data
  prob_actual_greater_than_mock <- mean(posterior_samples_actual > posterior_samples_mock)
  
  return(prob_actual_greater_than_mock)
}

# apply bayesian test for each call type using individual preferences
combined_frequencies <- combined_frequencies %>%
  rowwise() %>%
  mutate(Posterior_Prob = perform_bayesian_test_individual(Actual_Count, Mock_Count)) %>%
  ungroup()

# view results
print(combined_frequencies %>% arrange(desc(Posterior_Prob)))

# store table of posterior probabilities in the environment
posterior_table <- combined_frequencies

# visualise step 1: Select the top 20 and bottom 20 call types based on posterior probability
top_20 <- posterior_table %>%
  arrange(desc(Posterior_Prob)) %>%
  slice(1:20)

bottom_20 <- posterior_table %>%
  arrange(Posterior_Prob) %>%
  slice(1:20)

# combine the top and bottom call types
top_and_bottom_40 <- bind_rows(top_20, bottom_20)

# visualise step 2: add column to differentiate between top and bottom 20
top_and_bottom_40 <- top_and_bottom_40 %>%
  mutate(Group = ifelse(Call.type %in% top_20$Call.type, "Overrepresented", "Underrepresented"))

# visualise posterior probabilities with markers at the end of each bar
ggplot(top_and_bottom_40, aes(x = reorder(Call.type, Posterior_Prob), y = Posterior_Prob, color = Group)) +
  geom_bar(stat = "identity", fill = "lightgray") +
  geom_point(aes(shape = Group), size = 3, show.legend = TRUE) +
  scale_color_manual(values = c("Overrepresented" = "black", "Underrepresented" = "black")) +
  scale_shape_manual(values = c("Overrepresented" = 17, "Underrepresented" = 16)) +  # 17 is triangle, 16 is circle
  coord_flip() +
  labs(title = "Top 20 Overrepresented and Underrepresented Call Types by Posterior Probability",
       x = "Call Type", y = "Posterior Probability", color = "Group", shape = "Group") +
  theme_minimal() +
  theme(legend.position = "bottom")

# save table to a CSV file
write.csv(posterior_table, "xxxx.csv", row.names = FALSE)



