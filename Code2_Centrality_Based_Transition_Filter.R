# ------------------------------------------------------------------
# This code defines the functions used to calculate p-values for each transition 
# between call units under different conditions (the within condition was used for this study). 
# The p-values reflect how observed transition proportions differ from what would be 
# expected by chance, based on randomised simulations under the null hypothesis.
# ------------------------------------------------------------------

require(rlang)
require(cli)
require(dplyr)
require(tidyr)
require(tibble)
require(pbapply)
require(purrr)
require(dplyr)

data <- read.csv('xxx.csv')

#' Retrieves frequencies of transitions between unit-elements of call
#' sequences.Transitions can be evaluated within the same sequence or between
#' neighbouring sequences. 
#' 
#' There is no overlap between responses in that sequences are grouped by 
#' their unique individual-model combinations, e.g. a unit from Kelly-Spots will 
#' not be analysed in the same transition pair as a unit from Kelly-Tiger
#' 
#' Transitions describe shifts between call-units in a given position of the
#' 'reference sequence' to a given position in the 'contrast sequence'. For
#' within-sequence transitions, reference and contrast sequences are identical
#' (`direction = "within"`). On between-sequence transitions, the contrast
#' sequence can be either the sequence immediately preceding (`direction = "preceding"`) 
#' or succeeding (`direction = "succeeding"`) the reference sequence.
#'  
#'  Function allows the evaluation of multiple pairs of outbound-inbound
#'  positions in sequences under contrast (parameters `from_pos` and `to_pos`,
#'  which must be of equal length).
#' 
#' @param data a data frame, comprising vocal calls ordered chronologically
#' @param call_col <data-masking> the name of the column in `data` containing
#'   sequences of vocal calls
#'!!! @param from_length a integer vector, specifying the desired sequence length of 
#'    the reference sequence. for example, if `to_length = c(1, 4)` then only 
#'    from/reference sequence lengths of 1 and 4 are extracted. Results are 
#'    aggregated over the chosen sequence lengths.
#'!!! @param to_length a integer vector, specifying the desired sequence length of 
#'    the contrast sequence. for example, if `to_length = c(2, 3)` then only 
#'    to/contrast sequence lengths of 2 and 3 are extracted. Results are 
#'    aggregated over the chosen sequence lengths.
#' @param direction character string, defining the type of contrast. For each
#'   given sequence in the input data, three options of evaluation:
#'    - "within": transitions within the same sequence (default).
#'    - "preceding": transitions to the immediately prior sequence
#'    - "succeeding": transitions to the immediately posterior sequence
#' @param from_pos an integer vector, defining the outbound unit position(s)
#'   in the reference sequence under evaluation
#' @param to_pos an integer vector, defining the inbound unit position(s) in
#'   the contrast sequence under evaluation#'   
#' @param random_run logical, whether to randomize the input data. If `TRUE`,
#'   sequence order is shuffled, and call units are permuted across all call
#'   sequences. In other words, observed call units are randomly re-allocated to
#'   sequences regardless of the temporal and positional structure underlying
#'   the original input dataset.#'   
#' @param call_unit_codes string vector, specifying the expected call codes of
#'   call unit-elements. Function will return an error if the input data
#'   includes a code not comprised in this vector
#' @param by_pairwise_pos logical, whether to output results by pairwise
#'   outbound-inbound positions. If `TRUE`, outputs are provided as lists, one
#'   for each pairwise position 
#' @param output_type character string, specifying the type of output. Four options:
#'    - "matrix": the transition matrix, i.e. a squared matrix providing the number 
#'    of transitions between each possible pair of unit-elements in call sequences.
#'    - "trans_dets": data frame with the detailed information on each of the 
#'    retrieved transitions
#'    - "trans_freqs": similar to option "matrix", but provided in long format 
#'    and comprising only non-zero transitions.
#'!!!    - "prop_matrix": a matrix containing the proportions of transitions 
#'    relative to the whole matrix
#'    - "all": a list containing all the above objects
#'    
#'
#'
#'
get_transition_matrix_expanded <- function(data, 
                                           call_col,
                                           from_length = NULL,
                                           to_length = NULL,
                                           direction = c("within", "preceding", "succeeding"),
                                           from_pos, 
                                           to_pos,
                                           random_run = FALSE,
                                           call_unit_codes = c('CC', 'GR', 'RC', 'KS', 'SH', 'SI', 'RO', 'KH'),
                                           by_pairwise_pos = FALSE,
                                           output_type = c("matrix", "trans_dets", "trans_freqs", "prop_matrix", "within"),
                                           individuals = NULL,
                                           models = NULL
){
  # filter by specified individuals and models if provided
  if (!is.null(individuals)) {
    data <- data |> dplyr::filter(Individual %in% individuals)
  }
  
  if (!is.null(models)) {
    data <- data |> dplyr::filter(Model %in% models)
  }
  
  # check if there is data left after filtering
  if (nrow(data) == 0) {
    stop("There is no data for the specified individual-model combination.")
  }
  
  # data pre-processing ----------------------------------------------
  data <- data |> 
    dplyr::as_tibble() |> 
    dplyr::rename(call_seq = {{call_col}}) |> 
    tibble::add_column(call_id = 1:nrow(data), .before = 1) |> 
    dplyr::mutate(call_seq = strsplit(call_seq, "_"))
  
  # group by Individual and Model to remove between unique individual-model transition extraction
  data <- data |> 
    dplyr::group_by(Individual, Model)
  
  # inputs validation -----------------------------------------------
  
  # validate optional arguments
  direction <- rlang::arg_match(direction)
  output_type <- rlang::arg_match(output_type)
  
  # check consistency in unit codes between data and those expected under `call_unit_codes`
  unexpect <- data |> 
    pull(call_seq) |> 
    unlist() |> 
    unique() |> 
    setdiff(call_unit_codes)
  
  if(length(unexpect) > 1){
    cli::cli_abort(
      c(
        "Call unit codes in input data must be elements of vector specified in argument {.code call_unit_codes}",
        "x" = "{cli::qty(unexpect)} Unit code{?s} {.str {unexpect}} {?is/are} not in {cli::qty(call_unit_codes)} ({.str {call_unit_codes}})"
      )
    )
  }
  
  # check if specified from-to vectors have equal length
  if(length(from_pos) != length(to_pos)){
    cli::cli_abort("Input parameters {.code from_pos} and {.code to_pos} must be integer vectors of equal length")
  }
  
  # randomisation step ------------------------------------------------
  if(random_run){
    
    # shuffle order of sequences
    data$call_seq <- sample(data$call_seq)
    
    # shuffle unit codes within and across sequences 
    data <- data |> 
      # expand to one row per code
      tidyr::unnest_longer(call_seq) |> 
      # shuffle code positions
      dplyr::mutate(call_seq = sample(call_seq)) |> 
      # compress to original format
      dplyr::group_by(across(-call_seq)) |> 
      dplyr::summarise(call_seq = list(call_seq), .groups = "drop")
  }
  
  # derive transitions -------------------------------------------------
  
  # set up data based on choice of direction, by binding sequence to compare to
  data <- data |>
    dplyr::mutate(
      call_seq_contrast = dplyr::case_when(
        direction == "within" ~ call_seq,
        direction == "preceding" ~ dplyr::lag(call_seq),
        direction == "succeeding" ~ dplyr::lead(call_seq)
      ))
  
  # calculate transitions, for each specified pair of from-to positions
  trans <- purrr::map2(from_pos, to_pos, \(from, to){
    
    data |> 
      mutate(
        # add columns with lengths of reference and contrast sequences
        dplyr::across(c(call_seq, call_seq_contrast), ~purrr::map_int(.x, length), .names = "{.col}_lt"),
        from_to = paste0(from, "-", to),
        # extract unit codes from sought positions
        unit_from = purrr::map(call_seq, from),
        unit_to = purrr::map(call_seq_contrast, to)
      ) |> 
      tidyr::unnest(c(unit_from, unit_to), keep_empty = TRUE) |> 
      tidyr::drop_na()
  })
  
  # filter to specified reference sequence lengths, if defined
  if(!is.null(from_length)){
    trans <- trans |> 
      purrr::map(~dplyr::filter(.x, call_seq_lt %in% from_length))
  }
  if(!is.null(to_length)){
    trans <- trans |> 
      purrr::map(~dplyr::filter(.x, call_seq_contrast_lt %in% to_length))
  }
  
  # transition frequencies
  trans_freqs <- trans |> purrr::map(~dplyr::count(.x, from_to, unit_from, unit_to))
  
  # name list elements as the pairwise positions under evaluation
  names(trans_freqs) <- purrr::map2_chr(from_pos, to_pos, ~paste0(.x, "-", .y))
  
  # coerce as transition matrix ----------------------------------------
  
  # initialise template transition matrix
  m_size <- length(call_unit_codes)
  temp_trans_mtx <- matrix(0, nrow = m_size, ncol = m_size)
  rownames(temp_trans_mtx) <- call_unit_codes
  colnames(temp_trans_mtx) <- call_unit_codes
  
  # generate transition matrices, for each from-to pair
  trans_mtx <- purrr::map(
    trans_freqs,
    \(x){
      # set-up matrix for current from-to pair
      c_mtx <- temp_trans_mtx
      # populate matrix
      x |>
        purrr::pwalk(\(unit_from, unit_to, n, ...){
          c_mtx[unit_from, unit_to] <<- c_mtx[unit_from, unit_to] + n
        })
      return(c_mtx)
    })
  
  # sum all individual transition matrices into one matrix
  trans_mtx <- Reduce(`+`, trans_mtx)
  
  # calculate proportion matrix ----------------------------------------
  
  # calculate the total number of transitions
  total_transitions <- sum(trans_mtx)
  
  # compute proportion matrix
  if (total_transitions > 0) {
    prop_matrix <- trans_mtx / total_transitions
  } else {
    prop_matrix <- trans_mtx  # If there are no transitions, keep the matrix unchanged
  }
  
  # return results ---------------------------------------------------
  
  # aggregate transition frequencies over pairwise positions, if option activated
  if(!by_pairwise_pos){
    trans_freqs <- trans_freqs |> 
      purrr::list_rbind() |> 
      dplyr::ungroup() |>  # Ungroup before summarizing
      dplyr::summarise(n = sum(n), .by = c(unit_from, unit_to))
  }
  
  # output given chosen option
  switch(
    output_type,
    matrix = trans_mtx,
    trans_dets = trans,
    trans_freqs = trans_freqs,
    prop_matrix = prop_matrix,
    all = list(trans_dets = trans, trans_freqs = trans_freqs, trans_mtx = trans_mtx, prop_matrix = prop_matrix)
  )
}





#' ////////////////////////////////////////////////////////////////////////////////
#' Wrapper to generate mean and SD for each transition cell under the null
#' hypothesis of no relationship between neighbouring call units
#' Added generation of proportional matrix derived from mean matrix 
#' Inputs include from_length and to_length to specify reference sequence 
#'  length and contrast sequence length respectively
#'  
#'  
h0_trans_mean_sd <- function(N = 100, data, call_col, from_length, to_length, direction,
                             from_pos, to_pos, by_pairwise_pos) {
  
  # get random replicates under H0
  transition_reps <- pbapply::pbreplicate(
    n = N,
    get_transition_matrix_expanded(
      data = data,
      call_col = !!enquo(call_col),
      from_length = from_length,
      to_length = to_length,
      direction = direction,
      from_pos = from_pos,
      to_pos = to_pos,
      random_run = TRUE,
      by_pairwise_pos = by_pairwise_pos,
      output_type = "matrix"
    ),
    simplify = "array"
  )
  
  # check the structure of transition_reps
  if (length(dim(transition_reps)) != 3) {
    stop("transition_reps does not have the expected 3 dimensions.")
  }
  
  # calculate expected frequencies and associated SDs for each cell of the transition matrix
  if (by_pairwise_pos) {
    
    result_mean <- result_sd <- list()
    trans_dim <- dim(transition_reps)[1:2]
    trans_names <- dimnames(transition_reps[,,1])
    
    for (i in seq_len(dim(transition_reps)[3])) {
      
      trans_label <- rownames(transition_reps)[i]
      
      trans_mtx_rep <- array(
        unlist(transition_reps[,,i]),
        dim = trans_dim,
        dimnames = trans_names
      )
      
      result_mean[[trans_label]] <- apply(trans_mtx_rep, c(1,2), mean)
      result_sd[[trans_label]] <- apply(trans_mtx_rep, c(1,2), sd)
    }
    
  } else {
    result_mean <- apply(transition_reps, c(1,2), mean)
    result_sd <- apply(transition_reps, c(1,2), sd)
  }
  
  # generate random proportion matrix based on the result_mean
  generate_random_prop_mtx <- function(h0_result_mean, num_transitions = 1000) {
    # convert the result mean into a probability matrix
    prob_matrix <- h0_result_mean / sum(h0_result_mean)
    
    # flatten the probability matrix for sampling
    unit_codes <- rownames(prob_matrix)
    flattened_probs <- as.vector(prob_matrix)
    names(flattened_probs) <- as.vector(outer(unit_codes, unit_codes, paste, sep = "-"))
    
    # simulate transitions based on the probability matrix
    simulated_transitions <- sample(
      names(flattened_probs),
      size = num_transitions,
      replace = TRUE,
      prob = flattened_probs
    )
    
    # create an empty transition matrix
    random_trans_mtx <- matrix(0, nrow = length(unit_codes), ncol = length(unit_codes))
    rownames(random_trans_mtx) <- unit_codes
    colnames(random_trans_mtx) <- unit_codes
    
    # fill the transition matrix based on simulated transitions
    for (transition in simulated_transitions) {
      from_to <- strsplit(transition, "-")[[1]]
      random_trans_mtx[from_to[1], from_to[2]] <- random_trans_mtx[from_to[1], from_to[2]] + 1
    }
    
    # convert the transition matrix into a proportion matrix
    total_simulated_transitions <- sum(random_trans_mtx)
    if (total_simulated_transitions > 0) {
      random_prop_mtx <- random_trans_mtx / total_simulated_transitions
    } else {
      random_prop_mtx <- random_trans_mtx  # If no transitions, keep as zero matrix
    }
    
    return(random_prop_mtx)
  }
  
  # generate the random proportion matrix
  if (by_pairwise_pos) {
    random_prop_mtx <- purrr::map(result_mean, ~generate_random_prop_mtx(.x))
  } else {
    random_prop_mtx <- generate_random_prop_mtx(result_mean)
  }
  
  # output
  list(
    result_mean = result_mean,
    result_sd = result_sd,
    random_prop_mtx = random_prop_mtx
  )
}






#' ///////////////////////////////////////////////////////////////////////////////
#' 
#'perm_test integrates both get_transition_matrix_expanded and h0_trans_mean_sd 
#' so that both sets of outputs can be generated from one input.
#' 
#' #' It then:  
#'  (1) compares prop_matrix (the actual transitions) with 
#'      random_prop_matrix (the null hypothesis)
#'  (2) outputs p_values and test_stats in the form of matrices that correspond 
#'      with the rows and columns of the prop_matrix and random_prop_matrix 
#'      transitions
#'      
perm_test <- function(data, 
                      call_col,
                      from_length,
                      to_length,
                      direction,
                      from_pos,
                      to_pos,
                      random_run = FALSE,
                      output_type = c("matrix", "trans_dets", "trans_freqs", "prop_matrix", "all"),
                      individuals = NULL,
                      models = NULL,
                      N = 100,
                      by_pairwise_pos = FALSE) {
  
  # validate that shared arguments match between functions
  validate_inputs <- function(x1, x2, name) {
    if (!identical(x1, x2)) {
      stop(paste("Inputs do not match for", name))
    }
  }
  
  # validate the shared inputs for both functions
  validate_inputs(from_length, from_length, "from_length")
  validate_inputs(to_length, to_length, "to_length")
  validate_inputs(direction, direction, "direction")
  validate_inputs(from_pos, from_pos, "from_pos")
  validate_inputs(to_pos, to_pos, "to_pos")
  validate_inputs(by_pairwise_pos, by_pairwise_pos, "by_pairwise_pos")
  
  # run get_transition_matrix_expanded
  transition_matrix_results <- get_transition_matrix_expanded(
    data = data,
    call_col = call_col,
    from_length = from_length,
    to_length = to_length,
    direction = direction,
    from_pos = from_pos,
    to_pos = to_pos,
    random_run = random_run,
    call_unit_codes = c('CC', 'GR', 'RC', 'KS', 'SH', 'SI', 'RO', 'KH'),
    by_pairwise_pos = by_pairwise_pos,
    output_type = output_type,
    individuals = individuals,
    models = models
  )
  
  # run h0_trans_mean_sd
  h0_results <- h0_trans_mean_sd(
    N = N,
    data = data,
    call_col = call_col,
    from_length = from_length,
    to_length = to_length,
    direction = direction,
    from_pos = from_pos,
    to_pos = to_pos,
    by_pairwise_pos = by_pairwise_pos
  )
  
  # calculate the observed proportions and the test statistic
  observed_prop <- transition_matrix_results$prop_matrix
  difference_matrix <- observed_prop - h0_results$random_prop_mtx
  test_stat <- abs(difference_matrix)
  
  # generate signed test statistic
  sign_matrix <- sign(difference_matrix)
  signed_test_stat <- sign_matrix * test_stat
  
  # initialize matrix for null distribution of test statistics
  null_dist <- array(0, dim = c(N, nrow(test_stat), ncol(test_stat)))
  
  # generate null distribution
  for (i in 1:N) {
    # Shuffle data
    shuffled_data <- data[sample(nrow(data)), ]
    
    # recalculate transition matrix with shuffled data
    shuffled_results <- get_transition_matrix_expanded(
      data = shuffled_data,
      call_col = call_col,
      from_length = from_length,
      to_length = to_length,
      direction = direction,
      from_pos = from_pos,
      to_pos = to_pos,
      random_run = TRUE,
      call_unit_codes = c('CC', 'GR', 'RC', 'KS', 'SH', 'SI', 'RO', 'KH'),
      by_pairwise_pos = by_pairwise_pos,
      output_type = output_type,
      individuals = individuals,
      models = models
    )
    
    # calculate test statistic for shuffled data
    shuffled_prop <- shuffled_results$prop_matrix
    null_dist[i, , ] <- abs(shuffled_prop - h0_results$random_prop_mtx)
  }
  
  # calculate p-values
  p_values <- apply(test_stat, c(1, 2), function(obs_stat) {
    mean(null_dist >= obs_stat)
  })
  
  # extract necessary data to create trans_dets_updated
  trans_freqs_df <- transition_matrix_results$trans_freqs
  prop_matrix_df <- as.data.frame(as.table(observed_prop))
  stat_matrix_df <- as.data.frame(as.table(signed_test_stat))
  pval_matrix_df <- as.data.frame(as.table(p_values))
  
  # merge all extracted data
  trans_dets_updated <- trans_freqs_df %>%
    left_join(prop_matrix_df, by = c("unit_from" = "Var1", "unit_to" = "Var2")) %>%
    left_join(stat_matrix_df, by = c("unit_from" = "Var1", "unit_to" = "Var2")) %>%
    left_join(pval_matrix_df, by = c("unit_from" = "Var1", "unit_to" = "Var2")) %>%
    dplyr::rename(proportion = Freq.x, stat = Freq.y, pval = Freq) %>%
    dplyr::mutate(
      from_to = paste0(unit_from, "-", unit_to),
      Individual = paste(individuals, collapse = ", "),
      Model = paste(models, collapse = ", "),
      from_length = paste(from_length, collapse = ", "),
      to_length = paste(to_length, collapse = ", "),
      from_pos = paste(from_pos, collapse = ", "),
      to_pos = paste(to_pos, collapse = ", ")
    ) %>%
    dplyr::select(from_to, unit_from, unit_to, n, proportion, stat, pval, Individual, Model, from_length, to_length, from_pos, to_pos)
  
  # store the updated dataframe in the global environment
  trans_dets_updated_global <<- trans_dets_updated
  
  # combine and return all results including the signed test statistic and trans_dets_updated
  list(
    trans_mtx = transition_matrix_results$trans_mtx,
    trans_dets = transition_matrix_results$trans_dets,
    trans_freqs = transition_matrix_results$trans_freqs,
    prop_matrix = transition_matrix_results$prop_matrix,
    result_mean = h0_results$result_mean,
    result_sd = h0_results$result_sd,
    random_prop_mtx = h0_results$random_prop_mtx,
    test_stat = test_stat,
    signed_test_stat = signed_test_stat,
    p_values = p_values,
    trans_dets_updated = trans_dets_updated
  )
}




#' ////////////////////////////////////////////////////////////////////////////////
#'
#'Run all the function that extracts the transitions:
#'
get_transition_matrix_expanded(
  data = data,
  call_col = Call.type,
  from_length = c(5),
  to_length = c(5),
  direction = "succeeding",
  from_pos = c(3),
  to_pos = c(3),
  random_run = FALSE,
  by_pairwise_pos = FALSE,
  output_type = "within",
  models = c("Tiger")
)

#'Run all the function that creates the null transitions:
#'
h0_trans_mean_sd(
  N = 100,
  data = data,
  call_col = Call.type,
  from_length = c(5),
  to_length = c(5),
  direction = "succeeding",
  from_pos = c(3),
  to_pos = c(3),
  by_pairwise_pos = FALSE
)

#'Run the combined function that combines both previous steps in one:
#'
perm_test(
  data = data,
  call_col = "Call.type",
  from_length = c(5),
  to_length = c(5),
  direction = "succeeding",
  from_pos = c(3),
  to_pos = c(3),
  random_run = FALSE, 
  output_type = "within",
  individuals = c("Kelly ", "Kelly", "Elisa", "Chris", "Sina", "Yet", "Puji"),
  models = c("Tiger"),
  N = 100,
  by_pairwise_pos = FALSE
)



# ------------------------------------------------------------------
# This code performs n-gram transition analysis for vocal sequences.
# It supports a statistical investigation of center-embedding in sequences,
# as predicted for palindromic structures in animal communication.
# 
# 1. It calculates frequencies of unique "Call type" sequences from CSV input.
# 2. It filters significant transitions using a custom R-based model.
# 3. It visualises transition patterns in network graphs across various 
#    positions relative to the center of sequences (e.g., center ±1, ±2, etc.).
# 
# These analyses are designed to test for symmetrical transition biases—
# specifically, whether transitions are more likely to stem from the central unit
# of an n-gram (e.g., "C" in A-B-C-D-E) versus non-central units.
#
# Outputs include frequency tables and position-specific network graphs,
# allowing assessment of structure within communicative sequences.
# ------------------------------------------------------------------


# read CSV file containing call transition data
data <- read.csv("xxx.csv")
# check structure of the data
str(data)
# calculate frequencies of unique sequences in the "Call type" column
sequence_freq <- table(data$Call.type)
# convert table to a data frame for better manipulation (optional)
sequence_freq_df <- as.data.frame(sequence_freq)
# rename the columns for clarity (optional)
colnames(sequence_freq_df) <- c("Call type", "Frequency")
# print the resulting data frame
print(sequence_freq_df)
# optionally, sort the data frame by frequency
sorted_sequence_freq_df <- sequence_freq_df[order(sequence_freq_df$Frequency, decreasing = TRUE), ]
print(sorted_sequence_freq_df)


# load required library
library(openxlsx)
# define file path for the Excel file containing call transition data
excel_file <- "xxx.xlsx"
# write data frame to an Excel file
write.xlsx(sequence_freq_df, excel_file)
# confirm the file has been created
if (file.exists(excel_file)) {
  cat("Excel file exported successfully.\n")
} else {
  cat("Export failed.\n")
}







library(igraph)
library(dplyr)
library(scales) 


### Network Map: any sequence length and model

generate_network_maps <- function(model, sequence_length, csv_file) {
  # read the CSV file
  data <- read.csv(csv_file)
  # filter the data based on criteria
  filtered_data <- data %>%
    filter(Sig == "Yes", 
           Sequence_length == sequence_length,
           Model == model)
  # determine the maximum and minimum frequency across the entire dataset
  max_freq <- max(filtered_data$Frequency)
  min_freq <- min(filtered_data$Frequency)
  # iterate over unique pairs of x_position and y_position
  unique_positions <- unique(paste(filtered_data$x_position, filtered_data$y_position))
  for (pos in unique_positions) {
    # subset data for the current position
    subset_data <- filtered_data %>% filter(paste(x_position, y_position) == pos)
    # create a graph object
    g <- graph.data.frame(subset_data[, c("Transition1", "Transition2")], directed = TRUE)
    # set edge attributes
    edge_attrs <- subset_data$Frequency
    names(edge_attrs) <- paste(subset_data$Transition1, subset_data$Transition2, sep = "_")
    g <- set_edge_attr(g, "frequency", value = edge_attrs)
    # rescale edge weights based on maximum and minimum frequency
    scaled_weights <- rescale(edge_attrs, to = c(0.3, 5), from = c(min_freq, max_freq))
    # get x_position, y_position, and sequence length
    position_title <- paste("x_position:", unique(subset_data$x_position), "y_position:", unique(subset_data$y_position))
    sequence_title <- paste("Sequence_length:", sequence_length)
    # concatenate model name, sequence length, and position title using comma
    title <- paste("Model:", model, ",", sequence_title, ",", position_title)
    # plot the graph with title
    plot(g, layout = layout_in_circle(g), 
         main = title,
         edge.curved = FALSE, edge.arrow.size = 0.4, 
         edge.color = "black", edge.width = scaled_weights, 
         vertex.size = 40, vertex.label = V(g)$name, 
         vertex.label = TRUE)
  }
}
# usage example
generate_network_maps(model = "Tiger", sequence_length = 5, csv_file = "xxx.csv")







### NETWORK MAP: CENTRAL +/- 1

generate_network_map <- function(model, csv_file) {
  # read the CSV file
  data <- read.csv(csv_file)
  # filter the data based on criteria
  filtered_data <- data %>%
    filter(Sig == "Yes", 
           Model == model,
           ((Sequence_length == 3 & x_position == 2 & y_position %in% c(1, 3)) |
              (Sequence_length == 4 & x_position == "2_3" & y_position %in% c(1, 4)) |
              (Sequence_length == 5 & x_position == 3 & y_position %in% c(2, 4)) |
              (Sequence_length == 6 & x_position == "3_4" & y_position %in% c(2, 5)) |
              (Sequence_length == 7 & x_position == 4 & y_position %in% c(3, 5)) |
              (Sequence_length == 8 & x_position == "4_5" & y_position %in% c(3, 6)) |
              (Sequence_length == 9 & x_position == 5 & y_position %in% c(4, 6)) |
              (Sequence_length == 10 & x_position == "5_6" & y_position %in% c(4, 7))))
  # check if there is any filtered data
  if (nrow(filtered_data) == 0) {
    cat("No data matching the provided criteria.\n")
    return(NULL)
  }
  # determine the maximum and minimum frequency across the entire dataset
  max_freq <- max(filtered_data$Frequency, na.rm = TRUE)
  min_freq <- min(filtered_data$Frequency, na.rm = TRUE)
  # create a graph object
  g <- graph.data.frame(filtered_data[, c("Transition1", "Transition2")], directed = TRUE)
  # set edge attributes
  edge_attrs <- filtered_data$Frequency
  names(edge_attrs) <- paste(filtered_data$Transition1, filtered_data$Transition2, sep = "_")
  g <- set_edge_attr(g, "frequency", value = edge_attrs)
  # rescale edge weights based on maximum and minimum frequency
  scaled_weights <- rescale(edge_attrs, to = c(0.3, 5), from = c(min_freq, max_freq))
  # calculate Fruchterman-Reingold layout
  layout <- layout_with_fr(g)
  # get model name as a string
  model_str <- paste("Model:", model)
  # concatenate model name in the title
  title <- model_str
  # define colors for Transition1 and Transition2 nodes
  node_colors <- ifelse(V(g)$name %in% filtered_data$Transition1, "darkorange", "gold")
  # plot the graph with title and labeled nodes
  plot(g, layout = layout, 
       main = title,
       edge.curved = FALSE, edge.arrow.size = 0.4, 
       edge.color = "black", edge.width = scaled_weights, 
       vertex.size = 13, vertex.label.dist = 0, 
       vertex.color = node_colors,
       vertex.label = V(g)$name, 
       vertex.label.color = "black",
       vertex.label.font = 2)
}

# usage:
generate_network_map(model = "Tiger", 
                     csv_file = "xxx.csv")
generate_network_map(model = "Pattern", 
                     csv_file = "xxx.csv")
generate_network_map(model = "Spots", 
                     csv_file = "xxx.csv")
generate_network_map(model = "White", 
                     csv_file = "xxx.csv")




### NETWORK MAP: CENTRAL +/- 2

generate_network_map <- function(model, csv_file) {
  # read the CSV file
  data <- read.csv(csv_file)
  # filter the data based on criteria
  filtered_data <- data %>%
    filter(Sig == "Yes", 
           Model == model,
           ((Sequence_length == 5 & x_position == 3 & y_position %in% c(1, 5)) |
              (Sequence_length == 6 & x_position == "3_4" & y_position %in% c(1, 6)) |
              (Sequence_length == 7 & x_position == 4 & y_position %in% c(2, 6)) |
              (Sequence_length == 8 & x_position == "4_5" & y_position %in% c(2, 7)) |
              (Sequence_length == 9 & x_position == 5 & y_position %in% c(3, 7)) |
              (Sequence_length == 10 & x_position == "5_6" & y_position %in% c(3, 8))))
  # check if there is any filtered data
  if (nrow(filtered_data) == 0) {
    cat("No data matching the provided criteria.\n")
    return(NULL)
  }
  # determine the maximum and minimum frequency across the entire dataset
  max_freq <- max(filtered_data$Frequency, na.rm = TRUE)
  min_freq <- min(filtered_data$Frequency, na.rm = TRUE)
  # create a graph object
  g <- graph.data.frame(filtered_data[, c("Transition1", "Transition2")], directed = TRUE)
  # set edge attributes
  edge_attrs <- filtered_data$Frequency
  names(edge_attrs) <- paste(filtered_data$Transition1, filtered_data$Transition2, sep = "_")
  g <- set_edge_attr(g, "frequency", value = edge_attrs)
  # rescale edge weights based on maximum and minimum frequency
  scaled_weights <- rescale(edge_attrs, to = c(0.3, 5), from = c(min_freq, max_freq))
  # Calculate Fruchterman-Reingold layout
  layout <- layout_with_fr(g)
  # get model name as a string
  model_str <- paste("Model:", model)
  # concatenate model name in the title
  title <- model_str
  # define colors for Transition1 and Transition2 nodes
  node_colors <- ifelse(V(g)$name %in% filtered_data$Transition1, "darkorange", "gold")
  # plot the graph with title and labeled nodes
  plot(g, layout = layout, 
       main = title,
       edge.curved = FALSE, edge.arrow.size = 0.4, 
       edge.color = "black", edge.width = scaled_weights, 
       vertex.size = 13, vertex.label.dist = 0, 
       vertex.color = node_colors,
       vertex.label = V(g)$name, 
       vertex.label.color = "black",
       vertex.label.font = 2)
}

# usage:
generate_network_map(model = "Tiger", 
                     csv_file = "xxx.csv")
generate_network_map(model = "Pattern", 
                     csv_file = "xxx.csv")
generate_network_map(model = "Spots", 
                     csv_file = "xxx.csv")
generate_network_map(model = "White", 
                     csv_file = "xxx.csv")




# NETWORK MAP: CENTRAL +/- 3

generate_network_map <- function(model, csv_file) {
  # read the CSV file
  data <- read.csv(csv_file)
  # filter the data based on criteria
  filtered_data <- data %>%
    filter(Sig == "Yes", 
           Model == model,
           ((Sequence_length == 7 & x_position == 4 & y_position %in% c(1, 7)) |
              (Sequence_length == 8 & x_position == "4_5" & y_position %in% c(1, 8)) |
              (Sequence_length == 9 & x_position == 5 & y_position %in% c(2, 8)) |
              (Sequence_length == 10 & x_position == "5_6" & y_position %in% c(2, 9))))
  # check if there is any filtered data
  if (nrow(filtered_data) == 0) {
    cat("No data matching the provided criteria.\n")
    return(NULL)
  }
  # determine the maximum and minimum frequency across the entire dataset
  max_freq <- max(filtered_data$Frequency, na.rm = TRUE)
  min_freq <- min(filtered_data$Frequency, na.rm = TRUE)
  # create a graph object
  g <- graph.data.frame(filtered_data[, c("Transition1", "Transition2")], directed = TRUE)
  # set edge attributes
  edge_attrs <- filtered_data$Frequency
  names(edge_attrs) <- paste(filtered_data$Transition1, filtered_data$Transition2, sep = "_")
  g <- set_edge_attr(g, "frequency", value = edge_attrs)
  # rescale edge weights based on maximum and minimum frequency
  scaled_weights <- rescale(edge_attrs, to = c(0.3, 5), from = c(min_freq, max_freq))
  # calculate Fruchterman-Reingold layout
  layout <- layout_with_fr(g)
  # get model name as a string
  model_str <- paste("Model:", model)
  # concatenate model name in the title
  title <- model_str
  # define colors for Transition1 and Transition2 nodes
  node_colors <- ifelse(V(g)$name %in% filtered_data$Transition1, "darkorange", "gold")
  # plot the graph with title and labeled nodes
  plot(g, layout = layout, 
       main = title,
       edge.curved = FALSE, edge.arrow.size = 0.4, 
       edge.color = "black", edge.width = scaled_weights, 
       vertex.size = 13, vertex.label.dist = 0, 
       vertex.color = node_colors,
       vertex.label = V(g)$name, 
       vertex.label.color = "black",
       vertex.label.font = 2)
}

# usage:
generate_network_map(model = "Tiger", 
                     csv_file = "xxx.csv")
generate_network_map(model = "Pattern", 
                     csv_file = "xxx.csv")
generate_network_map(model = "Spots", 
                     csv_file = "xxx.csv")
generattie_network_map(model = "White", 
                       csv_file = "xxx.csv")


# NETWORK MAP: CENTRAL +/- 4

generate_network_map <- function(model, csv_file) {
  # read the CSV file
  data <- read.csv(csv_file)
  # filter the data based on criteria
  filtered_data <- data %>%
    filter(Sig == "Yes", 
           Model == model,
           ((Sequence_length == 9 & x_position == 5 & y_position %in% c(1, 9)) |
              (Sequence_length == 10 & x_position == "5_6" & y_position %in% c(1, 10))))
  # check if there is any filtered data
  if (nrow(filtered_data) == 0) {
    cat("No data matching the provided criteria.\n")
    return(NULL)
  }
  # determine the maximum and minimum frequency across the entire dataset
  max_freq <- max(filtered_data$Frequency, na.rm = TRUE)
  min_freq <- min(filtered_data$Frequency, na.rm = TRUE)
  # create a graph object
  g <- graph.data.frame(filtered_data[, c("Transition1", "Transition2")], directed = TRUE)
  # set edge attributes
  edge_attrs <- filtered_data$Frequency
  names(edge_attrs) <- paste(filtered_data$Transition1, filtered_data$Transition2, sep = "_")
  g <- set_edge_attr(g, "frequency", value = edge_attrs)
  # rescale edge weights based on maximum and minimum frequency
  scaled_weights <- rescale(edge_attrs, to = c(0.3, 5), from = c(min_freq, max_freq))
  # calculate Fruchterman-Reingold layout
  layout <- layout_with_fr(g)
  # get model name as a string
  model_str <- paste("Model:", model)
  # Concatenate model name in the title
  title <- model_str
  # define colors for Transition1 and Transition2 nodes
  node_colors <- ifelse(V(g)$name %in% filtered_data$Transition1, "darkorange", "gold")
  # plot the graph with title and labeled nodes
  plot(g, layout = layout, 
       main = title,
       edge.curved = FALSE, edge.arrow.size = 0.4, 
       edge.color = "black", edge.width = scaled_weights, 
       vertex.size = 13, vertex.label.dist = 0, 
       vertex.color = node_colors,
       vertex.label = V(g)$name, 
       vertex.label.color = "black",
       vertex.label.font = 2)
}

# usage:
generate_network_map(model = "Tiger", 
                     csv_file = "xxx.csv")
generate_network_map(model = "Pattern", 
                     csv_file = "xxx.csv")
generate_network_map(model = "Spots", 
                     csv_file = "xxx.csv")
generate_network_map(model = "White", 
                     csv_file = "xxx.csv")









