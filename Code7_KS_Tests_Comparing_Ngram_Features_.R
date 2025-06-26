
### datasets required: nonnatural = dataset containing columns:
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
#'                    

# KS test (fig 6a - duration)
ks_result_duration <- ks.test(nonnatural$RecordingDuration[nonnatural$Model == "Non-natural"],
                              nonnatural$RecordingDuration[nonnatural$Model == "Tiger"])
print(ks_result_duration)

# KS test (fig 6b - call rate)
ks_result_callrate <- ks.test(nonnatural$RecCallRate[nonnatural$Model == "Non-natural"],
                              nonnatural$RecCallRate[nonnatural$Model == "Tiger"])
print(ks_result_callrate)

# KS test - (6c - ngrams length)
ks_result_ngram <- ks.test(nonnatural$ngram[nonnatural$Model == "Non-natural"],
                           nonnatural$ngram[nonnatural$Model == "Tiger"])
print(ks_result_ngram)





### with bootstrap (confidence intervals):

library(boot)

### kolmogorov-Smirnov test – duration (fig. 6a: longer)
ks_result_duration <- ks.test(nonnatural$RecordingDuration[nonnatural$Model == "Non-natural"],
                              nonnatural$RecordingDuration[nonnatural$Model == "Tiger"])
print(ks_result_duration)

# bootstrap CI for D
duration_data <- nonnatural[nonnatural$Model %in% c("Non-natural", "Tiger"), c("Model", "RecordingDuration")]
set.seed(123)
boot_duration <- boot(data = duration_data, 
                      statistic = function(data, i) {
                        d <- data[i, ]
                        x <- d$RecordingDuration[d$Model == "Non-natural"]
                        y <- d$RecordingDuration[d$Model == "Tiger"]
                        as.numeric(ks.test(x, y)$statistic)
                      }, 
                      R = 1000)
boot.ci(boot_duration, type = "perc")


### KS TEST – call rate (fig. 6b: faster)
ks_result_callrate <- ks.test(nonnatural$RecCallRate[nonnatural$Model == "Non-natural"],
                              nonnatural$RecCallRate[nonnatural$Model == "Tiger"])
print(ks_result_callrate)

# bootstrap CI for D
callrate_data <- nonnatural[nonnatural$Model %in% c("Non-natural", "Tiger"), c("Model", "RecCallRate")]
set.seed(123)
boot_callrate <- boot(data = callrate_data, 
                      statistic = function(data, i) {
                        d <- data[i, ]
                        x <- d$RecCallRate[d$Model == "Non-natural"]
                        y <- d$RecCallRate[d$Model == "Tiger"]
                        as.numeric(ks.test(x, y)$statistic)
                      }, 
                      R = 1000)
boot.ci(boot_callrate, type = "perc")


### KS TEST – N-gram size (fig. 6c: larger n-grams)
ks_result_ngram <- ks.test(nonnatural$ngram[nonnatural$Model == "Non-natural"],
                           nonnatural$ngram[nonnatural$Model == "Tiger"])
print(ks_result_ngram)

# bootstrap CI for D
ngram_data <- nonnatural[nonnatural$Model %in% c("Non-natural", "Tiger"), c("Model", "ngram")]
set.seed(123)
boot_ngram <- boot(data = ngram_data, 
                   statistic = function(data, i) {
                     d <- data[i, ]
                     x <- d$ngram[d$Model == "Non-natural"]
                     y <- d$ngram[d$Model == "Tiger"]
                     as.numeric(ks.test(x, y)$statistic)
                   }, 
                   R = 1000)
boot.ci(boot_ngram, type = "perc")

