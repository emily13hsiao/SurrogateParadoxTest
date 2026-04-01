#################################################
## Calculates true paradox probabilties
#################################################

source("functions-with-se.R")
source("data_generation.R")


settings <- c(1, 2, 3, 4, 5, 6)
n.rep <- 1000

#################################################
## All settings, n.each=100, n.study=10
#################################################

n.each = 100
n.study = 10

true.p = c()
set.seed(1)

for (setting in settings) {
  is_paradox <- rep(NA, n.rep)
  for (i in 1:n.rep) {
    data <- generate_data(setting, n.study, n.each)
    delta <- mean(data$y1.B) - mean(data$y0.B)
    is_paradox[i] <- delta < 0
  }
  print(paste("Setting ", setting, "paradox: ", mean(is_paradox)))
  true.p = c(true.p,mean(is_paradox))
}

write.table(true.p, paste0("simulation_results/truep_", n.study, n.each, ".txt"), quote = FALSE, row.names = FALSE)


#################################################
## All settings, n.each=10, n.study=25
#################################################

n.each = 10
n.study = 25

true.p = c()
set.seed(1)

for (setting in settings) {
  is_paradox <- rep(NA, n.rep)
  for (i in 1:n.rep) {
    data <- generate_data(setting, n.study, n.each)
    delta <- mean(data$y1.B) - mean(data$y0.B)
    is_paradox[i] <- delta < 0
  }
  print(paste("Setting ", setting, "paradox: ", mean(is_paradox)))
  true.p = c(true.p,mean(is_paradox))
}

write.table(true.p, paste0("simulation_results/truep_", n.study, n.each, ".txt"), quote = FALSE, row.names = FALSE)


