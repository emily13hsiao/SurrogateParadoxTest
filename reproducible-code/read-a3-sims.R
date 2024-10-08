# Read in a3 sims

n <- 500
n.iter <- 10

###############################################################################
###########################      Read Results      ############################
###############################################################################

# Setting 1
p_vals1 <- c()
for (batch.num in (1:100)) {
  results <- readRDS(paste0("./a3-results/setting1n", toString(n), "batch", 
                            toString(batch.num), ".RDS"))
  p_vals1 <- c(p_vals1, sapply(1:n.iter, function(ii) results[[ii]]$p_value))
}
mean(p_vals1 > 0.05) # Ideally close to 1

# Setting 2
p_vals2 <- c()
for (batch.num in (1:100)) {
  results <- readRDS(paste0("./a3-results/setting2n", toString(n), "batch", 
                            toString(batch.num), ".RDS"))
  p_vals2 <- c(p_vals2, sapply(1:n.iter, function(ii) results[[ii]]$p_value))
}
mean(p_vals2 > 0.05) # Ideally close to 1


# Setting 3
p_vals3 <- c()
for (batch.num in (1:100)) {
  results <- readRDS(paste0("./a3-results/setting3n", toString(n), "batch", 
                            toString(batch.num), ".RDS"))
  p_vals3 <- c(p_vals3, sapply(1:n.iter, function(ii) results[[ii]]$p_value))
}
mean(p_vals3 > 0.05) # Ideally close to 1

# Setting 4
p_vals4 <- c()
for (batch.num in (1:100)) {
  results <- readRDS(paste0("./a3-results/setting4n", toString(n), "batch", 
                            toString(batch.num), ".RDS"))
  p_vals4 <- c(p_vals4, sapply(1:n.iter, function(ii) results[[ii]]$p_value))
}
mean(p_vals4 > 0.05) # Ideally close to 1

# Setting 5
p_vals5 <- c()
for (batch.num in (1:100)) {
  results <- readRDS(paste0("./a3-results/setting5n", toString(n), "batch", 
                            toString(batch.num), ".RDS"))
  p_vals5 <- c(p_vals5, sapply(1:n.iter, function(ii) results[[ii]]$p_value))
}
mean(p_vals5 > 0.05) # Ideally close to 1

# Setting 6
p_vals6 <- c()
for (batch.num in (1:100)) {
  results <- readRDS(paste0("./a3-results/setting6n", toString(n), "batch", 
                            toString(batch.num), ".RDS"))
  p_vals6 <- c(p_vals6, sapply(1:n.iter, function(ii) results[[ii]]$p_value))
}
mean(p_vals6 > 0.05) # Ideally close to 0

# Setting 7
p_vals7 <- c()
for (batch.num in (1:100)) {
  results <- readRDS(paste0("./a3-results/setting7n", toString(n), "batch", 
                            toString(batch.num), ".RDS"))
  p_vals7 <- c(p_vals7, sapply(1:n.iter, function(ii) results[[ii]]$p_value))
}
mean(p_vals7 > 0.05) # Ideally close to 0

# Setting 8
p_vals8 <- c()
for (batch.num in (1:100)) {
  results <- readRDS(paste0("./a3-results/setting8n", toString(n), "batch", 
                            toString(batch.num), ".RDS"))
  p_vals8 <- c(p_vals8, sapply(1:n.iter, function(ii) results[[ii]]$p_value))
}
mean(p_vals8 > 0.05) # Ideally close to 0

# Setting 9
p_vals9 <- c()
for (batch.num in (1:100)) {
  results <- readRDS(paste0("./a3-results/setting9n", toString(n), "batch", 
                            toString(batch.num), ".RDS"))
  p_vals9 <- c(p_vals9, sapply(1:n.iter, function(ii) results[[ii]]$p_value))
}
mean(p_vals9 > 0.05) # Ideally close to 0


###############################################################################
########################      Large Sample Size      ##########################
###############################################################################

# Setting 7, n = 1000
n <- 1000
p_vals1000 <- c()
for (batch.num in (1:100)) {
  results <- readRDS(paste0("./a3-results/setting7n", toString(n), "batch", 
                            toString(batch.num), ".RDS"))
  p_vals1000 <- c(p_vals1000, sapply(1:n.iter, function(ii) results[[ii]]$p_value))
}
mean(p_vals1000 > 0.05) # Ideally close to 0

# Setting 7, n = 2000
n <- 2000
p_vals2000 <- c()
for (batch.num in (1:100)) {
  results <- readRDS(paste0("./a3-results/setting7n", toString(n), "batch", 
                            toString(batch.num), ".RDS"))
  p_vals2000 <- c(p_vals2000, sapply(1:n.iter, function(ii) results[[ii]]$p_value))
}
mean(p_vals2000 > 0.05) # Ideally close to 0




