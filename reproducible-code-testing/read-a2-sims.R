# Read in simulation results for A2

setting <- "increasing_linear"
sim_sd <- 0.10
n.iter <- 10
n=500

if(n==500) {ss =""}
if(n==250) {ss="250"}
if(n==100) {ss="100"}

# For just reading one batch
batch.num <- 78
filename <- paste0("./a2-results/", toString(setting), 
                   "sd", toString(sim_sd), 
                   "batch", toString(batch.num), ".RDS")
results <- readRDS(filename)
p_vals <- sapply(1:n.iter, function(ii) results[[ii]]$p_val)

# For reading all results
p_vals <- c()
for (batch.num in 1:100) {
  filename <- paste0("./a2-results/", toString(setting), 
                     "sd", toString(sim_sd), 
                     "batch", ss,toString(batch.num), ".RDS")
  results <- readRDS(filename)
  p_vals <- c(p_vals, sapply(1:n.iter, function(ii) results[[ii]]$p_val))
}
mean(p_vals > 0.05)


