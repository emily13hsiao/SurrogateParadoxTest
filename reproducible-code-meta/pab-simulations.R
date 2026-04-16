# Simulations with SE
source("functions-with-se.R")
source("data_generation.R")

RNGkind("L'Ecuyer-CMRG")
set.seed(1)

# DATA GENERATION SETTINGS
setting = 5 # Options are 1, 2, 3, 4, 5, 6
use_spline = TRUE # Options are TRUE or FALSE
degree = 3 # For our sims, 1 or 3
use_pab = TRUE # Either use PAB or fully analytic method

n.studies = 10
n.each = 100

sim.reps = 1000

results_list <- vector(mode = "list", length = sim.reps)
se_list <- rep(NA, length = sim.reps)
p_list <- rep(NA, length = sim.reps)
for (iter in 1:sim.reps) {
  
  print(iter)
  data <- generate_data(setting, n.studies, n.each)
  
  r <- full_procedure(data$data, data$s0.B, data$s1.B, use_spline = use_spline, 
                      degree = degree, calculate_se = TRUE, try_analytic = use_pab)
  print(r$p)
  se_list[iter] <- r$se_p
  p_list[iter] <- r$p
  results_list[[iter]] <- r
}
mean(p_list)
sd(p_list)
mean(se_list)

saveRDS(results_list, file = paste0("setting", setting, "spline", use_spline, "degree", degree, "results.RDS"))








