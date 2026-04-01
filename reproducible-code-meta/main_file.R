# -------------------------------
# Main file for parallel
# -------------------------------

# -------------------------------
# Settings
# -------------------------------

setting = 1 # Options are 1, 2, 3, 4, 5, 6
use_pab = FALSE # Either use PAB or fully analytic method

n.study = 10
n.each = 100

sim.reps = 100

# -------------------------------
# Should not need to edit anything below
# -------------------------------




source("functions-with-se.R")
source("data_generation.R")

set.seed(100*round_file)


results_list <- c()

for (iter in 1:sim.reps) {
  
  print(iter)
  data <- generate_data(setting, n.study, n.each)
  print(system.time({
  r <- full_procedure(data$data, data$s0.B, data$s1.B, use_spline = use_spline, 
                      degree = degree, calculate_se = TRUE, try_analytic = use_pab, n_bootstrap = 200)
    print(unlist(r[-length(r)])) 
    results_list <- rbind(results_list,unlist(r[-length(r)]))
  }))
}

if(n.study == 25) {tp = ""}
if(n.study == 10) {tp = "100"}
write.csv(results_list, file = paste0("res_set", setting, "_sp", use_spline, "_d", degree, "_pab",use_pab,tp,"_",round_file,".csv"))








