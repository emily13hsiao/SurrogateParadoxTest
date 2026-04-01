# -------------------------------
# File to calculate comparison MV (Elliott)
# -------------------------------

# -------------------------------
# Should not need to edit
# -------------------------------

source("functions-with-se.R")
source("data_generation.R")


# -------------------------------
# K=25
# -------------------------------
set.seed(1)
sim.reps = 1000

n.study = 25
n.each = 10

mv_results = vector(length = 6)
for(uu in 1:6) {
  setting = uu # 
 
  probs_elliott = vector(length=sim.reps)

  for (iter in 1:sim.reps) {
  
    print(iter)
    data <- generate_data(setting, n.study, n.each)
    probs_elliott[iter] = elliott_prob(data$data, data$s0.B, data$s1.B) 
  }
  mv_results[uu] = mean(probs_elliott)
}

write.table(mv_results,paste0("mv_results_", n.study,".txt"), quote = FALSE, row.names = FALSE, col.names=FALSE)

# -------------------------------
# K=10
# -------------------------------
set.seed(1)
sim.reps = 1000

n.study = 10
n.each = 100

mv_results = vector(length = 6)
for(uu in 1:6) {
  setting = uu # 
  
  probs_elliott = vector(length=sim.reps)
  
  for (iter in 1:sim.reps) {
    
    print(iter)
    data <- generate_data(setting, n.study, n.each)
    probs_elliott[iter] = elliott_prob(data$data, data$s0.B, data$s1.B) 
  }
  mv_results[uu] = mean(probs_elliott)
}

write.table(mv_results,paste0("mv_results_", n.study,".txt"), quote = FALSE, row.names = FALSE, col.names=FALSE)

