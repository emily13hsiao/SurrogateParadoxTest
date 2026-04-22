# -------------------------------
# Reading in file
# -------------------------------

use_pab = FALSE

n.study = 10
n.each = 100

if(n.study == 25) {tp = ""}
if(n.study == 10) {tp = "100"}



true.p = read.table(paste0("truep_", n.study, n.each, ".txt"), header=T)

check_in = function(vec, true) {
 	return(1*(vec[1]<=true & vec[2] >= true))
}

# -------------------------------
# Read in linear
# -------------------------------

big.tab = c()
for(setting in c(1:6)) {
use_spline = FALSE
degree = 1

results_linear = c()
for(i in 1:10) {
  results.temp = read.csv(paste0("res_set", setting, "_sp", use_spline, "_d", degree, "_pab",use_pab,tp,"_",i,".csv"))
  results_linear = rbind(results_linear, results.temp)
}

truth = true.p[setting,1]

apply(results_linear, 2, mean)
row_linear = c(apply(results_linear, 2, mean)[2], apply(results_linear, 2, sd)[2], apply(results_linear, 2, mean)[3], mean(apply(results_linear[,c(4:5)] , 1,check_in, true =  truth)))


# -------------------------------
# Read in cubic
# -------------------------------


use_spline = FALSE
degree = 3

results_cubic = c()
for(i in 1:10) {
  results.temp = read.csv(paste0("res_set", setting, "_sp", use_spline, "_d", degree, "_pab",use_pab,tp,"_",i,".csv"))
  results_cubic = rbind(results_cubic, results.temp)
}

apply(results_cubic, 2, mean)
row_cubic= c(apply(results_cubic, 2, mean)[2], apply(results_cubic, 2, sd)[2], apply(results_cubic, 2, mean)[3], mean(apply(results_cubic[,c(4:5)] , 1,check_in, true =  truth)))


# -------------------------------
# Read in cubic spline
# -------------------------------


use_spline = TRUE
degree = 3

results_sp = c()
for(i in 1:10) {
  results.temp = read.csv(paste0("res_set", setting, "_sp", use_spline, "_d", degree, "_pab",use_pab,tp,"_",i,".csv"))
  results_sp = rbind(results_sp, results.temp)
}

apply(results_sp, 2, mean)
row_sp= c(apply(results_sp, 2, mean)[2], apply(results_sp, 2, sd)[2], apply(results_sp, 2, mean)[3], mean(apply(results_sp[,c(4:5)] , 1,check_in, true =  truth)))

big.tab = rbind(big.tab,round(c(truth,row_linear, row_cubic, row_sp),3))
}

elliott = read.table(paste0("mv_results_", n.study,  ".txt"), header=F)
elliott=as.vector(round(elliott,3))$V1
big.tab2 = cbind(big.tab, elliott)
big.tab3 = cbind(c("1","2","3","4","5","6"), as.matrix(big.tab2))


library(quantreg)
latex.table(big.tab3, paste0("results_", n.study), dcolumn=TRUE, caption = "")
