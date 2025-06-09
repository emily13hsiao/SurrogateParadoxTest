Instructions for reproducing simulation and data application results in Hsiao, E., Tian, L., & Parast, L. (2025). Resilience Measures for the Surrogate Paradox (Under Review).

There are two files in the aids-data folder. The data-cleaning.R file reads in the HIV clinical trial data (we are not authorized to share the raw data), and aids-application.R runs all functions and reproduces all results and figures in the paper. You will need to run the functions.R file in this folder before you run the aids-application.R file. 

To reproduce our simulation results, parallelization is recommended. Edit the code in run-simulations.R to reflect your desired setting, algorithm, sample size, and variance request:

setting <- 1 # Choices are: 1-10
method <- "gp_fixed" # gp_fixed, polynomial, fourier
n.A <- 400
n.B <- 200

As is, the above settings run simulation setting 1 with the Gaussian Process (algorithm), with a sample size of 400 in each group in Study A and 200 in each group in Study B. Then parallelize by running:
`nice R CMD BATCH run-simulations-1.R
nice r CMD BATCH run-simulations-2.R
...
nice R CMD BATCH run-simulations-20.R`

This runs 20 batches of 50 replications for a total of 1000 replications. When finished, read in the results of each simulation by running the code in the file read-in.R. This process will reproduce all simulation results in the paper. 

