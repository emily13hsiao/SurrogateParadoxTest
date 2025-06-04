Instructions for reproducing simulation results in Hsiao, E., Tian, L., & Parast, L. (2025). Avoiding the surrogate paradox: an empirical framework for assessing assumptions. Journal of Nonparametric Statistics, 1-22.

The files run-a1-sims.R, run-a2-sims.R, and run-a3-sims.R are used to reproduce the results in Tables 1-3. 

For table 1, simply run the code in the file run-a1-sims.R. Then, read in the results using the file read-a1-sims.R. This produces Table 1.

For tables 2 and 3, parallelization is required. Edit the code in run-a2-sims.R and run-a3-sims.R to the desired sample size and data generation settings, then parallelize by running (for assumption 2 for example) :
`nice R CMD BATCH run-a2-sims-1.R
nice r CMD BATCH run-a2-sims-2.R
...
nice R CMD BATCH run-a2-sims-100.R`

Read in the results of each simulation by running the code in the files read-a2-sims.R (or read-a3-sims.R for assumption 3). This produces Table 2 and 3.

Note that the simulations for A2 are quite computationally intensive; each setting will likely take several days to complete.






