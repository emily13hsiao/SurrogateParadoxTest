Instructions for reproducing simulation results

The files run-a1-sims.R, run-a2-sims.R, and run-a3-sims.R are used to reproduce the results in Tables 1-3. 

For table 1, simply run the code in the file run-a1-sims.R

For tables 2 and 3, parallelization is required. Edit the code in run-a2-sims.R and run-a3-sims.R to the desired sample size and data generation settings, then parallelize by running:
`nice R CMD BATCH run-a2-sims-1.R
nice r CMD BATCH run-a2-sims-2.R
...
nice R CMD BATCH run-a2-sims-100.R`

Read in the results of each simulation by running the code in the files in read-a1-sims.R, read-a2-sims.R, and read-a3-sims.R.

Note that the simulations for A2 are quite computationally intensive; each setting will likely take several days to complete.






