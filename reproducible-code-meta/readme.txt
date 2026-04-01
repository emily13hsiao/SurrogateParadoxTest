Instructions for reproducing simulation and data application results in Hsiao, E. & Parast, L. (2026). Functional-Class
Meta-Analytic Framework for Quantifying Surrogate Resilience (Under Review).



Files you will need in the same directory: 
true-probabilities.R
functions-with-se.R
data_generation.R


First, run the file true-probabilities.R - this calculated the "true" resilience probability using a large sample truth calculation. It writes the true probabilities to a text file that is read in below when constructing the tables. 


# -------------------------------
# K=25, n=10
# -------------------------------


To reproduce the table, you will need to run one setting at a time. To run Setting 1, first open the main_file.R and change setting to 1. If using PAB set use_pab to TRUE, otherwise set it to FALSE. Set:
n.study = 25
n.each = 10
Next, run in parallel the following 30 files sims_em_linX.R where X is 1 through 10, sims_em_cubX.R where X is 1 through 10, and sims_em_spX.R where X is 1 through 10. Each file will write results to a csv file. 
Do this for all 6 settings.


# -------------------------------
# K=10, n=100
# -------------------------------

Repeat the process above but set:
n.study = 10
n.each = 100
in both the main_file. R and the readin_file.R



# -------------------------------
# MV results (Elliott comparison) 
# -------------------------------

Run the MV_file.R to produce the MV (last column) comparison estimates. No edits to the file should be needed.


# -------------------------------
# Final tables 
# -------------------------------

Open the readin_file.R. Set use_pab, n.study, and n.each to match your specifications for the table you would like to produce. Make sure that the following are in your directory 1) the true probabilities file 2) the 30*6 files you created from the main simulation files 3) the output from the MV_file.R script. Then run this script and this will produce the results table and write the table to a .tex file.  




