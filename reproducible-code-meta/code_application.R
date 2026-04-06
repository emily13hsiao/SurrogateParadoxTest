# Reproduce application done by Layla

library(Surrogate)
library(tidyverse)
library(dplyr)
source("functions-with-se.R")

# Set seed
RNGkind("L'Ecuyer-CMRG")
set.seed(1)

#################################################
## Schizo data
#################################################

# Load data
data("Schizo")

# Prepare data
Schizo <- Schizo %>%
  mutate(
    Treat = as.factor(Treat),
    InvestId = as.factor(InvestId)
  )

#################################################

# For both PANSS and BPRS higher is worse, so we hope for a reduction. 
# Looking at PANSS binary, this must be coded as measurement at end minus 
# measurement at the beginning because negative numbers are GOOD and are coded as reductions. 
# making higher numbers indicate better

Schizo_flip = Schizo
Schizo_flip$PANSS = -1*Schizo$PANSS
Schizo_flip$BPRS = -1*Schizo$BPRS
Schizo_flip$Treat = 1*(Schizo$Treat == 1)
Schizo_flip = na.omit(Schizo_flip)

# Rename InvestId to PracticeID
Schizo_flip = Schizo_flip %>% mutate(PracticeID = InvestId)
# Just going to select the relevant columns
Schizo_flip = Schizo_flip %>% dplyr::select(Id, PracticeID, Treat, PANSS, BPRS)

# Look at how many patients per practice
table(Schizo_flip$PracticeID)

# Identify practices with at least 6 patients

big = which(table(Schizo_flip$PracticeID)>=5)
data_sub = Schizo_flip[Schizo_flip$PracticeID %in% big,]
data_sub

# Summarize treatment effect on S and on Y, as well as sample size in each arm, 
# for each practice

info_mat = c()
id_list = sort(unique(data_sub$PracticeID))
for(j in 1:length(id_list)){
  data_sub_id = data_sub[data_sub$PracticeID == id_list[j],]
  treat.s = mean(data_sub_id$BPRS[data_sub_id$Treat == 1] ) - mean(data_sub_id$BPRS[data_sub_id$Treat == 0] ) 
  treat.y = mean(data_sub_id$PANSS[data_sub_id$Treat == 1] ) - mean(data_sub_id$PANSS[data_sub_id$Treat == 0] ) 
  
  info_mat = rbind(info_mat, c(id_list[j], treat.s, treat.y, length(data_sub_id$BPRS[data_sub_id$Treat == 0]), length(data_sub_id$BPRS[data_sub_id$Treat == 1])))
}
colnames(info_mat) <- c("PracticeID", "DeltaS", "DeltaY", "n0", "n1")
info_mat

# Set minumum number per practice per treatment group and subset

min_number = 6
studies_want = info_mat[info_mat[,4] >=min_number & info_mat[,5] >=min_number, ]
studies_want <- data.frame(studies_want)

# We are going to use Study 50 and Study 3 as "Study B"

studies_b = c(50, 3)
studies_a = setdiff(studies_want[,1],studies_b)

# Make full dataset
dataset = data.frame("S" = Schizo_flip$BPRS, 
                     "Y" = Schizo_flip$PANSS, 
                     "G" = Schizo_flip$Treat, 
                     "study" = Schizo_flip$PracticeID)
dataset = dataset[dataset$study %in% studies_a,]

#################################################
## Setup with study B choice #1
#################################################

# SET UP DATA FORMAT
study3data <- Schizo_flip %>% filter(PracticeID == 3) %>% data.frame()
s0.B <- filter(study3data, Treat == 0)$BPRS
s1.B <- filter(study3data, Treat == 1)$BPRS

# Run procedure
linear_model <- full_procedure(dataset, s0.B, s1.B, use_spline = FALSE, 
                               degree = 1, try_analytic = FALSE)
linear_data <- c(linear_model$p, linear_model$se_p, 
                 linear_model$p - 1.96 * linear_model$se_p,
                 linear_model$p + 1.96 * linear_model$se_p)

cubic_model <- full_procedure(dataset, s0.B, s1.B, use_spline = FALSE, 
                              degree = 3, try_analytic = FALSE)
cubic_data <- c(cubic_model$p, cubic_model$se_p, 
                 cubic_model$p - 1.96 * cubic_model$se_p,
                 cubic_model$p + 1.96 * cubic_model$se_p)

spline_model <- full_procedure(dataset, s0.B, s1.B, use_spline = TRUE, 
                               degree = 3, try_analytic = FALSE)
spline_data <- c(spline_model$p, spline_model$se_p, 
                 spline_model$p - 1.96 * spline_model$se_p,
                 spline_model$p + 1.96 * spline_model$se_p)

# Now to make a table of: p, SE_p, lower endpoint, upper endpoint
study3row <- c(linear_data, cubic_data, spline_data)

#################################################
## Setup with study B choice #2
#################################################

# SET UP DATA FORMAT
study50data <- Schizo_flip %>% filter(PracticeID == 50) %>% data.frame()
s0.B <- filter(study50data, Treat == 0)$BPRS
s1.B <- filter(study50data, Treat == 1)$BPRS

# Run procedure
linear_model <- full_procedure(dataset, s0.B, s1.B, use_spline = FALSE, 
                               degree = 1, try_analytic = FALSE)
linear_data <- c(linear_model$p, linear_model$se_p, 
                 linear_model$p - 1.96 * linear_model$se_p,
                 linear_model$p + 1.96 * linear_model$se_p)

cubic_model <- full_procedure(dataset, s0.B, s1.B, use_spline = FALSE, 
                              degree = 3, try_analytic = FALSE)
cubic_data <- c(cubic_model$p, cubic_model$se_p, 
                cubic_model$p - 1.96 * cubic_model$se_p,
                cubic_model$p + 1.96 * cubic_model$se_p)

spline_model <- full_procedure(dataset, s0.B, s1.B, use_spline = TRUE, 
                               degree = 3, try_analytic = FALSE)
spline_data <- c(spline_model$p, spline_model$se_p, 
                 spline_model$p - 1.96 * spline_model$se_p,
                 spline_model$p + 1.96 * spline_model$se_p)

# Now to make a table of: p, SE_p, lower endpoint, upper endpoint
study50row <- c(linear_data, cubic_data, spline_data)


# study3row
# 0.0250000  0.3173252 -0.5969573  0.6469573  0.0050000  0.3371543 -0.6558224  0.6658224  0.0150000  0.3252365 -0.6224636  0.6524636

# study50row
# 0.1400000  0.3335185 -0.5136962  0.7936962  0.1750000  0.3239794 -0.4599997  0.8099997  0.1400000  0.3383623 -0.5231901  0.8031901











