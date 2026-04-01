#################################################
## Main data application
#################################################

# Load required packages

library(Surrogate)
library(ggplot2)
library(tidyverse)
library(dplyr)


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

# Plot PANSS vs BPRS by Treat, color by InvestId
ggplot(Schizo, aes(x = BPRS, y = PANSS, color = InvestId)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  facet_wrap(~Treat, ncol = 2, labeller = label_both) +
  labs(
    y = "PANSS score",
    x = "BPRS score",
    title = "Relationship between PANSS and BPRS by Treat group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

#################################################

###run functions file
source("functions-with-se.R")
#for both PANSS and BPRS higher is worse, so we hope for a reduction. Looking at PANSS binary, this must be coded as measurement at end minus measurement at the beginning because negative numbers are GOOD and are coded as reductions. 
#coding such that higher numbers indicate better



Schizo_flip = Schizo
Schizo_flip$PANSS = -1*Schizo$PANSS
Schizo_flip$BPRS = -1*Schizo$BPRS
Schizo_flip$Treat = 1*(Schizo$Treat == 1)
Schizo_flip = na.omit(Schizo_flip)

# look at how many patients per practice

table(Schizo_flip$InvestId)

#summarize treatment effect on s and on y, as well as sample size in each arm, for each practice

data_sub = Schizo_flip
info_mat = c()
id_list = sort(unique(data_sub$InvestId))
for(j in 1:length(id_list)){
  data_sub_id = data_sub[data_sub$InvestId == id_list[j],]
  treat.s = mean(data_sub_id$BPRS[data_sub_id$Treat == 1] ) - mean(data_sub_id$BPRS[data_sub_id$Treat == 0] ) 
  treat.y = mean(data_sub_id$PANSS[data_sub_id$Treat == 1] ) - mean(data_sub_id$PANSS[data_sub_id$Treat == 0] ) 
  
  info_mat = rbind(info_mat, c(id_list[j], treat.s, treat.y, length(data_sub_id$BPRS[data_sub_id$Treat == 0]), length(data_sub_id$BPRS[data_sub_id$Treat == 1])))
	}

info_mat

#set minumum number per practice per treatment group

min_number = 6

#identify subset of studies to include

studies_want = info_mat[info_mat[,4] >=min_number & info_mat[,5] >=min_number, ]

#decide on two different studies we want to pull out as study B

studies_b = c(50,3)
studies_a = setdiff(studies_want[,1],studies_b)

#################################################
## Setup with study B choice #1
#################################################

example_data = list(data =c(), s0.B = Schizo_flip$BPRS[Schizo_flip$InvestId == studies_b[1] & Schizo_flip$Treat == 0], s1.B = Schizo_flip$BPRS[Schizo_flip$InvestId == studies_b[1] & Schizo_flip$Treat == 1])

dataset = data.frame("S" = Schizo_flip$BPRS, "Y" = Schizo_flip$PANSS, "G" = Schizo_flip$Treat, "study" = Schizo_flip$InvestId)
dataset = dataset[dataset$study %in% studies_a,]

example_data$data = dataset

#run procedure
set.seed(1)
prob_degree3spline <- full_procedure(example_data$data, example_data$s0.B, example_data$s1.B, degree=3, use_spline = TRUE, calculate_se = TRUE, try_analytic = FALSE) 

#results
prob_degree3spline

##################################################
## Setup with study B choice #2
#################################################

example_data = list(data =c(), s0.B = Schizo_flip$BPRS[Schizo_flip$InvestId == studies_b[2] & Schizo_flip$Treat == 0], s1.B = Schizo_flip$BPRS[Schizo_flip$InvestId == studies_b[2] & Schizo_flip$Treat == 1])

dataset = data.frame("S" = Schizo_flip$BPRS, "Y" = Schizo_flip$PANSS, "G" = Schizo_flip$Treat, "study" = Schizo_flip$InvestId)
dataset = dataset[dataset$study %in% studies_a,]

example_data$data = dataset

#run procedure
set.seed(100)
prob_degree3spline <- full_procedure(example_data$data, example_data$s0.B, example_data$s1.B, degree=3, use_spline = TRUE, calculate_se = TRUE, try_analytic = FALSE) 

#results
prob_degree3spline


##################################################
## Make plot for paper
#################################################


# Prepare data
Schizo_flip <- Schizo_flip %>%
  mutate(
    Treat = as.factor(Treat),
    InvestId = as.factor(InvestId)
  )

Schizo_flip_studies_a = Schizo_flip[Schizo_flip$InvestId %in% studies_a,]

# Plot PANSS vs BPRS by Treat, color by InvestId
nw_smooth <- function(x, y, bw, trim = 0.02) {
  lo <- quantile(x, trim)
  hi <- quantile(x, 1 - trim)
  grid <- seq(lo, hi, length.out = 200)

  K <- function(u) exp(-0.5 * u^2) / sqrt(2 * pi)

  yhat <- sapply(grid, function(g) {
    w <- K((g - x) / bw)
    sum(w * y) / sum(w)
  })

  data.frame(x = grid, y = yhat)
}

library(purrr)

bw <- 8   

smoothed_df <- Schizo_flip_studies_a %>%
  group_by(Treat, InvestId) %>%
  group_modify(~ {
    df <- .
    nw_smooth(df$BPRS, df$PANSS, bw = bw)
  })


ggplot() +
  geom_point(data = Schizo_flip_studies_a,
             aes(x = BPRS, y = PANSS, color = InvestId),
             alpha = 0.7) +
  geom_line(data = smoothed_df,
            aes(x = x, y = y, color = InvestId),
            linewidth = 1) +
  facet_wrap(
    ~Treat,
    ncol = 2,
    labeller = labeller(Treat = c(`0`="Control", `1`="Treatment"))
  ) +
  labs(
    y = "PANSS score",
    x = "BPRS score"
   # title = "Relationship between PANSS and BPRS by Treatment group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )
