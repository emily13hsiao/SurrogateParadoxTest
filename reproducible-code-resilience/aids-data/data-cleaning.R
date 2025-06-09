
############################################
#### USING ACTG 320 AS STUDY A
#### Reference is Hammer et al. (1997) NEJM paper
#### Y = change in RNA from baseline to 24 weeks; coded as baseline - 24 weeks because lower RNA is better
#### S = change in CD4 from baseline to 24 weeks; coded as 24 weeks - baseline because higher CD4 is better
################################################

tmpdir <- "./aids-data/"
aids.base <- read.table(paste(tmpdir, "actg320.base.dat", sep=""))
aids.rna <- read.table(paste(tmpdir, "actg320.rna", sep=""))
aids.cd4 <- read.table(paste(tmpdir, "actg320.cd4", sep=""))
aids.trt <- read.table(paste(tmpdir, "actg320.trt", sep=""))

n.total = dim(aids.base)[1]
aids.base$TREAT = vector(length = n.total)
aids.base$CD4BASE = vector(length = n.total)
aids.base$CD424 = vector(length = n.total)
aids.base$RNABASE = vector(length = n.total)
aids.base$RNA24 = vector(length = n.total)

for(i in 1:n.total) {
	aids.base$TREAT[i] = aids.trt[aids.trt$PATID == aids.base$PATID[i],2]
	if(dim(aids.cd4[aids.cd4$PATID == aids.base$PATID[i] & aids.cd4$WEEK == 0,])[1] > 0) {aids.base$CD4BASE[i] = aids.cd4[aids.cd4$PATID == aids.base$PATID[i] & aids.cd4$WEEK == 0,]$CD4} else {aids.base$CD4BASE[i] = NA}
	if(dim(aids.cd4[aids.cd4$PATID == aids.base$PATID[i] & aids.cd4$WEEK == 24,])[1] > 0) {aids.base$CD424[i] = aids.cd4[aids.cd4$PATID == aids.base$PATID[i] & aids.cd4$WEEK == 24,]$CD4} else {aids.base$CD424[i] = NA}
	if(dim(aids.rna[aids.rna$PATID == aids.base$PATID[i] & aids.rna$WEEK == 0,])[1] >0) {aids.base$RNABASE[i] = aids.rna[aids.rna$PATID == aids.base$PATID[i] & aids.rna$WEEK == 0,]$ULOGRNA} else {aids.base$RNABASE[i] = NA}
	if(dim(aids.rna[aids.rna$PATID == aids.base$PATID[i] & aids.rna$WEEK == 24,])[1] >0) {aids.base$RNA24[i] = aids.rna[aids.rna$PATID == aids.base$PATID[i] & aids.rna$WEEK == 24,]$ULOGRNA} else {aids.base$RNA24[i] = NA}
}

aids.base$CD4CHANGE = aids.base$CD424 - aids.base$CD4BASE
aids.base$RNACHANGE = aids.base$RNABASE - aids.base$RNA24

aids.base = aids.base[!(is.na(aids.base$RNACHANGE) | is.na(aids.base$CD4CHANGE)),]

y1 = aids.base$RNACHANGE[aids.base$TREAT == 2]
y0 = aids.base$RNACHANGE[aids.base$TREAT == 1]
s1 = aids.base$CD4CHANGE[aids.base$TREAT == 2]
s0 = aids.base$CD4CHANGE[aids.base$TREAT == 1]


############################################
#### USING ACTG 193A AS STUDY B
################################################

#use the data that ACTG gave me with DOES have RNA, but only for a small random subset. Not ideal but nothing else is working.
#setwd("/Users/parastlm/Dropbox/RAND/Proposals/R21_Surrogate_Marker/AIDS data/ACTG 193A Data")
#Zidovudine and didanosine (same time) . (2 NRTIs) = CONTROL, GROUP 0, ZDV+ddI
#Zidovudine and didanosine and nevirapine . (2 NRTIs plus NNRTI), TREAT, GROUP 1 ZDV+ddI+NVP
#OK SO GROUP 0 = ZDV+ddI; GROUP 1 = ZDV+ddI+NVP

aids.193 = read.csv("./aids-data/193A_Base.csv", header = T)
aids.193 = aids.193[aids.193$TRT == "ZDV+ddI" | aids.193$TRT == "ZDV+ddI+NVP",]
aids.193$TRT = 1*(aids.193$TRT == "ZDV+ddI+NVP")
aids.193.cd4 = read.csv("./aids-data/193A_CD4.csv", header = T)

n.total = dim(aids.193)[1]
aids.193$CD4BASE = vector(length = n.total)
aids.193$CD424 = vector(length = n.total)


#get baseline and 24 week CD4 

for(i in 1:n.total) {
	if(dim(aids.193.cd4[aids.193.cd4$PATID == aids.193$PATID[i] & aids.193.cd4$VISIT == "Baseline",])[1] > 0) {aids.193$CD4BASE[i] = aids.193.cd4[aids.193.cd4$PATID == aids.193$PATID[i] & aids.193.cd4$VISIT == "Baseline",]$CD4} else {aids.193$CD4BASE[i] = NA}
	if(dim(aids.193.cd4[aids.193.cd4$PATID == aids.193$PATID[i] & aids.193.cd4$VISIT == "Week 24",])[1] > 0) {aids.193$CD424[i] = aids.193.cd4[aids.193.cd4$PATID == aids.193$PATID[i] & aids.193.cd4$VISIT == "Week 24",]$CD4} else {aids.193$CD424[i] = NA}
	}

#calculate change
aids.193$CD4CHANGE = aids.193$CD424 - aids.193$CD4BASE
#huge drop here
aids.193 = aids.193[!is.na(aids.193$CD4CHANGE),]

studyb.s1 = aids.193$CD4CHANGE[aids.193$TRT == 1]
studyb.s0 = aids.193$CD4CHANGE[aids.193$TRT == 0]

# Save the cleaned stuff
saveRDS(s0, "s0_A.RDS")
saveRDS(y0, "y0_A.RDS")
saveRDS(s1, "s1_A.RDS")
saveRDS(y1, "y1_A.RDS")
saveRDS(studyb.s0, "s0_B.RDS")
saveRDS(studyb.s1, "s1_B.RDS")

