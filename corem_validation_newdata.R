# how do corems perform on DREAM5 data
# which corems so signs of activity in new data set
# are the number of active conditions greater than random
# how does distribution of active conditions for corems in DREAM5 data set compare to Distiller data set

# starting dir: /home/abrooks/Documents/EGRIN2/Eco_ensemble_2/validation/DREAM5

# source scripts
source("~/Documents/git/Corems/processEGRIN.R")

# load data
load("DREAM5_normed_ratios.RData")

# make a random condition filehash for DREAM5 data 

