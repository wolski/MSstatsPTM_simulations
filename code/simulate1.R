
## Load packages
#library(MSstatsPTM)
library(data.table)
library(tidyverse)
library(limma)

## Source simulation code
source("PTMsimulateExperiment.R")

## Define Modeling Functions ---------------------------------------------------
convert_to_msstats_format <- function(df){
  ## Format into MSstatsPTM input
  df$PTM$protein <- paste(df$PTM$protein, df$PTM$site, sep = "_")
  df$PTM$PrecursorCharge <- NA
  df$PTM$FragmentIon <- NA
  df$PTM$ProductCharge <- NA
  df$PTM$IsotopeLabelType <- "L"
  df$PTM$BioReplicate <- paste0(df$PTM$run, "rep")
  df$PTM$log2inty <- 2**df$PTM$log2inty
  df$PTM$PeptideSequence <- paste(df$PTM$protein, df$PTM$feature, sep = "_")

  setnames(df$PTM, c("protein", "group", "run", "log2inty"),
           c("ProteinName", "Condition", "Run", "Intensity"))

  df$PROTEIN$PrecursorCharge <- NA
  df$PROTEIN$FragmentIon <- NA
  df$PROTEIN$ProductCharge <- NA
  df$PROTEIN$IsotopeLabelType <- "L"
  df$PROTEIN$log2inty <- 2**df$PROTEIN$log2inty
  df$PROTEIN$BioReplicate <- paste0(df$PROTEIN$run, "rep")
  df$PROTEIN$PeptideSequence <- paste(df$PROTEIN$protein, df$PROTEIN$feature, sep = "_")
  setnames(df$PROTEIN, c("protein", "group", "run", "log2inty"),
           c("ProteinName", "Condition", "Run", "Intensity"))

  return(df)
}



## Simulation 1 Start ----------------------------------------------------------
## Specify types of data to simulate
s <- c(.2,.3)
reps <- c(2,3,5,10)
cond <- c(2,3,4)

param_combos <- expand.grid(s, reps, cond)
#param_combos <- param_combos[1,]
all_data <- list()
i <- 1

for (row in seq_len(nrow(param_combos))){
  ## Change in conditions
  if (param_combos[row, 3] == 2){
    del_arr <- c(0., 1.)
    del_arr_no_change <- c(0, 0)
  } else if (param_combos[row, 3] == 3) {
    del_arr <- c(0., 1., 2.)
    del_arr_no_change <- c(0, 0, 0)
  } else if (param_combos[row, 3] == 4) {
    del_arr <- c(0., 1., 2., 3.)
    del_arr_no_change <- c(0, 0, 0, 0)
  }

  if (param_combos[row, 1] == .2){
    prot_var <- .2
  } else {
    prot_var <- .3
  }

  ## Sim
  undebug(PTMsimulateExperiment)
  sim <- PTMsimulateExperiment(
    nGroup=param_combos[row, 3], nRep=param_combos[row, 2], nProtein=250, nSite=1, nFeature=2, nFeature_prot = 10,
    logAbundance=list(
      PTM=list(mu=25, delta = del_arr, sRep=param_combos[row, 1], sPeak=.25),
      PROTEIN=list(mu=25, delta = del_arr_no_change, sRep=param_combos[row, 1], sPeak=0.25))
  )
  sim_no_change1 <- PTMsimulateExperiment(
    nGroup=param_combos[row, 3], nRep=param_combos[row, 2], nProtein=250, nSite=1, nFeature=2, nFeature_prot = 10,
    logAbundance=list(
      PTM=list(mu=25, delta = del_arr, sRep=param_combos[row, 1], sPeak=0.25),
      PROTEIN=list(mu=25, delta = del_arr, sRep=param_combos[row, 1], sPeak=0.25))
  )
  sim_no_change2 <- PTMsimulateExperiment(
    nGroup=param_combos[row, 3], nRep=param_combos[row, 2], nProtein=500, nSite=1, nFeature=2, nFeature_prot = 10,
    logAbundance=list(
      PTM=list(mu=25, delta = del_arr_no_change, sRep=param_combos[row, 1], sPeak=0.25),
      PROTEIN=list(mu=25, delta = del_arr_no_change, sRep=param_combos[row, 1], sPeak=0.25))
  )
  sim_no_change1$PTM$protein <- paste0(sim_no_change1$PTM$protein, "|NoChange1")
  sim_no_change1$PROTEIN$protein <- paste0(sim_no_change1$PROTEIN$protein, "|NoChange1")
  sim_no_change2$PTM$protein <- paste0(sim_no_change2$PTM$protein, "|NoChange2")
  sim_no_change2$PROTEIN$protein <- paste0(sim_no_change2$PROTEIN$protein, "|NoChange2")

  sim_no_change_ptm <- rbindlist(list(sim_no_change1$PTM, sim_no_change2$PTM))
  sim_no_change_prot <- rbindlist(list(sim_no_change1$PROTEIN, sim_no_change2$PROTEIN))
  sim_no_change <- list(PTM = sim_no_change_ptm, PROTEIN = sim_no_change_prot)
  simOO <- convert_to_msstats_format(sim)
  sim <- convert_to_msstats_format(sim)
  sim_no_change <- convert_to_msstats_format(sim_no_change)

  sim_PTM <- rbindlist(list(sim$PTM, sim_no_change$PTM))
  sim_PROTEIN <- rbindlist(list(sim$PROTEIN, sim_no_change$PROTEIN))

  # sim_PTM[sim_PTM$Run == "R_1"]$Intensity <- NA
  # sim_PROTEIN[sim_PROTEIN$Run == "R_1"]$Intensity <- NA

  sim_combined <- list(PTM = sim_PTM, PROTEIN = sim_PROTEIN)

  ## Save to list
  all_data[[i]] <- sim_combined
  i <- i + 1
}

simulation1_data = all_data
save(simulation1_data, file = "../data/simulation1_data.rda")

