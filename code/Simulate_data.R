
## Load packages
library(MSstatsPTM)
library(data.table)
library(tidyverse)
library(limma)

## Source simulation code
source("PTMsimulateExperiment.R")

## Define Modeling Functions ---------------------------------------------------
## Format simulation data into MSstatsPTM input
convert_to_msstats_format <- function(df){
  
  ## modified data, edit column names and add missing
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

  ## unmodified data, edit column names and add missing
  df$PROTEIN$PrecursorCharge <- NA
  df$PROTEIN$FragmentIon <- NA
  df$PROTEIN$ProductCharge <- NA
  df$PROTEIN$IsotopeLabelType <- "L"
  df$PROTEIN$log2inty <- 2**df$PROTEIN$log2inty
  df$PROTEIN$BioReplicate <- paste0(df$PROTEIN$run, "rep")
  df$PROTEIN$PeptideSequence <- paste(df$PROTEIN$protein, df$PROTEIN$feature, 
                                      sep = "_")
  setnames(df$PROTEIN, c("protein", "group", "run", "log2inty"),
           c("ProteinName", "Condition", "Run", "Intensity"))

  return(df)
}

## Run ttest model using run level summarized data
run_ttest <- function(summarized_data){

  ptms <- unique(summarized_data$PTM)

  ttest_temp <- data.table()
  ttest_adj_temp <- data.table()
  ## Run ttest for each ptm
  for (p in seq_along(ptms)){
    temp_joined <- summarized_data %>% filter(PTM == ptms[[p]])
    groups <- (temp_joined %>% distinct(Condition))[[1]]
    t <- 2

    ## Loop over groups
    for (g in 1:(length(groups)-1)){
      for (g2 in 2:length(groups)){
        tryCatch({ttest_ptm <- t.test((temp_joined %>% 
                                         filter(Condition == groups[g]) %>% 
                                         select(Abundance.x))[[1]],
                                      (temp_joined %>% 
                                         filter(Condition == groups[g2]) %>% 
                                         select(Abundance.x))[[1]])
        ttest_temp <- rbindlist(list(ttest_temp, 
                                     data.table(ptm = ptms[[p]],
                                        label = paste(groups[g], "vs", 
                                                      groups[g2], sep = " "),
                                        pval = ttest_ptm$p.value,
                                        tstat = ttest_ptm$statistic[[1]],
                                        SE = ttest_ptm$stderr,
                                        df = ttest_ptm$parameter[[1]],
              estimate = ttest_ptm$estimate[[2]] - ttest_ptm$estimate[[1]])))},
                 error=function(e){cat("ERROR :", as.character(p), "\n")})
        tryCatch({ttest_adj_ptm <- t.test((temp_joined %>% 
                                             filter(Condition == groups[g]) %>% 
                                             select(Adj_Abundance))[[1]],
                                          (temp_joined %>% 
                                             filter(Condition == groups[g2]) %>% 
                                             select(Adj_Abundance))[[1]])
        ttest_adj_temp <- rbindlist(list(ttest_adj_temp, 
                                         data.table(ptm = ptms[[p]],
                                          label = paste(groups[g], "vs", 
                                                        groups[g2], sep = " "),
                                          pval = ttest_adj_ptm$p.value,
                                          tstat = ttest_adj_ptm$statistic[[1]],
                                          SE = ttest_adj_ptm$stderr,
                                          df = ttest_adj_ptm$parameter[[1]],
      estimate = ttest_adj_ptm$estimate[[2]] - ttest_adj_ptm$estimate[[1]])))},
                 error=function(e){cat("ERROR :", as.character(p), "\n")})


      }
    }
  }
  return(list(ttest_temp = ttest_temp, ttest_adj_temp = ttest_adj_temp))
}

## Limma pairwise function
design.pairs <- function(levels) {
    n <- length(levels)
    design <- matrix(0,n,choose(n,2))
    rownames(design) <- levels
    colnames(design) <- 1:choose(n,2)
    k <- 0
    for (i in 1:(n-1))
      for (j in (i+1):n) {
        k <- k+1
        design[i,k] <- 1
        design[j,k] <- -1
        colnames(design)[k] <- paste(levels[i],"-",levels[j],sep="")
      }
    design
}

## Fit limma model for given dataset
fit_limma <- function(summarized_data, conditions, runs){

  ## Convert into format required for limma
  input <- data.frame(summarized_data %>% 
                        select(PTM, Run, Condition, Abundance.x) %>%
                        pivot_wider(names_from = c(Condition, Run), 
                                    values_from = Abundance.x,
                                    names_sort = TRUE))
  rownames(input) <- input$PTM
  input <- input %>% select(-PTM)
  input_adj <- data.frame(summarized_data %>% 
                            select(PTM, Run, Condition, Adj_Abundance) %>%
                            pivot_wider(names_from = c(Condition, Run), 
                                        values_from = Adj_Abundance,
                                        names_sort = TRUE))
  rownames(input_adj) <- input_adj$PTM
  input_adj <- input_adj %>% select(-PTM)

  ## Create contrast matrix
  class <- c()
  for (x in seq_len(conditions)){
    cond <- rep(paste0("Condition_", as.character(x)), runs)
    class <- c(class, cond)
  }
  class <- as.factor(class)
  design <- model.matrix(~0+class)

  input.matrix <- as.matrix(input)
  input_adj.matrix <- as.matrix(input_adj)

  ## Run models
  fit <- lmFit(input.matrix, design=design)
  fit_no_adj <- lmFit(input_adj.matrix, design=design)

  contrast.matrix <- design.pairs(colnames(design))
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)

  fit2_no_adj <- contrasts.fit(fit_no_adj, contrast.matrix)
  fit2_no_adj <- eBayes(fit2_no_adj)

  ## Output to data.table
  comparisons <- data.table()
  comparisons_no_adj <- data.table()

  for (g in seq_along(colnames(fit2$coefficients))){
    comparisons <- rbindlist(list(comparisons,
                                   data.table(PTM = rownames(fit2$coefficients),
                              Label = colnames(fit2$coefficients)[g],
                              Log2FC = as.vector(fit2$coefficients[,g]),
                              pvalue = as.vector(fit2$p.value[,g]),
                              df = as.vector(fit2$df.residual),
                              se = as.vector(fit2$sigma))))
    comparisons_no_adj <- rbindlist(list(comparisons_no_adj,
                            data.table(PTM = rownames(fit2_no_adj$coefficients),
                               Label = colnames(fit2_no_adj$coefficients)[g],
                               Log2FC = as.vector(fit2_no_adj$coefficients[,g]),
                               pvalue = as.vector(fit2_no_adj$p.value[,g]),
                               df = as.vector(fit2_no_adj$df.residual),
                               se = as.vector(fit2_no_adj$sigma))))
  }

  return(list(limma_test <- comparisons, 
              limma_no_adj_test = comparisons_no_adj))
}

## Simulation Start ------------------------------------------------------------
make_simulations <- function(missing_data, ptm_features){
  ## Specify types of data to simulate
  ## variance
  s <- c(.2,.3)
  ## Bioreps
  reps <- c(2,3,5,10)
  ## Conditions
  cond <- c(2,3,4)
  
  param_combos <- expand.grid(s, reps, cond)
  all_data <- list()
  i <- 1
  
  for (row in seq_len(nrow(param_combos))){
  
    ## Change in conditions
    if (param_combos[row, 3] == 2){
      del_arr <- c(0, .75)
      del_adj <- c(0, .25)
      del_arr_no_change <- c(0, 0)
    } else if (param_combos[row, 3] == 3) {
      del_arr <- c(0, .75, 1.5)
      del_adj <- c(0, .25, 1.25)
      del_arr_no_change <- c(0, 0, 0)
    } else if (param_combos[row, 3] == 4) {
      del_arr <- c(0, .75, 1.5, 2.25)
      del_adj <- c(0, .25, 1.25, 2.25)
      del_arr_no_change <- c(0, 0, 0, 0)
    }
  
    if (param_combos[row, 1] == .2){
      prot_var <- .2
    } else {
      prot_var <- .3
    }
  
    ## Sim
    sim <- PTMsimulateExperiment(
      nGroup=param_combos[row, 3], nRep=param_combos[row, 2], nProtein=500, 
      nSite=1, nFeature=ptm_features, nFeature_prot = 10,
      logAbundance=list(
        PTM=list(mu=25, delta = del_arr, sRep=param_combos[row, 1], sPeak=.25),
        PROTEIN=list(mu=25, delta = del_arr_no_change, sRep=param_combos[row, 1], 
                     sPeak=0.25))
    )
    sim_no_change1 <- PTMsimulateExperiment(
      nGroup=param_combos[row, 3], nRep=param_combos[row, 2], nProtein=250, 
      nSite=1, nFeature=ptm_features, nFeature_prot = 10,
      logAbundance=list(
        PTM=list(mu=25, delta = del_arr, sRep=param_combos[row, 1], sPeak=0.25),
        PROTEIN=list(mu=25, delta = del_arr, sRep=param_combos[row, 1], 
                     sPeak=0.25))
    )
    sim_no_change2 <- PTMsimulateExperiment(
      nGroup=param_combos[row, 3], nRep=param_combos[row, 2], nProtein=250, 
      nSite=1, nFeature=ptm_features, nFeature_prot = 10,
      logAbundance=list(
        PTM=list(mu=25, delta = del_arr_no_change, sRep=param_combos[row, 1], 
                 sPeak=0.25),
        PROTEIN=list(mu=25, delta = del_arr_no_change, sRep=param_combos[row, 1], 
                     sPeak=0.25))
    )
    sim_no_change1$PTM$protein <- paste0(sim_no_change1$PTM$protein, 
                                         "|NoChange1")
    sim_no_change1$PROTEIN$protein <- paste0(sim_no_change1$PROTEIN$protein, 
                                             "|NoChange1")
    sim_no_change2$PTM$protein <- paste0(sim_no_change2$PTM$protein, "|NoChange2")
    sim_no_change2$PROTEIN$protein <- paste0(sim_no_change2$PROTEIN$protein, 
                                             "|NoChange2")
  
    sim_no_change_ptm <- rbindlist(list(sim_no_change1$PTM, sim_no_change2$PTM))
    sim_no_change_prot <- rbindlist(list(sim_no_change1$PROTEIN, 
                                         sim_no_change2$PROTEIN))
    sim_no_change <- list(PTM = sim_no_change_ptm, PROTEIN = sim_no_change_prot)
  
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
  
  ## Remove random rows
  if (missing_data){
    for (d in seq_along(all_data)){
      
      dump_ptm <- sample(1:nrow(all_data[[d]]$PTM), 
                         nrow(all_data[[d]]$PTM) * 2)
      dump_prot <- sample(1:nrow(all_data[[d]]$PROTEIN), 
                          nrow(all_data[[d]]$PROTEIN) * 2)
      
      all_data[[d]]$PTM[dump_ptm]$Intensity <- NA
      all_data[[d]]$PROTEIN[dump_prot]$Intensity <- NA
    }
  }
  return(all_data)
}


## Simulation 1 ----------------------------------------------------------------
## Gen data
all_data <- make_simulations(FALSE, 10)

ptm_models <- list()
adjusted_models <- list()
ttest <- list()
adj_ttest <- list()
limma_results <- list()
adj_limma <- list()

ttest_msstats <- list()
adj_ttest_msstats <- list()
limma_msstats <- list()
adj_limma_msstats <- list()

## Loop to summarize and model each dataset
for (i in seq_along(all_data)){

  ## Run MSstatsPTM
  temp_sum <- dataSummarizationPTM(all_data[[i]], normalization.PTM = FALSE, 
                                   normalization = FALSE,
                                   MBimpute = FALSE, MBimpute.PTM = FALSE)
  temp_model <- groupComparisonPTM(temp_sum, data.type = "LabelFree")
  ptm_models[[i]] <- temp_model$PTM.Model
  temp_model$ADJUSTED.Model <- temp_model$ADJUSTED.Model %>% 
    filter(!is.na(Protein))
  adjusted_models[[i]] <- temp_model$ADJUSTED.Model

  ## Run ttest and limma
  ## Merge datasets
  ptm_df <- all_data[[i]]$PTM
  protein_df <- all_data[[i]]$PROTEIN

  ## Summarize using log sum of runs
  summarized_ptm <- data.table()
  summarized_proteins <- data.table()
  runs <- unique(ptm_df$Run)
  for (r in seq_along(runs)){
    sum_runs <- ptm_df %>% filter(Run == runs[[r]]) %>%
      group_by(ProteinName, Condition, Run) %>%
      summarize(Abundance = log2(sum(Intensity)))
    summarized_ptm <- rbindlist(list(summarized_ptm, sum_runs))

    sum_runs_prot <- protein_df %>% filter(Run == runs[[r]]) %>%
      group_by(ProteinName, Condition, Run) %>%
      summarize(Abundance = log2(sum(Intensity, na.rm=TRUE)))
    summarized_proteins <- rbindlist(list(summarized_proteins, sum_runs_prot))
  }

  ## Merge data
  summarized_ptm$PTM <- summarized_ptm$ProteinName
  summarized_ptm$ProteinName <- sapply(summarized_ptm$PTM, function(x) {paste(
    str_split(x, "_",3)[[1]][1:2], collapse = "_")})
  joined <- merge(summarized_ptm, summarized_proteins, 
                  by = c("ProteinName", "Run", "Condition"), all.x = TRUE)
  joined$Adj_Abundance <- joined$Abundance.x - joined$Abundance.y

  ## Run ttest
  ttest_list <- run_ttest(joined)

  ttest[[i]] <- ttest_list$ttest_temp
  adj_ttest[[i]] <- ttest_list$ttest_adj_temp

  ## Run Limma
  limma_test_res <- fit_limma(joined, param_combos[i, 3], param_combos[i, 2])
  limma_results[[i]] <- limma_test_res[[1]]
  adj_limma[[i]] <- limma_test_res$limma_no_adj_test

  print(paste0("Dataset ", as.character(i), " completed"))
}

ptm_models_sim1 <- ptm_models
adjusted_models_sim1 <- adjusted_models
ttest_sim1 <- ttest
adj_ttest_sim1 <- adj_ttest
limma_results_sim1 <- limma_results
adj_limma_sim1 <- adj_limma

# save(ptm_models_sim1, file = "../data/ptm_models_sim1.rda")
# save(adjusted_models_sim1, file = "../data/adjusted_models_sim1.rda")
# save(ttest_sim1, file = "../data/ttest_models_sim1.rda")
# save(adj_ttest_sim1, file = "../data/adj_ttest_models_sim1.rda")
# save(limma_results_sim1, file = "../data/limma_models_sim1.rda")
# save(adj_limma_sim1, file = "../data/adj_limma_models_sim1.rda")

## Simulation 2 ----------------------------------------------------------------
## Gen data
all_data <- make_simulations(TRUE, 2)

ptm_models <- list()
adjusted_models <- list()
ttest <- list()
adj_ttest <- list()
limma_results <- list()
adj_limma <- list()

ttest_msstats <- list()
adj_ttest_msstats <- list()
limma_msstats <- list()
adj_limma_msstats <- list()

## Loop to summarize and model each dataset
for (i in seq_along(all_data)){
  
  ## Run MSstatsPTM
  temp_sum <- dataSummarizationPTM(all_data[[i]], normalization.PTM = FALSE, 
                                   normalization = FALSE,
                                   MBimpute = FALSE, MBimpute.PTM = FALSE)
  temp_model <- groupComparisonPTM(temp_sum, data.type = "LabelFree")
  ptm_models[[i]] <- temp_model$PTM.Model
  temp_model$ADJUSTED.Model <- temp_model$ADJUSTED.Model %>% 
    filter(!is.na(Protein))
  adjusted_models[[i]] <- temp_model$ADJUSTED.Model
  
  ## Run ttest and limma
  ## Merge datasets
  ptm_df <- all_data[[i]]$PTM
  protein_df <- all_data[[i]]$PROTEIN
  
  ## Summarize using log sum of runs
  summarized_ptm <- data.table()
  summarized_proteins <- data.table()
  runs <- unique(ptm_df$Run)
  for (r in seq_along(runs)){
    sum_runs <- ptm_df %>% filter(Run == runs[[r]]) %>%
      group_by(ProteinName, Condition, Run) %>%
      summarize(Abundance = log2(sum(Intensity)))
    summarized_ptm <- rbindlist(list(summarized_ptm, sum_runs))
    
    sum_runs_prot <- protein_df %>% filter(Run == runs[[r]]) %>%
      group_by(ProteinName, Condition, Run) %>%
      summarize(Abundance = log2(sum(Intensity, na.rm=TRUE)))
    summarized_proteins <- rbindlist(list(summarized_proteins, sum_runs_prot))
  }
  
  ## Merge data
  summarized_ptm$PTM <- summarized_ptm$ProteinName
  summarized_ptm$ProteinName <- sapply(summarized_ptm$PTM, function(x) {paste(
    str_split(x, "_",3)[[1]][1:2], collapse = "_")})
  joined <- merge(summarized_ptm, summarized_proteins, 
                  by = c("ProteinName", "Run", "Condition"), all.x = TRUE)
  joined$Adj_Abundance <- joined$Abundance.x - joined$Abundance.y
  
  ## Run ttest
  ttest_list <- run_ttest(joined)
  
  ttest[[i]] <- ttest_list$ttest_temp
  adj_ttest[[i]] <- ttest_list$ttest_adj_temp
  
  ## Run Limma
  limma_test_res <- fit_limma(joined, param_combos[i, 3], param_combos[i, 2])
  limma_results[[i]] <- limma_test_res[[1]]
  adj_limma[[i]] <- limma_test_res$limma_no_adj_test
  
  print(paste0("Dataset ", as.character(i), " completed"))
}

ptm_models_sim2 <- ptm_models
adjusted_models_sim2 <- adjusted_models
ttest_sim2 <- ttest
adj_ttest_sim2 <- adj_ttest
limma_results_sim2 <- limma_results
adj_limma_sim2 <- adj_limma

# save(ptm_models_sim2, file = "../data/ptm_models_sim2.rda")
# save(adjusted_models_sim2, file = "../data/adjusted_models_sim2.rda")
# save(ttest_sim2, file = "../data/ttest_models_sim2.rda")
# save(adj_ttest_sim2, file = "../data/adj_ttest_models_sim2.rda")
# save(limma_results_sim2, file = "../data/limma_models_sim2.rda")
# save(adj_limma_sim2, file = "../data/adj_limma_models_sim2.rda")
