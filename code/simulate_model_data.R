
## Load packages
library(MSstatsPTM)
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

## Run anova using run level summared data
fit_anova = function(data, contrast){

  ptms = unique(data$PTM)
  test_results = data.table()
  unadj_test_results = data.table()
  for (i in seq_along(ptms)){#
    # print(ptms[[i]])
    ptm_data_adj = data %>% filter(PTM == ptms[[i]] & is.finite(Adj_Abundance))
    ptm_data_unadj = data %>% filter(PTM == ptms[[i]] & is.finite(Abundance.x))

    temp_contrast = contrast


    tryCatch({

      ## Adj model
      res_aov_adj = lm(Adj_Abundance~Condition, data = ptm_data_adj)
      model_summary_adj = summary(res_aov_adj)

      coefs_adj = res_aov_adj$coefficients
      coefs_adj[1] = 0
      df_adj = res_aov_adj$df.residual
      sigma_adj = model_summary_adj$sigma

      ## unadj model
      res_aov_unadj = lm(Abundance.x~Condition, data = ptm_data_unadj)
      model_summary_unadj = summary(res_aov_unadj)

      coefs_unadj = res_aov_unadj$coefficients
      coefs_unadj[1] = 0
      df_unadj = res_aov_unadj$df.residual
      sigma_unadj = model_summary_unadj$sigma

      for (j in seq_len(nrow(temp_contrast))){

        comp = temp_contrast[j,]
        label = rownames(temp_contrast)[j]

        ## Adj
        beta_adj = sum(coefs_adj*comp)
        tstat_adj = beta_adj / sigma_adj

        pval_adj = 2*pt(-abs(tstat_adj),df=df_adj)

        test_results = rbindlist(list(test_results,
                                      data.table(ptm = ptms[[i]], label = label, log2FC = beta_adj,
                                                 se = sigma_adj, tstat = tstat_adj, pvalue = pval_adj)))

        ## unadj
        beta_unadj = sum(coefs_unadj*comp)
        tstat_unadj = beta_unadj / sigma_unadj

        pval_unadj = 2*pt(-abs(tstat_unadj),df=df_unadj)

        unadj_test_results = rbindlist(list(unadj_test_results,
                                            data.table(ptm = ptms[[i]], label = label, log2FC = beta_unadj,
                                                       se = sigma_unadj, tstat = tstat_unadj, pvalue = pval_unadj)))

        }

    },
    error=function(e){cat("ERROR : fitting error", as.character(i))})
  }

  return(list(anova_unadj = unadj_test_results, anova_adj = test_results))
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

## Fit limma model (both not and with adjust) for given dataset
fit_limma <- function(summarized_data, conditions, runs){

  ## Convert into format required for limma
  input <- data.frame(summarized_data %>% select(PTM, Run,
                                                 Condition, Abundance.x) %>%
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
  # conditions = length(colnames(contrast))
  for (x in seq_len(conditions)){
    cond <- rep(paste0("G_", as.character(x)), runs)
    class <- c(class, cond)
  }
  class <- as.factor(class)
  design <- model.matrix(~0+class)
  # colnames(design) <- colnames(contrast)

  input.matrix <- as.matrix(input)
  input_adj.matrix <- as.matrix(input_adj)

  ## Run models
  fit <- lmFit(input.matrix, design=design)
  fit_adj <- lmFit(input_adj.matrix, design=design)

  contrast.matrix <- design.pairs(colnames(design))
  # contrast.matrix <- makeContrasts(mix2-mix1, mix3-mix1,
  #                                  mix4-mix1, mix3-mix2,
  #                                  mix4-mix2, mix4-mix3, levels=design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)

  fit2_adj <- contrasts.fit(fit_adj, contrast.matrix)
  fit2_adj <- eBayes(fit2_adj)

  ## Output to data.table
  comparisons <- data.table()
  comparisons_adj <- data.table()

  for (g in seq_along(colnames(fit2$coefficients))){
    comparisons <- rbindlist(list(comparisons,
                                  data.table(PTM = rownames(fit2$coefficients),
                                             Label = colnames(fit2$coefficients)[g],
                                             Log2FC = as.vector(fit2$coefficients[,g]),
                                             pvalue = as.vector(fit2$p.value[,g]),
                                             df = as.vector(fit2$df.residual),
                                             se = as.vector(fit2$sigma))))
    comparisons_adj <- rbindlist(list(comparisons_adj,
                                      data.table(PTM = rownames(fit2_adj$coefficients),
                                                 Label = colnames(fit2_adj$coefficients)[g],
                                                 Log2FC = as.vector(fit2_adj$coefficients[,g]),
                                                 pvalue = as.vector(fit2_adj$p.value[,g]),
                                                 df = as.vector(fit2_adj$df.residual),
                                                 se = as.vector(fit2_adj$sigma))))
  }

  return(list(limma_test = comparisons, limma_adj_test = comparisons_adj))
}

## Simulation 1 Start ----------------------------------------------------------
## Specify types of data to simulate
s <- c(.2,.3)
reps <- c(2,3,5,10)
cond <- c(2,3,4)

param_combos <- expand.grid(s, reps, cond)
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

## Run Models ------------------------------------------------------------------
ptm_models <- list()
adjusted_models <- list()
anova <- list()
adj_anova <- list()
limma_results <- list()
adj_limma <- list()

# Contrast
comparison1 = matrix(c(-1,1),nrow=1, byrow=TRUE)
row.names(comparison1) = c("G_2 vs G_1")
colnames(comparison1) = c("G_1", "G_2")

comparison2 = matrix(c(-1,1,0,
                       0,-1,1),nrow=2, byrow=TRUE)
row.names(comparison2) = c("G_2 vs G_1", "G_3 vs G_2")
colnames(comparison2) = c("G_1", "G_2", "G_3")

comparison3 = matrix(c(-1,1,0,0,
                       0,-1,1,0,
                       0,0,-1,1),nrow=3, byrow=TRUE)
row.names(comparison3) = c("G_2 vs G_1", "G_3 vs G_2",
                           "G_4 vs G_3")
colnames(comparison3) = c("G_1", "G_2", "G_3", "G_4")

## Loop to summarize and model each dataset
for (i in seq_along(all_data)){

  if (i <= 8){
    contrast = comparison1
  } else if (i <= 16){
    contrast = comparison2
  } else {
    contrast = comparison3
  }

  ## Run MSstatsPTM
  temp_sum <- dataSummarizationPTM(all_data[[i]], normalization.PTM = FALSE,
                                   normalization = FALSE, MBimpute = FALSE,
                                   MBimpute.PTM = FALSE)
  temp_model <- groupComparisonPTM(temp_sum, data.type = "LabelFree",
                                   contrast.matrix = contrast)
  ptm_models[[i]] <- temp_model$PTM.Model
  temp_model$ADJUSTED.Model <- temp_model$ADJUSTED.Model %>% filter(!is.na(Protein))
  adjusted_models[[i]] <- temp_model$ADJUSTED.Model

  ## Run anova and limma
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
  summarized_ptm$ProteinName <- sapply(summarized_ptm$PTM, function(x) {paste(str_split(x, "_",3)[[1]][1:2], collapse = "_")})
  joined <- merge(summarized_ptm, summarized_proteins, by = c("ProteinName", "Run", "Condition"), all.x = TRUE)
  joined$Adj_Abundance <- joined$Abundance.x - joined$Abundance.y

  ## Run anova
  anova_list <- fit_anova(joined, contrast)

  anova[[i]] <- anova_list$anova_unadj
  adj_anova[[i]] <- anova_list$anova_adj

  ## Run Limma
  limma_test_res <- fit_limma(joined, param_combos[i, 3], param_combos[i, 2])
  limma_test_res$limma_test$Log2FC = limma_test_res$limma_test$Log2FC*-1
  limma_results[[i]] <- limma_test_res$limma_test
  limma_test_res$limma_adj_test$Log2FC = limma_test_res$limma_adj_test$Log2FC*-1
  adj_limma[[i]] <- limma_test_res$limma_adj_test

  print(paste0("Dataset ", as.character(i), " completed"))
}

ptm_models_sim1 <- ptm_models
adjusted_models_sim1 <- adjusted_models
anova_sim1 <- anova
adj_anova_sim1 <- adj_anova
limma_results_sim1 <- limma_results
adj_limma_sim1 <- adj_limma

save(ptm_models_sim1, file = "../data/ptm_models_sim1.rda")
save(adjusted_models_sim1, file = "../data/adjusted_models_sim1.rda")
save(anova_sim1, file = "../data/anova_models_sim1.rda")
save(adj_anova_sim1, file = "../data/adj_anova_models_sim1.rda")
save(limma_results_sim1, file = "../data/limma_models_sim1.rda")
save(adj_limma_sim1, file = "../data/adj_limma_models_sim1.rda")


## Simulation 2 Start ----------------------------------------------------------
## Specify types of data to simulate
s <- c(.2,.3)
reps <- c(2,3,5,10)
cond <- c(2,3,4)

param_combos <- expand.grid(s, reps, cond)
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
  sim <- PTMsimulateExperiment(
    nGroup=param_combos[row, 3], nRep=param_combos[row, 2], nProtein=250, nSite=1, nFeature=2, nFeature_prot = 10,
    logAbundance=list(
      PTM=list(mu=25, delta = del_arr, sRep=param_combos[row, 1], sPeak=.25),#0.05),
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
for (d in seq_along(all_data)){

  dump_ptm <- sample(1:nrow(all_data[[d]]$PTM), nrow(all_data[[d]]$PTM) * .2)
  dump_prot <- sample(1:nrow(all_data[[d]]$PROTEIN), nrow(all_data[[d]]$PROTEIN) * .2)

  all_data[[d]]$PTM[dump_ptm]$Intensity <- NA
  all_data[[d]]$PROTEIN[dump_prot]$Intensity <- NA

}

simulation2_data = all_data
save(simulation2_data, file = "../data/simulation2_data.rda")

## Run Models ------------------------------------------------------------------
ptm_models <- list()
adjusted_models <- list()
anova <- list()
adj_anova <- list()
limma_results <- list()
adj_limma <- list()

# Contrast
comparison1 = matrix(c(-1,1),nrow=1, byrow=TRUE)
row.names(comparison1) = c("G_2 vs G_1")
colnames(comparison1) = c("G_1", "G_2")

comparison2 = matrix(c(-1,1,0,
                       0,-1,1),nrow=2, byrow=TRUE)
row.names(comparison2) = c("G_2 vs G_1", "G_3 vs G_2")
colnames(comparison2) = c("G_1", "G_2", "G_3")

comparison3 = matrix(c(-1,1,0,0,
                       0,-1,1,0,
                       0,0,-1,1),nrow=3, byrow=TRUE)
row.names(comparison3) = c("G_2 vs G_1", "G_3 vs G_2",
                           "G_4 vs G_3")
colnames(comparison3) = c("G_1", "G_2", "G_3", "G_4")

## Loop to summarize and model each dataset
for (i in seq_along(all_data)){

  if (i <= 8){
    contrast = comparison1
  } else if (i <= 16){
    contrast = comparison2
  } else {
    contrast = comparison3
  }

  ## Run MSstatsPTM
  temp_sum <- dataSummarizationPTM(all_data[[i]], normalization.PTM = FALSE,
                                   normalization = FALSE, MBimpute = FALSE,
                                   MBimpute.PTM = FALSE)
  temp_model <- groupComparisonPTM(temp_sum, data.type = "LabelFree",
                                   contrast.matrix = contrast)
  ptm_models[[i]] <- temp_model$PTM.Model
  temp_model$ADJUSTED.Model <- temp_model$ADJUSTED.Model %>% filter(!is.na(Protein))
  adjusted_models[[i]] <- temp_model$ADJUSTED.Model

  ## Run anova and limma
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
  summarized_ptm$ProteinName <- sapply(summarized_ptm$PTM, function(x) {paste(str_split(x, "_",3)[[1]][1:2], collapse = "_")})
  joined <- merge(summarized_ptm, summarized_proteins, by = c("ProteinName", "Run", "Condition"), all.x = TRUE)
  joined$Adj_Abundance <- joined$Abundance.x - joined$Abundance.y

  ## Run anova
  anova_list <- fit_anova(joined, contrast)

  anova[[i]] <- anova_list$anova_unadj
  adj_anova[[i]] <- anova_list$anova_adj

  ## Run Limma
  limma_test_res <- fit_limma(joined, param_combos[i, 3], param_combos[i, 2])
  limma_test_res$limma_test$Log2FC = limma_test_res$limma_test$Log2FC*-1
  limma_results[[i]] <- limma_test_res$limma_test
  limma_test_res$limma_adj_test$Log2FC = limma_test_res$limma_adj_test$Log2FC*-1
  adj_limma[[i]] <- limma_test_res$limma_adj_test

  print(paste0("Dataset ", as.character(i), " completed"))
}

ptm_models_sim2 <- ptm_models
adjusted_models_sim2 <- adjusted_models
anova_sim2 <- anova
adj_anova_sim2 <- adj_anova
limma_results_sim2 <- limma_results
adj_limma_sim2 <- adj_limma

save(ptm_models_sim2, file = "../data/ptm_models_sim2.rda")
save(adjusted_models_sim2, file = "../data/adjusted_models_sim2.rda")
save(anova_sim2, file = "../data/anova_models_sim2.rda")
save(adj_anova_sim2, file = "../data/adj_anova_models_sim2.rda")
save(limma_results_sim2, file = "../data/limma_models_sim2.rda")
save(adj_limma_sim2, file = "../data/adj_limma_models_sim2.rda")

