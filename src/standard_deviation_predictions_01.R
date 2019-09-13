#!/usr/bin/env/ Rscript

if(F){ # if(TRUE) run only this script
  library(deSolve)
  library(pracma)
  library(ggplot2) # used to make plot
  library(gridExtra) #used to make plot
  setwd('../')
  source('src/optimization_fmut.R')
  source('src/Human_model_R_ORIGINAL_vDp_2012.R')
}


SDpredictions <- function(){
  T.S <- read.table(file = "data/Other_genes_database.txt", header=TRUE, sep = '\t') # training set
  fmut.table <- read.table(file = './results/fmut_optimized_Levenberg_Mar.txt', sep = '\t', header = T)
  fmut_v <- rep(1, 21)
  CL <<- VdP_model(fmut_v) # control cholesterol levels
  sd.table <- matrix(ncol = 4, nrow = 10) # table to store sd
  colnames(sd.table) <- c("GENE", "HDL", "LDL", "TC")
  
  
  # make function to compute SD
  Comp_SD <- function(fmut.table = fmut.table, gene = "LDLR", rate = c(5,7)) {
    mut_eff <- fmut.table$mean[grep(gene, as.character(fmut.table$GENE))]
    var.tab <- cbind(mut_eff, rep(NA, length(mut_eff)), rep(NA, length(mut_eff)),  rep(NA, length(mut_eff)))
    fmut_v <- rep(1, 21)
    for (i in (1:nrow(var.tab))) {
      fmut_v[rate] <- mut_eff[i]
      mut_chol <- VdP_model(fmut_v)
      var.tab[i, c(2,3,4)] <- c(mut_chol/CL, (sum(mut_chol)/sum(CL))) # [HDL, LDL]/[HDL_control, LDL_control], TC/TC_control
    }
    return(c(round(sd(var.tab[, 2]), digits = 2), round(sd(var.tab[, 3]), digits = 2), round(sd(var.tab[, 4]), digits = 2)))
  }
  
  # LDLR
  sd.table[1, ] <- c("LDLR", Comp_SD(fmut.table = fmut.table, gene = "LDLR", rate = c(5,7)))
  
  # APOB HET
  sd.table[2, ] <- c("APOB het", Comp_SD(fmut.table = fmut.table, gene = "^APOB$", rate = c(5,7)))
  
  # APOB HOM
  sd.table[3, ] <- c("APOB hom", 0, 0, 0) # only one value in training set
  
  # ABCA1
  sd.table[4, ] <- c("ABCA1", Comp_SD(fmut.table = fmut.table, gene = "ABCA1", rate = c(8, 16, 17)))
  
  # APOE
  sd.table[5, ] <- c("APOE", Comp_SD(fmut.table = fmut.table, gene = "APOE", rate = c(5, 7)))
  
  # CETP
  sd.table[6, ] <- c("CETP", 0, 0, 0)
  
  # LCAT het
  sd.table[7, ] <- c("LCAT het", Comp_SD(fmut.table = fmut.table, gene = "^LCAT$", rate = c(9)))
  
  # LCAT hom
  sd.table[8, ] <- c("LCAT hom", Comp_SD(fmut.table = fmut.table, gene = "LCAT_HOM", rate = c(9)))
  
  # DHCR7
  table.dhcr <- read.table(file = "results/fmut_optimized_Levenberg_Mar_TC.txt", sep = "\t", header = T)
  colnames(table.dhcr) <- c("GENE", "mean")
  table.dhcr[, 1] <- rep("DHCR7", nrow(table.dhcr))
  sd.table[9, ] <- c("DHCR7", Comp_SD(fmut.table = table.dhcr, gene = "DHCR7", rate = c(1, 2, 3)))
  
  # CYP7A1
  sd.table[10, ] <- c("CYP7A1", Comp_SD(fmut.table = fmut.table, gene = "CYP7A1", rate = c(18)))
  
  write.table(file = './results/paper/sd_model_predictions.txt', sd.table, sep = '\t', row.names = F, col.names = T, quote = F)
}


