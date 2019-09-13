#!/usr/bin/env/ Rscript

# define functions of scoring indices

RMSE <- function(P, freq = rep(1,length(P))){ # freq of all 1 means that each element has equal probability to be taken during bootstrap, attribute required by boot
  RMSE <- sqrt(sum((sample(predicted.P, Nvar)-real.P)^2)/Nvar)
  return(RMSE)
}

MAE <- function(P, freq = rep(1,length(P))){ 
  MAE <- sum(abs(sample(predicted.P, Nvar)-real.P))/Nvar
  return(MAE)
}

PCC <- function(P, freq = rep(1,length(P))){ 
  PCC <- cor(sample(predicted.P, Nvar), real.P, method = "pearson" )
  return(PCC)
}

KCC <- function(P, freq = rep(1,length(P))){ 
  KCC <- cor(sample(predicted.P, Nvar), real.P, method = "kendall", use = "pairwise" )
  return(KCC)
}

Rsquared <- function(P, freq = rep(1,length(P))){ 
  Rsquared <- (cor(sample(predicted.P, Nvar), real.P, method = "pearson" ))^2
  return(Rsquared)
}




Bootstrap <- function(){
  #setwd("../")
  library(boot)
  results.table <- read.table("./results/paper/score.table_LM_completo.txt", header = T, sep ="\t", stringsAsFactors = F)
  predicted.table <- read.table("./results/paper/predicted_values.txt", header = T, sep ="\t", stringsAsFactors = F)
  # remove mutation and gene column
  real.data <- predicted.table[, c(3,4,5)]
  real.data <- cbind(real.data,real.data) # NB this reduce complexity of next for cycle
  predicted.table <- predicted.table[, (6:11)]
  P.bts.table <- results.table # create table to store bootstrap results
  # 
  nSamples <<- 10000
  for(i in (1:ncol(predicted.table))){
    real.P <<-  real.data[, i] [which(is.na(real.data[, i])==F)] # declare as global variable to use it inside boot function
    Nvar <<- length(real.P) # number of cells
    predicted.P <<- predicted.table[, i][which(is.na(real.data[, i])==F)] # remove predicted values without corresponding real ones
    PCC.rand <- boot(predicted.P, PCC, nSamples, parallel = "multicore", ncpus = 2)$t
    PCC.real <- results.table[i, 1]
    PCC.P <- sum( PCC.rand >= PCC.real)/nSamples # P of obtaining a value of PCC > than the original one,  after a random shuffle of predicted values
    P.bts.table[i, 1] <- round(PCC.P, digits = 2)
    
    KCC.rand <- boot(predicted.P, KCC, nSamples, parallel = "multicore", ncpus = 2)$t
    KCC.real <- results.table[i, 2]
    KCC.P <- sum( KCC.rand >= KCC.real)/nSamples
    P.bts.table[i, 2] <- round(KCC.P, digits = 2)
    
    MAE.rand <- boot(predicted.P, MAE, nSamples, parallel = "multicore", ncpus = 2)$t # compute bootstrap and obtain v
    MAE.real <- results.table[i, 3]
    MAE.P <- sum( MAE.rand <= MAE.real)/nSamples # P of obtaining a value of MAE < than the original one,  after a random shuffle of predicted values
    P.bts.table[i, 3] <- round(MAE.P, digits = 2)
    
    RMSE.rand <- boot(predicted.P, RMSE, nSamples, parallel = "multicore", ncpus = 2)$t # compute bootstrap and obtain v
    RMSE.real <- results.table[i, 4]
    RMSE.P <- sum( RMSE.rand <= RMSE.real)/nSamples # P of obtaining a value of RMSE < than the original one,  after a random shuffle of predicted values
    P.bts.table[i, 4] <- round(RMSE.P, digits = 2)
    
    Rsquared.rand <- boot(predicted.P, Rsquared, nSamples, parallel = "multicore", ncpus = 2)$t
    Rsquared.real <- results.table[i, 5]
    Rsquared.P <- sum( Rsquared.rand >= Rsquared.real)/nSamples # P of obtaining a value of Rsquared > than the original one,  after a random shuffle of predicted values
    P.bts.table[i, 5] <- round(Rsquared.P, digits = 2)
    write.table(P.bts.table, file='./results/paper/bootstrap_table_assessment.txt', sep = "\t" , quote = F, row.names = T)
  }
 
}

Bootstrap()
