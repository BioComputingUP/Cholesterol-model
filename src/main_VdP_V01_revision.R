#!/usr/bin/env/ Rscript

library(deSolve)
library(pracma)
library(ggplot2) # used to make plot
library(gridExtra) #used to make plot

setwd('../')

source('src/optimization_fmut.R')
source('src/Human_model_R_ORIGINAL_vDp_2012.R')
source('src/make_simulations_V01_revision.R')
source('src/create_boxplot_patients.R')
source('src/bootstrap_02_revision.R')
source('src/standard_deviation_predictions_01.R')

T.S <- read.table(file = "data/Other_genes_database.txt", header=TRUE, sep = '\t') # training set

if(FALSE){ # training phase
  LDL.fmut <- Compute_f_mut_LDL() # compute optimal fmut considering HDL, LDL level 
  fmut.table <- as.data.frame(cbind(LDL.fmut, as.character(T.S$GENE), as.character(T.S$Mutation)), stringsAsFactors = F)
  fmut.table[, 1] <- as.numeric(LDL.fmut)
  colnames(fmut.table) <- c("mean", "GENE", "mutation")
  write.table(file = './results/fmut_optimized_Levenberg_Mar.txt', fmut.table, sep = '\t', row.names = F, col.names = T, quote = F)
  if(FALSE){
    TC.fmut <- Compute_f_mut_TC() # compute optimal fmut considering TC level, only done with DHCR7 fmut (1,2,3)
    fmut.TC.results <- as.data.frame(cbind(seq(1:length(TC.fmut)),TC.fmut))
    colnames(fmut.TC.results) <- c("ID", "fmut")
    write.table(file = './results/fmut_optimized_Levenberg_Mar_TC.txt', fmut.TC.results, sep = '\t', row.names = F, col.names = T, quote = F)
    }
  }else{
  # load trained pararameters
  fmut.table <- read.table(file = './results/fmut_optimized_Levenberg_Mar.txt', sep = '\t', header = T)
}

# create matrix to store validation results
HDL.matrix <- matrix("empty", nrow = 9, ncol = 3)
colnames(HDL.matrix) <- c("Nomi", "type", "mean")
LDL.matrix <- HDL.matrix
TC.matrix <- HDL.matrix


#### create table with abundance of patients fo each gene 
pos.g.tab <- cbind(c(4, 2, 3, 5, 6, 10, 7, 8, 1), names(summary(T.S$GENE)), as.numeric(summary(T.S$GENE)))
pos.g.tab <- pos.g.tab[order(as.numeric(pos.g.tab[, 1])), ]
n.mut <- sapply(1:nrow(pos.g.tab), function(x){length(unique(T.S$Variant[grep(paste("^", pos.g.tab[x, 2], "$", sep = ""), T.S$GENE)])) } )
pos.g.tab <- cbind(pos.g.tab, n.mut)
colnames(pos.g.tab) <-c("ID", "GENE", "# Patients", "# Mutations")

write.table(file = './results/paper/number_pat_gene.txt', pos.g.tab, sep = '\t', row.names = F, col.names = T, quote = F)

# check predictions with VDP dataset
VdP_table <- read.table(file = "data/VdP_results_table.txt", header = TRUE, sep = '\t') # load vdP table 4
fmut <- rep(1, 21) # fmut of selected rate 1, 2, 3, 5, 7, 8, 9, 16, 17, 18, 21 
Control_Chol_levels <- VdP_model(fmut)







#LDLR FH part

LDLR_effect <- mean(fmut.table$mean[grep("LDLR", as.character(fmut.table$GENE))])
fmut_FH <- fmut
fmut_FH[c(5,7)]<- LDLR_effect
results_TL <- VdP_model(fmut_FH)
fmut_FH[c(5,7)] <- VdP_table$fmut[1]
HDL.matrix[2, ] <- c('Predicted (genotype)', '1', results_TL[1]/Control_Chol_levels[1])
LDL.matrix[2, ] <- c('Predicted (genotype)', '1', results_TL[2]/Control_Chol_levels[2])
TC.matrix[2, ] <- c('Predicted (genotype)', '1', sum(results_TL)/sum(Control_Chol_levels))
results_VdP <- VdP_model(fmut_FH)
HDL.matrix[3, ] <- c('Predicted (van de Pas)', '1', results_VdP[1]/Control_Chol_levels[1])
LDL.matrix[3, ] <- c('Predicted (van de Pas)', '1', results_VdP[2]/Control_Chol_levels[2])
TC.matrix[3, ] <- c('Predicted (van de Pas)', '1', sum(results_VdP)/sum(Control_Chol_levels))
HDL.matrix[1, ] <- c('Experimental', '1', VdP_table$HDL.C..Relative.to.Control.[1])
LDL.matrix[1, ] <- c('Experimental', '1', VdP_table$non.HDL.C.Relative.to.Control.[1])
TC.matrix[1, ] <- c('Experimental', '1', VdP_table$TC.Relative.to.Control.[1])

# ApoB_het part

ApoB_Het <- mean(fmut.table$mean[grep("^APOB$", as.character(fmut.table$GENE))])
fmut_APOB <- fmut
fmut_APOB[c(5,7)] <- ApoB_Het
results_TL <- VdP_model(fmut_APOB)
HDL.matrix[5, ] <- c('Predicted (genotype)', '2',results_TL[1]/Control_Chol_levels[1])
LDL.matrix[5, ] <- c('Predicted (genotype)', '2', results_TL[2]/Control_Chol_levels[2])
TC.matrix[5, ] <- c('Predicted (genotype)', '2', sum(results_TL)/sum(Control_Chol_levels))
fmut_APOB[c(5,7)] <- VdP_table$fmut[2]
results_VdP <- VdP_model(fmut_APOB)
HDL.matrix[6, ] <- c('Predicted (van de Pas)', '2',results_VdP[1]/Control_Chol_levels[1])
LDL.matrix[6, ] <- c('Predicted (van de Pas)', '2', results_VdP[2]/Control_Chol_levels[2])
TC.matrix[6, ] <- c('Predicted (van de Pas)', '2', sum(results_VdP)/sum(Control_Chol_levels))
HDL.matrix[4, ] <- c('Experimental', '2', VdP_table$HDL.C..Relative.to.Control.[2])
LDL.matrix[4, ] <- c('Experimental', '2', VdP_table$non.HDL.C.Relative.to.Control.[2])
TC.matrix[4, ] <- c('Experimental', '2', VdP_table$TC.Relative.to.Control.[2])

# ApoB_hom part

fmut_APOB_hom <- fmut
APOB_hom <- fmut.table$mean[grep("APOB_HOM", as.character(fmut.table$GENE))] 
fmut_APOB_hom[c(5,7)]<- APOB_hom # use previuosly computed fmut APOB effect
results_TL <- VdP_model(fmut_APOB_hom)
HDL.matrix[8, ] <- c('Predicted (genotype)', '3', results_TL[1]/Control_Chol_levels[1])
LDL.matrix[8, ] <- c('Predicted (genotype)', '3', results_TL[2]/Control_Chol_levels[2])
TC.matrix[8, ] <- c('Predicted (genotype)', '3', sum(results_TL)/sum(Control_Chol_levels))
fmut_APOB_hom[c(5,7)] <- VdP_table$fmut[2]
results_VdP <- VdP_model(fmut_APOB_hom)
HDL.matrix[9, ] <- c('Predicted (van de Pas)', '3', results_VdP[1]/Control_Chol_levels[1])
LDL.matrix[9, ] <- c('Predicted (van de Pas)', '3', results_VdP[2]/Control_Chol_levels[2])
TC.matrix[9, ] <- c('Predicted (van de Pas)', '3', sum(results_VdP)/sum(Control_Chol_levels))
HDL.matrix[7, ] <- c('Experimental', '3',  VdP_table$HDL.C..Relative.to.Control.[3])
LDL.matrix[7, ] <- c('Experimental', '3', VdP_table$non.HDL.C.Relative.to.Control.[3])
TC.matrix[7, ] <- c('Experimental', '3', VdP_table$TC.Relative.to.Control.[3])
HDL.matrixI <- HDL.matrix
LDL.matrixI <- LDL.matrix
TC.matrixI <- TC.matrix


# start second part with other variants

HDL.matrix <- matrix("empty", nrow = (nrow(VdP_table)*3), ncol = 3)
colnames(HDL.matrix) <- c("Nomi", "type", "mean")
LDL.matrix <- HDL.matrix
TC.matrix <- HDL.matrix

# ABCA1 validation

fmut_ABCA1 <- fmut # fmut of selected rate 1, 2, 3, 5, 7, 8, 9, 16, 17, 18, 21 
ABCA1_mut <- mean(fmut.table$mean[grep("ABCA1", as.character(fmut.table$GENE))])
fmut_ABCA1[c(8,16,17)] <- ABCA1_mut  # I use whole mean of fmut computed in training set because I don't know mutations of patients used to validate vdP model in the original publication
results_TL <- VdP_model(fmut_ABCA1)
HDL.matrix[2, ] <- c('Predicted (genotype)', VdP_table$No.[4], results_TL[1]/Control_Chol_levels[1])
LDL.matrix[2, ] <- c('Predicted (genotype)', VdP_table$No.[4], results_TL[2]/Control_Chol_levels[2])
TC.matrix[2, ] <- c('Predicted (genotype)', VdP_table$No.[4], sum(results_TL)/sum(Control_Chol_levels))
fmut_ABCA1[c(8,16,17)] <- VdP_table$fmut[4]
results_VdP <- VdP_model(fmut_ABCA1)
HDL.matrix[3, ] <- c('Predicted (van de Pas)', VdP_table$No.[4], results_VdP[1]/Control_Chol_levels[1])
LDL.matrix[3, ] <- c('Predicted (van de Pas)', VdP_table$No.[4], results_VdP[2]/Control_Chol_levels[2])
TC.matrix[3, ] <- c('Predicted (van de Pas)', VdP_table$No.[4], sum(results_VdP)/sum(Control_Chol_levels))
HDL.matrix[1, ] <- c('Experimental', VdP_table$No.[4], VdP_table$HDL.C..Relative.to.Control.[4])
LDL.matrix[1, ] <- c('Experimental', VdP_table$No.[4], VdP_table$non.HDL.C.Relative.to.Control.[4])
TC.matrix[1, ] <- c('Experimental', VdP_table$No.[4], VdP_table$TC.Relative.to.Control.[4])

# APOE validation

fmut_APOE <- fmut # fmut of selected rate 1, 2, 3, 5, 7, 8, 9, 16, 17, 18, 21 
APOE_mut <-  mean(fmut.table$mean[grep("APOE", as.character(fmut.table$GENE))])
fmut_APOE[c(5,7)] <- APOE_mut # I use whole mean of fmut computed in training set because type of mutations where not specified in original paper
results_TL <- VdP_model(fmut_APOE)
HDL.matrix[5, ] <- c('Predicted (genotype)', VdP_table$No.[5], results_TL[1]/Control_Chol_levels[1])
LDL.matrix[5, ] <- c('Predicted (genotype)', VdP_table$No.[5], results_TL[2]/Control_Chol_levels[2])
TC.matrix[5, ] <- c('Predicted (genotype)', VdP_table$No.[5], sum(results_TL)/sum(Control_Chol_levels))
fmut_APOE[c(5,7)] <- VdP_table$fmut[5]
results_VdP <- VdP_model(fmut_APOE)
HDL.matrix[6, ] <- c('Predicted (van de Pas)', VdP_table$No.[5], results_VdP[1]/Control_Chol_levels[1])
LDL.matrix[6, ] <- c('Predicted (van de Pas)', VdP_table$No.[5], results_VdP[2]/Control_Chol_levels[2])
TC.matrix[6, ] <- c('Predicted (van de Pas)', VdP_table$No.[5], sum(results_VdP)/sum(Control_Chol_levels))
HDL.matrix[4, ] <- c('Experimental', VdP_table$No.[5], VdP_table$HDL.C..Relative.to.Control.[5])
LDL.matrix[4, ] <- c('Experimental', VdP_table$No.[5], VdP_table$non.HDL.C.Relative.to.Control.[5])
TC.matrix[4, ] <- c('Experimental', VdP_table$No.[5], VdP_table$TC.Relative.to.Control.[5])

# CETP validation

fmut_CETP <- fmut # fmut of selected rate 1, 2, 3, 5, 7, 8, 9, 16, 17, 18, 21 
CETP_mut <- fmut.table$mean[grep("CETP_D442G", as.character(fmut.table$GENE))]
fmut_CETP[21] <-CETP_mut #  D442G mutation only people het for exon 15 mutations are considered in vdP original model validation
results_TL <- VdP_model(fmut_CETP)
HDL.matrix[8, ] <- c('Predicted (genotype)', VdP_table$No.[6], results_TL[1]/Control_Chol_levels[1])
LDL.matrix[8, ] <- c('Predicted (genotype)', VdP_table$No.[6], results_TL[2]/Control_Chol_levels[2])
TC.matrix[8, ] <- c('Predicted (genotype)', VdP_table$No.[6], sum(results_TL)/sum(Control_Chol_levels))
fmut_CETP[21] <- VdP_table$fmut[6]
results_VdP <- VdP_model(fmut_CETP)
HDL.matrix[9, ] <- c('Predicted (van de Pas)', VdP_table$No.[6], results_VdP[1]/Control_Chol_levels[1])
LDL.matrix[9, ] <- c('Predicted (van de Pas)', VdP_table$No.[6], results_VdP[2]/Control_Chol_levels[2])
TC.matrix[9, ] <- c('Predicted (van de Pas)', VdP_table$No.[6], sum(results_VdP)/sum(Control_Chol_levels))
HDL.matrix[7, ] <- c('Experimental', VdP_table$No.[6], VdP_table$HDL.C..Relative.to.Control.[6])
LDL.matrix[7, ] <- c('Experimental', VdP_table$No.[6], VdP_table$non.HDL.C.Relative.to.Control.[6])
TC.matrix[7, ] <- c('Experimental', VdP_table$No.[6], VdP_table$TC.Relative.to.Control.[6])

#LCAT_het validation

fmut_LCAT_het <- fmut # fmut of selected rate 1, 2, 3, 5, 7, 8, 9, 16, 17, 18, 21
LCAT_het_mut <- mean(fmut.table$mean[grep("^LCAT$", as.character(fmut.table$GENE))])
fmut_LCAT_het[9] <- LCAT_het_mut
results_TL <- VdP_model(fmut_LCAT_het)
HDL.matrix[11, ] <- c('Predicted (genotype)', VdP_table$No.[7], results_TL[1]/Control_Chol_levels[1])
LDL.matrix[11, ] <- c('Predicted (genotype)', VdP_table$No.[7], results_TL[2]/Control_Chol_levels[2])
TC.matrix[11, ] <- c('Predicted (genotype)', VdP_table$No.[7], sum(results_TL)/sum(Control_Chol_levels))
fmut_LCAT_het[9] <- VdP_table$fmut[7]
results_VdP <- VdP_model(fmut_LCAT_het)
HDL.matrix[12, ] <- c('Predicted (van de Pas)', VdP_table$No.[7], results_VdP[1]/Control_Chol_levels[1])
LDL.matrix[12, ] <- c('Predicted (van de Pas)', VdP_table$No.[7], results_VdP[2]/Control_Chol_levels[2])
TC.matrix[12, ] <- c('Predicted (van de Pas)', VdP_table$No.[7], sum(results_VdP)/sum(Control_Chol_levels))
HDL.matrix[10, ] <- c('Experimental', VdP_table$No.[7], VdP_table$HDL.C..Relative.to.Control.[7])
LDL.matrix[10, ] <- c('Experimental', VdP_table$No.[7], VdP_table$non.HDL.C.Relative.to.Control.[7])
TC.matrix[10, ] <- c('Experimental', VdP_table$No.[7], VdP_table$TC.Relative.to.Control.[7])

#LCAT_hom validation

fmut_LCAT_hom <- fmut # fmut of selected rate 1, 2, 3, 5, 7, 8, 9, 16, 17, 18, 21
LCAT_hom_mut <- mean(fmut.table$mean[grep("LCAT_HOM", as.character(fmut.table$GENE))])
fmut_LCAT_hom[9] <- LCAT_hom_mut 
results_TL <- VdP_model(fmut_LCAT_hom)
HDL.matrix[14, ] <- c('Predicted (genotype)', VdP_table$No.[8], results_TL[1]/Control_Chol_levels[1])
LDL.matrix[14, ] <- c('Predicted (genotype)', VdP_table$No.[8], results_TL[2]/Control_Chol_levels[2])
TC.matrix[14, ] <- c('Predicted (genotype)', VdP_table$No.[8], sum(results_TL)/sum(Control_Chol_levels))
fmut_LCAT_hom[9] <- VdP_table$fmut[8]
results_VdP <- VdP_model(fmut_LCAT_hom)
HDL.matrix[15, ] <- c('Predicted (van de Pas)', VdP_table$No.[8], results_VdP[1]/Control_Chol_levels[1])
LDL.matrix[15, ] <- c('Predicted (van de Pas)', VdP_table$No.[8], results_VdP[2]/Control_Chol_levels[2])
TC.matrix[15, ] <- c('Predicted (van de Pas)', VdP_table$No.[8], sum(results_VdP)/sum(Control_Chol_levels))
HDL.matrix[13, ] <- c('Experimental', VdP_table$No.[8], VdP_table$HDL.C..Relative.to.Control.[8])
LDL.matrix[13, ] <- c('Experimental', VdP_table$No.[8], VdP_table$non.HDL.C.Relative.to.Control.[8])
TC.matrix[13, ] <- c('Experimental', VdP_table$No.[8], VdP_table$TC.Relative.to.Control.[8])

#DHCR7 validation
fmut_DHCR7 <- fmut # fmut of selected rate 1, 2, 3, 5, 7, 8, 9, 16, 17, 18, 21 
fmut_DHCR7[c(1,2,3)] <- 0 # as computed during the training phase, no mutation information in the original paper, is an autosomal recessive disorder
results_TL <- VdP_model(fmut_DHCR7)
HDL.matrix[17, ] <- c('Predicted (genotype)', VdP_table$No.[9], results_TL[1]/Control_Chol_levels[1])
LDL.matrix[17, ] <- c('Predicted (genotype)', VdP_table$No.[9], results_TL[2]/Control_Chol_levels[2])
TC.matrix[17, ] <- c('Predicted (genotype)', VdP_table$No.[9], sum(results_TL)/sum(Control_Chol_levels))
fmut_DHCR7[c(1,2,3)] <- VdP_table$fmut[9]
results_VdP <- VdP_model(fmut_DHCR7)
HDL.matrix[18, ] <- c('Predicted (van de Pas)', VdP_table$No.[9], results_VdP[1]/Control_Chol_levels[1])
LDL.matrix[18, ] <- c('Predicted (van de Pas)', VdP_table$No.[9], results_VdP[2]/Control_Chol_levels[2])
TC.matrix[18, ] <- c('Predicted (van de Pas)', VdP_table$No.[9], sum(results_VdP)/sum(Control_Chol_levels))
HDL.matrix[16, ] <- c('Experimental', VdP_table$No.[9], VdP_table$HDL.C..Relative.to.Control.[9])
LDL.matrix[16, ] <- c('Experimental', VdP_table$No.[9], VdP_table$non.HDL.C.Relative.to.Control.[9])
TC.matrix[16, ] <- c('Experimental', VdP_table$No.[9], VdP_table$TC.Relative.to.Control.[9])

#CYP7A1 validation

fmut_CYP7A1 <- fmut # fmut of selected rate 1, 2, 3, 5, 7, 8, 9, 16, 17, 18, 21 
CYP7A1_mut <- mean(fmut.table$mean[grep("CYP7A1", as.character(fmut.table$GENE))])
fmut_CYP7A1[18] <- CYP7A1_mut # no mutation information in the original paper
results_TL <- VdP_model(fmut_CYP7A1)
HDL.matrix[20, ] <- c('Predicted (genotype)', VdP_table$No.[10], results_TL[1]/Control_Chol_levels[1])
LDL.matrix[20, ] <- c('Predicted (genotype)', VdP_table$No.[10], results_TL[2]/Control_Chol_levels[2])
TC.matrix[20, ] <- c('Predicted (genotype)', VdP_table$No.[10], sum(results_TL)/sum(Control_Chol_levels))
fmut_CYP7A1[18] <- VdP_table$fmut[10]
results_VdP <- VdP_model(fmut_CYP7A1)
HDL.matrix[21, ] <- c('Predicted (van de Pas)', VdP_table$No.[10], results_VdP[1]/Control_Chol_levels[1])
LDL.matrix[21, ] <- c('Predicted (van de Pas)', VdP_table$No.[10], results_VdP[2]/Control_Chol_levels[2])
TC.matrix[21, ] <- c('Predicted (van de Pas)', VdP_table$No.[10], sum(results_VdP)/sum(Control_Chol_levels))
HDL.matrix[19, ] <- c('Experimental', VdP_table$No.[10], VdP_table$HDL.C..Relative.to.Control.[10])
LDL.matrix[19, ] <- c('Experimental', VdP_table$No.[10], VdP_table$non.HDL.C.Relative.to.Control.[10])
TC.matrix[19, ] <- c('Experimental', VdP_table$No.[10], VdP_table$TC.Relative.to.Control.[10])


final.TC.matrix <- rbind(TC.matrixI,TC.matrix[1:21, ])
final.HDL.matrix <- rbind(HDL.matrixI,HDL.matrix[1:21, ])
final.LDL.matrix <- rbind(LDL.matrixI, LDL.matrix[1:21, ])


# print results

write.table(file = './results/paper/HDL_total_chol.txt', final.HDL.matrix, sep = '\t', row.names = F, col.names = T)
write.table(file = './results/paper/LDL_total_chol.txt', final.LDL.matrix, sep = '\t', row.names = F, col.names = T)
write.table(file = './results/paper/total_chol.txt', final.TC.matrix, sep = '\t', row.names = F, col.names = T)

# print fmut table

f.matrix <- matrix(data = rep(NA, 10), ncol = 1)

f.matrix[1, 1] <- LDLR_effect
f.matrix[2, 1] <- ApoB_Het
f.matrix[3, 1] <- APOB_hom
f.matrix[4, 1] <- ABCA1_mut
f.matrix[5, 1] <- APOE_mut
f.matrix[6, 1] <- CETP_mut
f.matrix[7, 1] <- LCAT_het_mut
f.matrix[8, 1] <- LCAT_hom_mut
f.matrix[9, 1] <- 0 # as computed during the training phase, no mutation information in the original paper, is an autosomal recessive disorder
f.matrix[10, 1] <- CYP7A1_mut

f.matrix <- cbind(as.character(VdP_table$Gene), round(as.numeric(f.matrix), digits = 2), VdP_table$fmut)
colnames(f.matrix) <- c("Gene", "trained fmut", "vdP et al fmut")
write.table(file = './results/paper/fmut_paper.txt', f.matrix, sep = '\t', row.names = F, col.names = T, quote = F)


# make results plot 
png("./results/paper/Plot_total_chol.png", width = 7, height = 7, units = 'in', res = 300)
myData=read.table('./results/paper/total_chol.txt', header=TRUE)
limits <- aes(ymax = myData$mean + myData$se, ymin = myData$mean - myData$se)
p <- ggplot(data = myData, aes(x = factor(type), y = mean, fill = factor(Nomi)))
p + geom_bar(stat = "identity",
             position = position_dodge(0.9)) +
  labs(x = "Mutation", y = "Blood total cholesterol level case/control") +
  ggtitle("Total cholesterol predicted levels") +
  scale_fill_discrete(name = "Labels") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

png("./results/paper/HDL_chol_prediction.png", width = 7, height = 7, units = 'in', res = 300)
myData=read.table('./results/paper/HDL_total_chol.txt', header=TRUE)
limits <- aes(ymax = myData$mean + myData$se, ymin = myData$mean - myData$se)
p <- ggplot(data = myData, aes(x = factor(type), y = mean, fill = factor(Nomi)))
p + geom_bar(stat = "identity",
             position = position_dodge(0.9)) +
  labs(x = "Mutation", y = "Blood HDL level case/control") +
  ggtitle("HDL cholesterol predicted levels") +
  scale_fill_discrete(name = "Labels") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

png("./results/paper/LDL_chol_prediction.png", width = 7, height = 7, units = 'in', res = 300)
myData=read.table('./results/paper/LDL_total_chol.txt', header=TRUE)
limits <- aes(ymax = myData$mean + myData$se, ymin = myData$mean - myData$se)
p <- ggplot(data = myData, aes(x = factor(type), y = mean, fill = factor(Nomi)))
p + geom_bar(stat = "identity",
             position = position_dodge(0.9)) +
  labs(x = "Mutation", y = "Blood LDL level case/control") +
  ggtitle("LDL cholesterol predicted levels") +
  scale_fill_discrete(name = "Labels") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()



# table of error of predicted - real for each gene
vdP.HDL.predicted <- as.numeric(final.HDL.matrix[grep('Pas', final.HDL.matrix), 3])
vdP.LDL.predicted <- as.numeric(final.LDL.matrix[grep('Pas', final.LDL.matrix), 3])
vdP.TC.predicted <- as.numeric(final.TC.matrix[grep('Pas', final.TC.matrix), 3])
PM.HDL.predicted <- as.numeric(final.HDL.matrix[grep('genotype', final.HDL.matrix), 3])
PM.LDL.predicted <- as.numeric(final.LDL.matrix[grep('genotype', final.LDL.matrix), 3])
PM.TC.predicted <- as.numeric(final.TC.matrix[grep('genotype', final.TC.matrix), 3])
exp.HDL.predicted <- as.numeric(final.HDL.matrix[grep('Experimental', final.HDL.matrix), 3])
exp.LDL.predicted <- as.numeric(final.LDL.matrix[grep('Experimental', final.LDL.matrix), 3])
exp.TC.predicted <- as.numeric(final.TC.matrix[grep('Experimental', final.TC.matrix), 3])


error.matrix <- matrix(data = 0, nrow = length(vdP.HDL.predicted), ncol = 8)
error.matrix[, 1] <- seq(1:nrow(error.matrix))
error.matrix[, 3] <-  vdP.HDL.predicted-exp.HDL.predicted 
error.matrix[, 4] <-  vdP.LDL.predicted-exp.LDL.predicted 
error.matrix[, 5] <-  vdP.TC.predicted-exp.TC.predicted
error.matrix[, 6] <-  PM.HDL.predicted-exp.HDL.predicted 
error.matrix[, 7] <-  PM.LDL.predicted-exp.LDL.predicted
error.matrix[, 8] <-  PM.TC.predicted-exp.TC.predicted
error.matrix <- round(error.matrix, digits = 2)
colnames(error.matrix) <- c("Mutation", "Gene", "HDL vdP", "LDL vdP", "TC vdP", "HDL GM", "LDL GM", "TC GM")
error.matrix[, 2] <- c("LDLR", "APOB het", "APOB hom", "ABCA1", "APOE", "CETP", "LCAT het", "LCAT hom", "DHCR7", "CYP7A1")
write.table(file = './results/paper/deviation_predicted_experimental_cholesterol.txt', error.matrix, sep = '\t', row.names = F, col.names = T, quote = F)

#  compute percentage of error

error.matrix <- matrix(data = 0, nrow = length(vdP.HDL.predicted), ncol = 8)
error.matrix[, 3] <-  (vdP.HDL.predicted-exp.HDL.predicted)/ exp.HDL.predicted
error.matrix[, 4] <-  (vdP.LDL.predicted-exp.LDL.predicted) / exp.LDL.predicted
error.matrix[, 5] <-  (vdP.TC.predicted-exp.TC.predicted) / exp.TC.predicted
error.matrix[, 6] <-  (PM.HDL.predicted-exp.HDL.predicted) / exp.HDL.predicted
error.matrix[, 7] <-  (PM.LDL.predicted-exp.LDL.predicted) / exp.LDL.predicted
error.matrix[, 8] <-  (PM.TC.predicted-exp.TC.predicted) / exp.TC.predicted
error.matrix <- round(error.matrix*100, digits = 1)
error.matrix[, 1] <- seq(1:nrow(error.matrix))
colnames(error.matrix) <- c("Mutation", "Gene", "HDL vdP", "LDL vdP", "TC vdP", "HDL GM", "LDL GM", "TC GM")
error.matrix[, 2] <- c("LDLR", "APOB het", "APOB hom", "ABCA1", "APOE", "CETP", "LCAT het", "LCAT hom", "DHCR7", "CYP7A1" )
write.table(file = './results/paper/percentage_deviation_predicted_experimental_cholesterol.txt', error.matrix, sep = '\t', row.names = F, col.names = T, quote = F)

# compute matrix of results

error.matrix <- matrix(data = 0, nrow = length(vdP.HDL.predicted), ncol = 11)
error.matrix[, 3] <-  exp.HDL.predicted
error.matrix[, 4] <-  exp.LDL.predicted
error.matrix[, 5] <-  exp.TC.predicted
error.matrix[, 6] <-  vdP.HDL.predicted
error.matrix[, 7] <-  vdP.LDL.predicted
error.matrix[, 8] <-  vdP.TC.predicted
error.matrix[, 9] <-  PM.HDL.predicted
error.matrix[, 10] <-  PM.LDL.predicted
error.matrix[, 11] <-  PM.TC.predicted
error.matrix <- round(error.matrix, digits = 2)
error.matrix[, 1] <- seq(1:nrow(error.matrix))
colnames(error.matrix) <- c("Mutation", "Gene", "HDL real", "LDL real", "TC real", "HDL vdP", "LDL vdP", "TC vdP", "HDL GM", "LDL GM", "TC GM")
error.matrix[, 2] <- c("LDLR", "APOB het", "APOB hom", "ABCA1", "APOE", "CETP", "LCAT het", "LCAT hom", "DHCR7", "CYP7A1" )
write.table(file = './results/paper/predicted_values.txt', error.matrix, sep = '\t', row.names = F, col.names = T, quote = F)







# final results

#remove NA rows 
id.rem <- which(is.na(final.LDL.matrix[grep('Experimental', final.LDL.matrix), 3])==FALSE) # in some elements of the test set only TC values are available


vdP.HDL.predicted <- as.numeric(final.HDL.matrix[grep('Pas', final.HDL.matrix), 3])[id.rem]
vdP.LDL.predicted <- as.numeric(final.LDL.matrix[grep('Pas', final.LDL.matrix), 3])[id.rem]
vdP.TC.predicted <- as.numeric(final.TC.matrix[grep('Pas', final.TC.matrix), 3])
PM.HDL.predicted <- as.numeric(final.HDL.matrix[grep('genotype', final.HDL.matrix), 3])[id.rem]
PM.LDL.predicted <- as.numeric(final.LDL.matrix[grep('genotype', final.LDL.matrix), 3])[id.rem]
PM.TC.predicted <- as.numeric(final.TC.matrix[grep('genotype', final.TC.matrix), 3])
exp.HDL.predicted <- as.numeric(final.HDL.matrix[grep('Experimental', final.HDL.matrix), 3])[id.rem]
exp.LDL.predicted <- as.numeric(final.LDL.matrix[grep('Experimental', final.LDL.matrix), 3])[id.rem]
exp.TC.predicted <- as.numeric(final.TC.matrix[grep('Experimental', final.TC.matrix), 3])

HDL.LDL.data <- as.data.frame(cbind(exp.HDL.predicted, exp.LDL.predicted, vdP.HDL.predicted,vdP.LDL.predicted,PM.HDL.predicted,PM.LDL.predicted))
TC.data <- as.data.frame(cbind(exp.TC.predicted, vdP.TC.predicted, PM.TC.predicted))

final.score.matrix <- matrix(0, nrow = 6, ncol = 5)
colnames(final.score.matrix) <- c('PCC', 'KCC', 'MAE', 'RMSD', 'Rsquared')
rownames(final.score.matrix) <- c('vdP predicted HDL ratio', 'vdP predicted LDL ratio', 'vdp predicted TC ratio', 'PM predicted HDL ratio', 'PM predicted LDL ratio', 'PM predicted TC ratio')
final.score.matrix[1, 1] <-  cor(vdP.HDL.predicted, exp.HDL.predicted, method ='pearson' )
final.score.matrix[1, 2] <-  cor(vdP.HDL.predicted, exp.HDL.predicted, method ='kendall', use = 'pairwise')
final.score.matrix[1, 3] <-  sum(abs(vdP.HDL.predicted-exp.HDL.predicted))/length(exp.HDL.predicted) 
final.score.matrix[1, 4] <-  sqrt(sum((vdP.HDL.predicted-exp.HDL.predicted)^2)/length(exp.HDL.predicted)) 
final.score.matrix[1, 5] <- summary(lm(vdP.HDL.predicted~exp.HDL.predicted, data = HDL.LDL.data ))$r.squared
final.score.matrix[2, 1] <-  cor(vdP.LDL.predicted, exp.LDL.predicted, method ='pearson' )
final.score.matrix[2, 2] <-  cor(vdP.LDL.predicted, exp.LDL.predicted, method ='kendall', use = 'pairwise')
final.score.matrix[2, 3] <-  sum(abs(vdP.LDL.predicted-exp.LDL.predicted))/length(exp.LDL.predicted) 
final.score.matrix[2, 4] <-  sqrt(sum((vdP.LDL.predicted-exp.LDL.predicted)^2)/length(exp.LDL.predicted)) 
final.score.matrix[2, 5] <- summary(lm(vdP.LDL.predicted~exp.LDL.predicted, data = HDL.LDL.data ))$r.squared
final.score.matrix[3, 1] <-  cor(vdP.TC.predicted, exp.TC.predicted, method ='pearson' )
final.score.matrix[3, 2] <-  cor(vdP.TC.predicted, exp.TC.predicted, method ='kendall', use = 'pairwise')
final.score.matrix[3, 3] <-  sum(abs(vdP.TC.predicted-exp.TC.predicted))/length(exp.TC.predicted)
final.score.matrix[3, 4] <-  sqrt(sum((vdP.TC.predicted-exp.TC.predicted)^2)/length(exp.TC.predicted))
final.score.matrix[3, 5] <- summary(lm(vdP.TC.predicted~exp.TC.predicted, data = TC.data ))$r.squared
final.score.matrix[4, 1] <-  cor(PM.HDL.predicted, exp.HDL.predicted, method ='pearson' )
final.score.matrix[4, 2] <-  cor(PM.HDL.predicted, exp.HDL.predicted, method ='kendall', use = 'pairwise')
final.score.matrix[4, 3] <-  sum(abs(PM.HDL.predicted-exp.HDL.predicted))/length(exp.HDL.predicted) 
final.score.matrix[4, 4] <-  sqrt(sum((PM.HDL.predicted-exp.HDL.predicted)^2)/length(exp.HDL.predicted)) 
final.score.matrix[4, 5] <- summary(lm(PM.HDL.predicted~exp.HDL.predicted, data = HDL.LDL.data ))$r.squared
final.score.matrix[5, 1] <-  cor(PM.LDL.predicted, exp.LDL.predicted, method ='pearson' )
final.score.matrix[5, 2] <-  cor(PM.LDL.predicted, exp.LDL.predicted, method ='kendall', use = 'pairwise')
final.score.matrix[5, 3] <-  sum(abs(PM.LDL.predicted-exp.LDL.predicted))/length(exp.LDL.predicted)
final.score.matrix[5, 4] <-  sqrt(sum((PM.LDL.predicted-exp.LDL.predicted)^2)/length(exp.LDL.predicted))
final.score.matrix[5, 5] <- summary(lm(PM.LDL.predicted~exp.LDL.predicted, data = HDL.LDL.data ))$r.squared
final.score.matrix[6, 1] <-  cor(PM.TC.predicted, exp.TC.predicted, method ='pearson' )
final.score.matrix[6, 2] <-  cor(PM.TC.predicted, exp.TC.predicted, method ='kendall', use = 'pairwise')
final.score.matrix[6, 3] <-  sum(abs(PM.TC.predicted-exp.TC.predicted))/length(exp.TC.predicted) 
final.score.matrix[6, 4] <-  sqrt(sum((PM.TC.predicted-exp.TC.predicted)^2)/length(exp.TC.predicted)) 
final.score.matrix[6, 5] <- summary(lm(PM.TC.predicted~exp.TC.predicted, data = TC.data ))$r.squared

write.table(file = './results/paper/score.table_LM_completo.txt', round(final.score.matrix, digits = 2), sep = '\t', row.names = T, col.names = T)

if(T){
  # make Boxplot dataset
  makePlotDataset(patient.dataset.name = "data/Other_genes_database.txt")
}
if(T){
  # make sensitivity plot
  MakeSimulations()
} 
if(T){
  # make bootstrap
  Bootstrap()
}
if(T){
  # compute standard deviation of model predictions
  SDpredictions()
}





