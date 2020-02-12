#!/usr/bin/env/ Rscript

makePlotDataset <- function(patient.dataset.name = "data/Other_genes_database.txt"){
  title.document <- read.table(file = patient.dataset.name, header=TRUE, sep = '\t', stringsAsFactors = F) # training set
  HDL <- title.document$C.HDL
  HDL[which(title.document$mmol.l!=T)] <- HDL[which(title.document$mmol.l!=T)]/38.67 # make all in mmol/L 
  LDL <- title.document$C.LDL
  LDL[which(title.document$mmol.l!=T)] <- LDL[which(title.document$mmol.l!=T)]/38.67 
  class.genes <- unique(title.document$GENE)
  class.column <- title.document$GENE
  
  for(i in (1:length(class.genes))){
    class.column[which(class.column == class.genes[i])] <- i
  }
  # change wrong class ID
  w.id <- which(class.column == "7") 
  class.column[which(class.column == "8")] <- "7"
  class.column[w.id] <- "8"
  class.column[which(class.column == "9")] <- "10"
  
  # create data.frame for plot 
  HDL <- c(HDL, 1.19)
  LDL <- c(LDL, 4.03)
  class.column <- c(class.column, "0")
  HDL <- HDL[order(as.numeric(class.column))] 
  LDL <- LDL[order(as.numeric(class.column))] 
  class.column <- class.column[order(as.numeric(class.column))]
  class.columnI <- class.column # make column to split boxplot
  class.columnI[which(as.numeric(class.column) < 4)] <- "Autosomal Dominant Hypercholesterolemia" 
  class.columnI[which(as.numeric(class.column) == 5)] <- "Autosomal Dominant Hypercholesterolemia" 
  class.columnI[which(class.columnI != "Autosomal Dominant Hypercholesterolemia")] <- "Other diseases affecting lipoprotein metabolism"
  class.columnI[1] <- "Control"
  class.column2 <- rep("Case", length(class.column))
  class.column2[1] <- c("Control HDL")
  class.column3 <- class.column2 
  class.column3[1] <- c("Control LDL")
  class.column[which(class.column == "0")] <- "Control"
  Gene.chol <- data.frame(rbind(cbind(HDL, rep("HDL-C",length(HDL)), class.column, class.columnI, class.column2),cbind(LDL, rep("LDL-C",length(LDL)), class.column, class.columnI, class.column3)) )
  colnames(Gene.chol) <- c("amount", "cholesterol", "gene", "dataset", "color_line")
  Gene.chol$cholesterol <- as.factor(Gene.chol$cholesterol)
  Gene.chol$gene <- factor(Gene.chol$gene, levels = c("Control", 1, 2, 3, 5, 4, 6, 7, 8, 10))
  Gene.chol$amount <- as.numeric(c(HDL,LDL))
  Gene.chol$dataset <- factor(Gene.chol$dataset, levels = c("Control", "Autosomal Dominant Hypercholesterolemia", "Other diseases affecting lipoprotein metabolism"))
  Gene.chol$color_line <- factor(Gene.chol$color_line)
  #png("./results/Boxplot_chol_levels.png", width = 7, height = 7, units = 'in', res = 300)
  tiff("results/paper/Boxplot_chol_levels.tiff", units="in", width=7, height=7, res=300)
  
  p <- ggplot(Gene.chol, aes(x = gene, y = amount, colour = color_line, fill = cholesterol)) +
    ggtitle("Training set cholesterol levels") +
    geom_boxplot()+ facet_grid(. ~ dataset, scales = "free", space = "free") + scale_color_manual(values=c("black", "indianred2", "turquoise3"), guide=F) + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2), legend.title = element_blank() ,  panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                      strip.text.x = element_text(size = 7),
                                      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
                                      axis.title.y = element_text(size = 14, face = "bold"),
                                      axis.text.x = element_text(angle = 300, hjust = 0.2, vjust = 0.8, face = "bold"),
                                      axis.title.x = element_text(size = 14, face = "bold"),
                                      axis.text = element_text(size = 8))
  p <- p + scale_x_discrete(breaks=c("Control", "1", "2", "3", "5", "4", "6", "7", "8", "10"),
                       labels = c("Control", "LDLR", "APOB", "APOB hom", "APOE", "ABCA1", "CETP", "LCAT", "LCAT hom", "CYP7A1")) + theme(axis.text = element_text(size = 12)) + labs(x = "Gene affected by the mutation", y = "Cholesterol levels [mmol/L]")
  print(p)
  dev.off()
  if(FALSE){
    png("./results/Boxplot_scatterplot_chol_levels.png", width = 7, height = 7, units = 'in', res = 300)
    p <- ggplot(Gene.chol, aes(x = gene, y = amount, fill = cholesterol)) +
      ggtitle("Training set cholesterol levels") +   geom_point(aes(color= cholesterol), position=position_dodge(width=0.5)) +
      theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                         plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
                         #axis.text.y = element_text(size = 14),
                         axis.title.y = element_text(size = 14, face = "bold"),
                         #axis.text.x = element_text(size = 14),
                         axis.title.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                         element_text(size = 14, face = "bold"),
                         axis.text = element_text(size = 8))
    p +  scale_x_discrete(limits=c("Control", "1", "2", "3", "4", "5", "6", "7", "8", "10")) + theme(axis.text = element_text(size = 12)) + labs(x = "Gene affected by the mutation", y = "mmol/L")
    dev.off()
    }
  }
