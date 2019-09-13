MakeSimulations<- function(){
  fmut <- rep(1, 21) # fmut of selected rate 1, 2, 3, 5, 7, 8, 9, 16, 17, 18, 21 
  Control_Chol_levels <- VdP_model(fmut)
  #fix number of times you want to repeat simulation (number fmut)
  fmut.factor <- rev(seq(0.1:1, by = 0.1))
  genes.list <- list(c(5,7), c(8,16,17), c(21), c(9), c(1,2,3), c(18))
  #create matrix to store HDL,LDL,TC change in time
  HDL.levels <- matrix(data = 0, nrow = length(genes.list), ncol = length(fmut.factor))
  LDL.levels <- matrix(data = 0, nrow = length(genes.list), ncol = length(fmut.factor))
  TC.levels <- matrix(data = 0, nrow = length(genes.list), ncol = length(fmut.factor))
  for(j in (1:length(genes.list))){
    fmut <- rep(1, 21) # restore fmut from previous for cycle
    for(i in (1:length(fmut.factor))){
      fmut[genes.list[[j]]] <- fmut.factor[i]
      chol.pred <- VdP_model(fmut)
      HDL.levels[j, i] <- chol.pred[1]
      LDL.levels[j, i] <- chol.pred[2]
      TC.levels[j, i] <- sum(chol.pred)
    }
  }
  # compute cholesterol rlevels relative to control
  HDL.levels.control <- HDL.levels/Control_Chol_levels[1]
  LDL.levels.control <- LDL.levels/Control_Chol_levels[2]
  TC.levels.control <- TC.levels/sum(Control_Chol_levels)
  
  # color genes
  color.chol <- c("red", "blue", "yellow", "green", "black", "purple")
  # legend genes
  legend.genes <-  c("LDLR, APOB, APOE (5,7)", "ABCA1 (8,16,17)",
                     "CETP (21)", "LCAT (9)", "DHCR7 (1,2,3)", "CYP7A1 (18)")
  
  # plot HDL
  
  if(F){ # previous version plot
    tiff("results/paper/HDL.levels.mutation.tiff", units="in", width=15, height=10, res=300)
    par(mfrow = c(2,2), xpd = NA, mar = par()$mar + c(0,0,0,12) )#, pty= "s") # increase size of right margin for legend, pty = s force plot to be squared
    plot(-c(min(fmut.factor),max(fmut.factor)), c(min(HDL.levels.control), max(HDL.levels.control)), pch = "",
         main = "Effect of reduced model rates on HDL cholesterol ",
         xlab = "fmut",
         ylab = "HDL cholesterol case / control ", xaxt = "n" )
    axis(1, at = rev(-seq(0.1, 1, by = 0.1)), labels = rev(seq(0.1, 1, by = 0.1)))
    for(i in (1:nrow(HDL.levels.control))){
      lines(-fmut.factor, HDL.levels.control[i,], col = color.chol[i])
    }
    coord <- par("usr") # get plot coordinates
    legend(x = coord[2] * 0.95, y = coord[4], inset=c(-0.20,0.1), legend = legend.genes,
           col = color.chol, lty = 1, cex=1, title = "Gene (model rate)")
    plot(-c(min(fmut.factor),max(fmut.factor)), c(min(LDL.levels.control), max(LDL.levels.control)), pch = "",
         main = "Effect of reduced model rates on LDL cholesterol ",
         xlab = "fmut",
         ylab = "LDL cholesterol case / control ", xaxt = "n" )
    axis(1, at = rev(-seq(0.1, 1, by = 0.1)), labels = rev(seq(0.1, 1, by = 0.1)))
    for(i in (1:nrow(LDL.levels.control))){
      lines(-fmut.factor, LDL.levels.control[i,], col = color.chol[i])
    }
    coord <- par("usr") # get plot coordinates
    legend(x = coord[2] * 0.95, y = coord[4], inset=c(-0.20,0.1), legend = legend.genes,
           col = color.chol, lty = 1, cex=1, title = "Gene (model rate)")
    plot(-c(min(fmut.factor),max(fmut.factor)), c(min(TC.levels.control), 
                                                  max(TC.levels.control)), 
         pch = "", main = "Effect of reduced model rates on total cholesterol ",
         xlab = "fmut",
         ylab = "Total cholesterol case / control ", xaxt = "n" )
    axis(1, at = rev(-seq(0.1, 1, by = 0.1)), labels = rev(seq(0.1, 1, by = 0.1)))
    color.chol <- c("red", "blue", "yellow", "green", "black", "purple")
    for(i in (1:nrow(TC.levels.control))){
      lines(-fmut.factor, TC.levels.control[i,], col = color.chol[i])
    }
    coord <- par("usr") # get plot coordinates
    legend(x = coord[2] * 0.95, y = coord[4], inset=c(-0.20,0.1), legend = legend.genes,
           col = color.chol, lty = 1, cex=1, title = "Gene (model rate)")
    dev.off()
  } else { # last version plot
    legend.genes <-  c("LDLR, APOB, APOE (5,7)                          ", "ABCA1 (8,16,17)",
                       "CETP (21)", "LCAT (9)", "DHCR7 (1,2,3)", "CYP7A1 (18)")
    tiff("results/paper/HDL.levels.mutation.tiff", units="in", width=15, height=10, res=300)
    color.chol <- c("red", "blue", "orange", "green", "black", "purple")
    par(mfrow = c(1,3), cex = 0.65, xpd = NA, mar = par()$mar + c(0,0,0,0), pty= "s") # increase size of right margin for legend, pty = s force plot to be squared
    plot(-c(min(fmut.factor),max(fmut.factor)), c(min(HDL.levels.control), max(HDL.levels.control)), pch = "",
         main = "Effect of reduced model rates on HDL cholesterol ",
         xlab = "fmut",
         ylab = "HDL cholesterol case / control ", xaxt = "n" )
    axis(1, at = rev(-seq(0.1, 1, by = 0.1)), labels = rev(seq(0.1, 1, by = 0.1)))
    for(i in (1:nrow(HDL.levels.control))){
      lines(-fmut.factor, HDL.levels.control[i,], col = color.chol[i])
    }
    plot(-c(min(fmut.factor),max(fmut.factor)), c(min(LDL.levels.control), max(LDL.levels.control)), pch = "",
         main = "Effect of reduced model rates on LDL cholesterol ",
         xlab = "fmut",
         ylab = "LDL cholesterol case / control ", xaxt = "n" )
    axis(1, at = rev(-seq(0.1, 1, by = 0.1)), labels = rev(seq(0.1, 1, by = 0.1)))
    for(i in (1:nrow(LDL.levels.control))){
      lines(-fmut.factor, LDL.levels.control[i,], col = color.chol[i])
    }
    legend("bottom", legend = legend.genes, inset=c(-10,-0.3),
           col = color.chol, lwd = 3, lty = 1, cex=1, title = expression(bold("Gene (model rate)")), horiz = T,
           xjust = 0.5, yjust = 0.5, x.intersp = 0.5, text.font = 2)
    plot(-c(min(fmut.factor),max(fmut.factor)), c(min(TC.levels.control), 
                                                  max(TC.levels.control)), 
         pch = "", main = "Effect of reduced model rates on total cholesterol ",
         xlab = "fmut",
         ylab = "Total cholesterol case / control ", xaxt = "n" )
    axis(1, at = rev(-seq(0.1, 1, by = 0.1)), labels = rev(seq(0.1, 1, by = 0.1)))
    for(i in (1:nrow(TC.levels.control))){
      lines(-fmut.factor, TC.levels.control[i,], col = color.chol[i])
    }
    dev.off()
  }
}




