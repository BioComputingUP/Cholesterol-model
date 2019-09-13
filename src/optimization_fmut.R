# find the value of fmut that minimize the difference between the value of cholesterol predicted by the model and the patient level. 

library(deSolve)
library(pracma)
library('minpack.lm')


l<-function(f_fact){
  fmut <- rep(1,21)
  fmut[reaction_model] <- f_fact
  x0 = c(1.99, 0.25, 8.00, 5.30, 0.30, 0.89, 4.03, 1.20)
  times <- seq(0, 10000, by=1 )
  
  V_liver = 1.8
  V_plasma = 2.79
  V_int = 0.64
  V_per = 64.8
  V_tot = V_liver + V_plasma + V_int + V_per
  bw = 70
  statin_eff=1
  
  # [C_6 C_8 C_1 C_7 C_2 C_3 C_4 C_5]
  C_6 = x0[1]
  C_8 = x0[2]
  C_1 = x0[3]
  C_7 = x0[4]
  C_2 = x0[5]
  C_3 = x0[6]
  C_4 = x0[7]
  C_5 = x0[8]
  E_LDLR = exp(-3.92*(C_1+C_7)/13.3)/exp(-3.92)
  
  
  
  Y=c(v5 = fmut[5]*1.0157*C_4*V_plasma, v6 = 0.4990*C_7*V_liver, v7 = fmut[7]*0.1165*C_4*V_plasma, v8 = fmut[8]*0.0319*C_5*V_per, v9 = fmut[9]*9.3907*C_2*V_plasma, v10 = 1.18*C_3*V_plasma, v11 = 19.0*C_8*V_int, v12 = 0.0337*C_5*V_per, v13 = 1.8757*C_2*V_plasma, v14 = (0.2542*C_1*0.5 + 2.0333*0.5)*V_liver, v15 = 1.4526*C_6*V_int, v16 = fmut[16]*(0.0314*C_6*0.75 + 0.0625*0.25)*V_int, v19 = (0.3306*C_1*0.25 + 2.6444*0.75)*V_liver, v20 = (2.3869*C_6*0.75 + 4.75*0.25)*V_int, v21 = fmut[21]*0.4927*C_3*C_4*V_plasma, C_1 = x0[3], C_2 = x0[5], C_3 = x0[6], C_4 = x0[7], C_5 = x0[8], C_6 = x0[1],   C_7 = x0[4], C_8 = x0[2])
  pars=c(V_liver = 1.8, V_plasma = 2.79, V_int = 0.64, V_per = 64.77, V_tot = V_liver + V_plasma + V_int + V_per, bw = 70, v1 = fmut[1]*(0.0063*V_tot)*statin_eff, v2 = fmut[2]*(0.0541*V_tot), v3 = fmut[3]*(0.0026*V_tot), v4 = 0.0156*V_tot, v17 = fmut[17]*3.8389*V_liver, v18 = fmut[18]*0.5722*V_liver)
  mouse_kinetic_model <- function(times, Y, pars) {# TimeStep #InitialValues #KonstantParameters
    
    with(as.list(c(Y, pars)), {
      
      
      
      # C
      
      
      v5 = fmut[5]*1.0157*C_4*V_plasma;
      v6 = 0.4990*C_7*V_liver;
      v7 = fmut[7]*0.1165*C_4*V_plasma;
      v8 = fmut[8]*0.0319*C_5*V_per;
      v9 = fmut[9]*9.3907*C_2*V_plasma;
      v10 = 1.18*C_3*V_plasma;
      v11 = 19.0*C_8*V_int;
      v12 = 0.0337*C_5*V_per;
      v13 = 1.8757*C_2*V_plasma;
      v14 = (0.2542*C_1*0.5 + 2.0333*0.5)*V_liver;
      v15 = 1.4526*C_6*V_int;
      v16 = fmut[16]*(0.0314*C_6*0.75 + 0.0625*0.25)*V_int;
      v19 = (0.3306*C_1*0.25 + 2.6444*0.75)*V_liver;
      v20 = (2.3869*C_6*0.75 + 4.75*0.25)*V_int;
      v21 = fmut[21]*0.4927*C_3*C_4*V_plasma;
      
      C_1 <- (v1+v13+v10+v5-v18-v19-v14-v17)/V_liver
      C_2 <- (v16+v17+v8-v9-v13)/V_plasma
      C_3 <- (v9-v10-v21)/V_plasma
      C_4 <- (v11+v6-v5-v7+v21)/V_plasma
      C_5 <- (v7+v2-v12-v8)/V_per
      C_6 <- (v4+v14+v3-v16-v20-v15)/V_int
      C_7 <- (v19-v6)/V_liver
      C_8 <- (v20-v11)/V_int
      
      
      
      
      
      
      
      return(list(c(v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,v19,v20,v21,C_1,C_2,C_3,C_4,C_5,C_6,C_7,C_8)))
    })
  }
  out <- ode(y = Y, 
             times = times, 
             func = mouse_kinetic_model, 
             parms = pars)
  res <- out[,17:24]
  vd=c(1:8)
  predicted_HDL_FC=out[nrow(out),18]
  predicted_HDL_CE=out[nrow(out),19]
  predicted_LDL=out[nrow(out),20]
  return(c(predicted_HDL_FC+predicted_HDL_CE, predicted_LDL))
} # par function(a)


Compute_f_mut_TC <- function(){ # used to optimize fmut for DHCR7 patients
  T.S <- read.table(file = "data/DHCR7_dataset.txt", header=TRUE, sep = '\t') # training set
  f_mut <- rep(1, nrow(T.S))
  for(i in 1:nrow(T.S)){
    ID=T.S[i,1]
    reaction_model <<-  as.numeric(unlist(strsplit(as.character(T.S[i, 4]), ","))) # set model reactions multiplied by fmut, set as global variable to work inside function
    if(isTRUE(T.S$mmol.l[i])){
      TC <- T.S[i, 3]
    } else {
      TC <-T.S[i, 3]/38.67 # convert from mg/dL to mmol/l
    }
    real.output <- TC
    n<-function(f_fact){y<-(real.output-sum(l(f_fact)))^2; return(y)} # objective function (TC_real_value-TC_predicted_value)
    v <- nls.lm(par = 0.5 ,fn = n , upper = 1 , lower = 0 , control = nls.control( maxiter = 200 ))
    #v<-fminbnd(n,0,1) # compute the value of the paramether a (fmut) that minimizes the objective function 
    f_mut[i] <- v$xmin
  }
  return(f_mut)
}


Compute_f_mut_LDL <- function(){
  T.S <- read.table(file = "data/Other_genes_database.txt", header=TRUE, sep = '\t') # training set
  f_mut <- rep(1, nrow(T.S))
  for(i in 1:nrow(T.S)){
    ID=T.S[i,1]
    reaction_model <<-  as.numeric(unlist(strsplit(as.character(T.S[i,5]), ","))) # set model reactions multiplied by fmut, set as global variable to work inside function
    if(isTRUE(T.S$mmol.l[i])){
      HDL=T.S[i,3]/1.19
      LDL=T.S[i,4]/4.03
    } else {
      HDL=(T.S[i,3]/38.67)/1.19 # convert from mg/dL to mmol/l
      LDL=(T.S[i,4]/38.67)/4.03
    }
    real.output <- c(HDL,LDL)
    n<-function(f_fact){y<-abs(real.output-(l(f_fact)/c(1.19, 4.03))); return(y)} # objective function
    v <- nls.lm(par = 0.5 ,fn = n , upper = 1 , lower = 0 , control = nls.control( maxiter = 200 ))
    f_mut[i] <- v$par
    
  }
  return(f_mut)
}
