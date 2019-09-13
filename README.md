# Scripts ro run in silico model for cholesterol level prediction 
 

## REQUIREMENTS
```
R:
>install.packages('deSolve')
>install.packages('pracma')
>install.packages('minpack.lm')
>install.packages('plotrix')
>install.packages('gridExtra)
>install.packages('ggplot2')
```


## USAGE

to run all analyses, you have just to run the following command inside the src folder: Rscript ./main_VdP_V01_revision.R
      This script will compute all statistics and make plots
      This script expects a folder structure like this to run:
            * ./src : contains all scripts
            * ./results : will contain all performance tables and plots
            * ./results/paper/ : will contain all paper tables and plots    
            * ./data : will contain all input tables 
      
      These are Input needed:
            * input tables in ./data


      Running the script these Output files will be generated in ./results
           * table with fmut optimized (HDL, LDL)
           * table with fmut optimized (TC), only for DHCR7 gene
           * boxplot of cholesterol levels (training set)
           * table with difference between predicted and experimental cholesterol
           * table with trained vs van de Pas et al. fmut
           * sensitivity analysis plot (HDL, LDL, TC)
           * barplot with HDL levels (predicted, real)
           * table with HDL levels (predicted, real)
           * barplot with LDL levels (predicted, real)
           * table with LDL levels (predicted, real)
           * table with number of mutations and patients for each gene
           * percentage of difference between predicted and experimental
           * barplot with TC levels (predicted, real)
           * table with predicted and real blood cholesterol values
           * table with PCC, KCC, RMSD and Rsquared of predicted blood cholesterol values
           * table with TC levels (predicted, real)
           * table with standard deviation of model predictions
           * table with bootstrap




