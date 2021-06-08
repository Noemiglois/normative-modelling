################################ IMPORTS ######################################
library(variancePartition) #install via bioconductor: https://bioconductor.org/packages/release/bioc/html/variancePartition.html
library(doParallel)
library(Rmisc)
library(dplyr)

################################ PREPARE DATA #################################
# Load dataset
dfvarpart <- read.csv("/data_J/data/lC_harmonized_thickness_reduced.csv")
reg_names <- colnames(dfvarpart%>%select(ends_with("thickness")))

# Rename predictors
names(dfvarpart)[names(dfvarpart) == "age"] <- "Age"
names(dfvarpart)[names(dfvarpart) == "scode"] <- "Sex"
names(dfvarpart)[names(dfvarpart) == "euler"] <- "Euler_number"
names(dfvarpart)[names(dfvarpart) == "dcode"] <- "Diagnosis"
names(dfvarpart)[names(dfvarpart) == "subID"] <- "Individual"

# Type coercion
dfvarpart$Sex <- as.factor(dfvarpart$Sex)
dfvarpart$Individual <- as.factor(dfvarpart$Individual)
dfvarpart$Diagnosis <- as.factor(dfvarpart$Diagnosis)
 
############################### VIOLIN PLOTS ##################################
n_tp<-c(1,2,3)
for (tp in n_tp){
  # Timepoint selection
  dfvarpart_tp <- dfvarpart %>%
    filter(timepoint==tp)
  
  # Select regions, transpose
  dfvarpart_base_ima <- as.data.frame(t(dfvarpart_tp[,reg_names]))
  
  # Create vectors for each predictor
  Age <- dfvarpart_tp$Age
  Sex <- dfvarpart_tp$Sex
  Diagnosis <- dfvarpart_tp$Diagnosis
  Euler_numer <- dfvarpart_tp$Euler_number
  Individual<-dfvarpart$Individual
  
  # Create formula
  # form <- ~Age+(1|Sex)+(1|Diagnosis)+Euler_numer+(1|Individual)
  form <- ~Age+Sex+Diagnosis+Euler_numer
  
  # Create dataframe with predictors
  info <- dfvarpart_tp[,c("Age", "Sex", "Diagnosis","Euler_number")]
  
  # Run variance partition
  varPart <- fitExtractVarPartModel(dfvarpart_base_ima, form, info)
  
  # Sort columns
  vp <- sortCols(varPart)

  # Plot sorted columns as violins, setting colors
  theme_set(theme_gray(base_size = 40))
  violins <- plotVarPart(vp,col = c("#5DC863FF", "#3B528BFF", "#21908CFF", "#440154FF", "lightgrey"))

  # Save as png
  png(file=paste0("/data_J/results/Plots/violin_tp",tp,"_",lab,".png"), width=5, height=2.5, units="in",res=600)
  print(violins)
  dev.off()
  
  # # Save as pdf
  # pdf(file = paste0("/data_J/results/Plots/violin_tp",tp,".pdf"), width = 12, height = 6)
  # multiplot(violins,cols = 1)
  # dev.off()
}




