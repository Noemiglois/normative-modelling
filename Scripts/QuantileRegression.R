####### Prepare data ##########
library("quantreg") # install.packages("quantreg")
setwd("/data_J/Scripts")
source("1_DataPreparation.R")
df_raw     <- read.csv("/data_J/Data/morphosim_lC_harmonized_data.csv", header = T)
df_NOmatch <- DataPreparation(df = df_raw, harmonization = F, match = F)
df <- df_NOmatch
reg_names  <- colnames(df%>%select(ends_with("thickness"))) # Select thickness features names
df <- df[,c(reg_names,"age","scode","subID","ID","euler","timepoint","dcode","dcode_age")]

df_tp1_pat <- df %>%
  filter(timepoint == 1)%>%
  filter(dcode == 1)

df_tp1_hc <- df %>%
  filter(timepoint == 1)%>%
  filter(dcode == 0)

# cov_hc <- df_tp1_hc[ , c("age","scode","subID","ID","euler","timepoint","dcode","dcode_age")]
# x_hc <- cbind(cov_hc$age + cov_hc$scode + cov_hc$euler) # x has to be a matrix, not a df
# y_hc   <- df_tp1_hc[ , reg_names]
attach(df_tp1_hc)
X <- cbind(age + scode + euler)
rm(df_raw, df_NOmatch, df, df_tp1_pat, DataPreparation)

############ Regression ###########
taus <- c(0.05, 0.5,0.9)
for (t in taus[1:1]){
  #cat('tau:',t)
  for(i in 1:length(reg_names[1:1])){   
    #cat('reg',i)
    Y <- cbind(formula(noquote(reg_names[i])))
    QR_orig = rq(Y ~ X , data=df_tp1_hc, tau=t, model=T) 
    summary(QR_orig)            # CI
    summary(QR_orig, se="boot") # SE
  }
}