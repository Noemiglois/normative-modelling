############################## IMPORTS #######################################
library(parallel)
library(lme4)
library(nlme)
library(pbmcapply)
library(JMbayes)
############################# DATASET  ########################################
df <- read.csv("/data_J/Data/morphosim_lC_harmonized_data.csv", header = T)
df <- transform(df, euler=(euler_lh+euler_rh)/2) # Create euler variable
df <- transform(df, dcode_age=dcode*age)         # Create dcode_age variable
df <- df[order(df$age), ]                        # Sort by age (not sure why)

reg_names <- colnames(df%>%select(ends_with("thickness"))) # Select thickness features names

idx_pat <- which(df$dcode== 1)  # 1: patient
idx_hc  <- which(df$dcode== 0)  # 0: healthy control
data_pat <- df[idx_pat, ]
data_hc  <- df[idx_hc, ]

pat_subIDs <- unique(data_pat$subID) # 169 patients
# hc_subIDs  <- unique(data_hc$subID)  # 298 healthy controls
# Total of 467 subjects
# hc_pat_subIDs <- unique(df$subID)

ncores = detectCores() # Number of cores to parallelize with mcapply

# vars<-c("age", "scode", "subID","ID", "euler", "timepoint","dcode", "dcode_age",reg_names)
# data_joost <- df[,vars]
# write.csv(data_joost,"/data_J/data/lC_harmonized_thickness_reduced.csv", row.names = FALSE)

######################### NORMATIVE MODEL ##################################

f_normative_model <- function(i){
  data_pat <- data_pat[, c(reg_names[i], "age", "scode", "subID", "euler","timepoint")]
  data_hc  <- data_hc[ , c(reg_names[i], "age", "scode", "subID", "euler","timepoint")]
  colnames(data_pat) <- c("reg", "age", "scode", "subID", "euler","timepoint")
  colnames(data_hc)  <- c("reg", "age", "scode", "subID", "euler","timepoint")
  
  # NORMATIVE MODEL
  model = lme(as.formula("reg ~ scode + age + euler"), random = ~1 + age|subID, 
              data = data_hc, control=lmeControl(opt='optim', msMaxIter=150), 
              method='REML')
  
  # PREDICTED VALUES FOR CONTROLS
  predictions <- predict(model, data_hc, na.action = na.omit, interval = "confidence")
  
  # COMPUTE Zs FOR CONTROLS
  for (hc in 1:nrow(data_hc)){
    Z_hc <- (data_hc$reg[hc] - predictions[hc]) / model$sigma
    df_Z1 = rbind(df_Z1, data.frame("region" = reg_names[i], 
                                    "subID" = data_hc$subID[hc], 
                                    "timepoint" = data_hc$timepoint[hc],
                                    "z" = Z_hc,
                                    "group" = 'control'))
  }
  
  # DYNAMIC PREDICTIONS FOR PATIENTS
  for (sub in pat_subIDs){
    data_sub <- subset(data_pat, subID == sub)
    prediction_sub <- IndvPred_lme(model, newdata = data_sub, all_times = F,
                                   timeVar = "age", M = 500, return_data = T)
    # COMPUTE Zs FOR PATIENTS
    n_tp <- nrow(data_sub)
    for (tp in 1:n_tp){
      Z <- (data_sub[tp,1] - prediction_sub$pred[tp]) / model$sigma
      df_Z2 = rbind(df_Z2, data.frame("region"=reg_names[i], "subID"=sub, "timepoint"=tp, "z"=Z,"group"='patient'))
    }
  }
  df_Z = rbind(df_Z1, df_Z2)
  return(df_Z)
}

# 11.5 min mclapply /vs/ +2h lapply
start_time <- Sys.time()
df_Z1 = NULL
df_Z2 = NULL
df_Z <- do.call(rbind, mclapply(1:length(reg_names), 
                                FUN= f_normative_model, 
                                mc.cores=ncores))
end_time <- Sys.time()
end_time - start_time

#write.csv(df_Z,"/data_J/results/Zs.csv", row.names = FALSE)
rm(df_Z, df_Z1, df_Z2, start_time, end_time, data_hc, data_pat, pat_subIDs, idx_pat, idx_hc)
Zs <- read.csv("/data_J/results/Zs.csv", header = T)
########################## COMPUTE DEVIANTS ####################################

# 1. Number of deviant subjects (>1.96 or <-1.96) (per region per timepoint)
# we create a new column to show if (subID, region, tp) is deviant (1) or not (0).
Z1 <- Zs %>% 
  mutate(deviant = ifelse(z < -1.96 | z > 1.96, 1, 0))%>%
  group_by(region, timepoint,group) %>% 
  summarise(dev= sum(deviant==1), non_dev= sum(deviant==0))

# 2. Number of subjects who are non-deviant at BL (tp=1) but deviant at FU (tp>1)

Z2 <- Zs %>% 
  mutate(deviant = ifelse(z < -1.96 | z > 1.96, 1, 0))%>%
  group_by(subID, region) %>% 
  summarise(nondev_bl=sum(deviant==0 & timepoint==1),
            dev_fu=sum(deviant==1 & timepoint>1)) %>% 
  group_by(subID)%>%
  summarise(nreg = sum(nondev_bl*dev_fu))

n <- sum(Z2$nreg)   # 3.119 timepoints
n <- sum(Z2$nreg>0) # 417 de 467 subjects

# Calculate # of times each subject is deviant for each region
Z3 <- Zs %>% 
  mutate(deviant = ifelse(z < -1.96 | z > 1.96, 1, 0))%>%
  group_by(region, subID) %>% 
  summarise(ndev=sum(deviant==1),non_dev=sum(deviant==0))
  # filter(ndev==3) # Nb of subjects deviant 3 times

# 3. Get subID and region of subjects non-deviant for all timepoints
# To create the "age*dcode without deviants" model
Z4 <- Z3 %>% 
  filter(ndev==0)

rm(Z1, Z2, Z3, n)


########################### AGE*DIAGNOSIS WITH ALL CASES #######################

f_age_dcode <- function(i){
  data_hc_pat <- df[, c(reg_names[i], "age", "scode", "subID", "euler", "dcode", "dcode_age")]
  colnames(data_hc_pat) <- c("reg", "age", "scode", "subID", "euler", "dcode", "dcode_age")
  
  model2 = lme(as.formula("reg ~ scode + age + euler + dcode + dcode_age"), 
               random = ~1 + age|subID, data = data_hc_pat,
               control=lmeControl(opt='optim', msMaxIter=150), method='REML')
  pval <- anova(model2)$p
  return(pval)
}
Func <- function(i)  p.adjust(p_values_1[,i], method='fdr')

# <30s
start_time <- Sys.time()
p_values_1 <- do.call(rbind, mclapply(1:length(reg_names), FUN= f_age_dcode, mc.cores=ncores))
colnames(p_values_1)<-c("reg", "scode","age", "euler", "dcode", "dcode_age")
end_time <- Sys.time()
end_time - start_time
rm(start_time, end_time)

# p-values with FDR correction
p_values_1_fdr <- t(do.call(rbind, lapply(1:ncol(p_values_1), Func)))
colnames(p_values_1_fdr)<-c("reg", "scode","age", "euler", "dcode", "dcode_age")
#sum(p_values_2[,5]<0.05)
# PLOT p-values
par(mfrow = c(2, 1), mai = c(0.8, 0.6, 0.2, 0.1))  # 2 rows and 1 column
x1<-p_values_1[,6]
x2<-p_values_1_fdr[,6]
hist(x1, main="pvalues", ylim=c(0,150), col="orange")
hist(x2, main="pvalues with FDR", ylim=c(0,150), col="green")

######################## AGE*DIAGNOSIS WITHOUT DEVIANTS #######################
# We take the non-deviant subjects (-1.96 -  1.96) for each region and timepoint

f_age_dcode_nondev <- function(i){
  subIDs <- Z4$subID[Z4$region==reg_names[i]] # (x subIDs)/467
  data_nondev <- df %>% 
    select(reg_names[i], age, scode, subID, euler, dcode, dcode_age) %>% 
    filter(subID %in% subIDs)
  colnames(data_nondev) <- c("reg", "age", "scode", "subID", "euler", "dcode", "dcode_age")
  
  model = lme(as.formula("reg ~ scode + age + euler + dcode + dcode_age"),
              random = ~1 + age|subID, data = data_nondev,
              control=lmeControl(opt='optim', msMaxIter=150), method='REML')
  
  pval <- anova(model)$p
  return(pval)
}
Func <- function(i)  p.adjust(p_values_2[,i], method='fdr')

# <30s
start_time <- Sys.time()
p_values_2 <- do.call(rbind, mclapply(1:length(reg_names), FUN= f_age_dcode_nondev, mc.cores=ncores))
colnames(p_values_2)<-c("reg", "scode","age", "euler", "dcode", "dcode_age")
end_time <- Sys.time()
end_time - start_time
rm(start_time, end_time)

p_values_2_fdr <- t(do.call(rbind, lapply(1:ncol(p_values_2), Func)))
colnames(p_values_2_fdr)<-c("reg", "scode","age", "euler", "dcode", "dcode_age")

# PLOT p-values
# reg, scode, age, euler, dcode, dcode_age
par(mfrow = c(2, 1), mai = c(0.8, 0.6, 0.2, 0.1))  # 2 rows and 1 column
x1<-p_values_2[,5]
x2<-p_values_2_fdr[,5]
hist(x1, main="pvalues", col="orange")
hist(x2, main="pvalues with FDR", col="green")

####################### GLOBAL RATIOS #######################
# 4. Global Z score ratio (per subject per timepoint): Σ∣Z∣>2 / Σ∣Z∣<2
Z5 <- Zs %>% 
  group_by(subID, timepoint)%>%
  summarise(num = sum(abs(z)>2), den = sum(abs(z)<2)) %>%
  summarise(timepoint = timepoint, subID=subID, ratio = num/den)

# 5. Global Z score ratio (per subject per timepoint): code each region as 
# > 1.96 (1), 
# 1.96 - -1.96 (0),
# < -1.96 (-1) 
# and calculate average over all regions

Z6 <- Zs %>% 
  mutate(z = ifelse(z > -1.96 & z<1.96, 0, z))%>%
  mutate(z = ifelse(z > 1.96, 1, z))%>%
  mutate(z = ifelse(z < -1.96, -1, z))%>% 
  group_by(subID, timepoint)%>%
  summarise(ratio = mean(z))

# 6. Percentage positive/negative deviation (per subject per timepoint): 
# (#regions >1.96 / 308) * 100 and (#regions <-1.96 / 308) * 100

Z7 <- Zs %>% 
  mutate(z = ifelse(z > -1.96 & z<1.96, 0, z))%>%
  mutate(z = ifelse(z > 1.96, 1, z))%>%
  mutate(z = ifelse(z < -1.96, -1, z))%>% 
  group_by(subID, timepoint)%>%
  summarise(zpos = sum(z==1), zneg = sum(z==-1), z0 = sum(z==0))%>%
  summarise(per_pos = (zpos/308)*100, per_neg = (zneg/308)*100, per_0 = (z0/308)*100)

# 7. Global score (per region) across subjects: (#deviating subjects / total #subjects) * 100
# DUDA: Se tienen en cuenta los timepoints: un sujeto puede ser devi para 1 tp y nondev para 2 tp.
# total #subjects = 1087 (unique subjects = 467)

########################## CORRELATIONS ##################################
# 9. Correlations between symptoms, IQ, and changes in global scores

# SYMPTOMS
df$PANSS_total

df$Total_P
df$Total_N
df$Total_G

df$P1delusions
df$N1blunted_affect
df$G1somatic_concern

# IQ 
df$IQ_date
df$IQ_study_code
df$IQ_total

###################### ONE SAMPLE T-TEST ######################################
# One sample t-test for each region (per timepoint per diagnostic group) 
# (across subjects): is Z on average non-zero?. FDR necessary.

install.packages("ggpubr")
t.test(data, mu = 0, alternative = "two.sided")
shapiro.test(data) # n>30 asique no hace falta mirar si tiene distr Normal




