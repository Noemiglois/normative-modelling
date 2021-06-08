library(lme4)
library(JMbayes)
library(abind)
library(dplyr)
library(lmerTest)
library(bestNormalize)

# DATASET
df <- read.csv("/data_J/data/morphosim_lC_harmonized_data.csv", header = T)
df <- transform(df, euler=(euler_lh+euler_rh)/2) # Create euler variable
df <- df[order(df$age), ]                        # Sort by age ¿WHY?
df <- transform(df, age_squared=age^2)           # Create age^2 variable
# test <- bestNormalize(df$euler)
# df$euler_norm <- predict(test)
# df$euler <- scale(df$euler, scale = F) # Scale euler 
# shapiro.test(df$euler)

df_Z = NULL

idx_pat <- which(df$dcode== 1) # 1 pat
idx_hc  <- which(df$dcode== 0) # 0 hc

reg_names <- colnames(df%>%select(ends_with("thickness")))[1:3]

start_time <- Sys.time()
for (i in 1:length(reg_names)){
  cat(i,'...')
  # vars <- c(reg_names[i], "age", "scode", "age_squared", "subID","euler")
  vars <- c(reg_names[i], "age", "scode", "subID", "euler")
  data_pat <- df[idx_pat, vars]
  data_hc  <- df[idx_hc, vars]

  feat <- noquote(reg_names[i])
  fixed_effects <- paste(feat, " ~ scode + age + euler")    
  # fixed_reduced <- paste(feat, " ~ scode + age")                # fixed effects 
  # fixed_full <- paste(feat, " ~ scode + age + age_squared")     # fixed effects 

  model = lme(formula(fixed_effects), random = ~1 + age|subID, data = data_hc, control=lmeControl(opt='optim', msMaxIter=150), method='ML') # opt='optim', msMaxIter=150
  # model_reduced = lme(formula(fixed_reduced), random = ~1 + age|subID, data = data_hc, control=lmeControl(opt='optim', msMaxIter=150), method='ML') # opt='optim', msMaxIter=150
  # model_full = lme(formula(fixed_full), random = ~1 + age + age_squared |subID, data = data_hc, control=lmeControl(opt='optim', msMaxIter=150), method='ML')
  
  # pval <- anova(model_reduced, model_full)$p[2]
  # if(pval<0.05) {
  #   m = 1  # Full model
  #   model = lme(formula(fixed_full), random = ~1+ poly(age, degree = 2) |subID, data = data_hc, control=lmeControl(opt='optim', msMaxIter=150), method='REML')
  # } else {
  #   m = 0 # Reduced model
  #   model = lme(formula(fixed_reduced), random = ~ 1 + age|subID, data = data_hc, control=lmeControl(opt='optim', msMaxIter=150), method='REML')
  # }

  #pred <- predict(model, data_pat, level = 0, na.action = na.omit, interval = "confidence")

  subIDs <- unique(data_pat$subID)
  for (sub in subIDs){
    data_sub <- subset(data_pat, subID == sub)
    n_tp <- nrow(data_sub)
    prediction_sub <- IndvPred_lme(model, newdata = data_sub, all_times = F,
                                   timeVar = "age", M = 500, return_data = T)
    for (tp in 1:n_tp){
      Z <- (data_sub[tp,1] - prediction_sub$pred[tp]) / model$sigma
      df_Z = rbind(df_Z, data.frame("region"=reg_names[i], "subID"=sub, "timepoint"=tp, "z"=Z))
      # df_Z = rbind(df_Z, data.frame("region"=reg_names[i], "subID"=sub, "timepoint"=tp, "z"=Z, "model"= m))
    }
  }
}
end_time <- Sys.time()
end_time - start_time

write.csv(df_Z,"/data_J/results/Zs_lme.csv", row.names = FALSE)