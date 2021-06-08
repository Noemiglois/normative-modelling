library(lme4)
library(JMbayes)
library(abind)
library(dplyr)
library(lmerTest)

# DATASET
df <- read.csv("/data_J/data/morphosim_lC_harmonized_data.csv", header = T)
# df <- transform(df, age_squared=age^2)           # Create age^2 variable
df <- transform(df, euler=(euler_lh+euler_rh)/2) # Create euler variable
df <- df[order(df$age), ]                        # Sort by age ¿WHY?
euler <- df$euler
test <- bestNormalize(euler)
p <- predict(test)
df$euler_norm <- p

df_Z = NULL

idx_pat <- which(df$dcode== 1) # 1 pat
idx_hc  <- which(df$dcode== 0) # 0 hc

reg_names <- colnames(df%>%select(ends_with("thickness")))
for (i in 1:length(reg_names)){
  cat(i,'...')
  vars <- c(reg_names[i], "age", "acode", "scode", "dcode", "subID", "euler_norm")
  data_pat <- df[idx_pat, vars] 
  data_hc  <- df[idx_hc, vars] 
  
  feat <- noquote(reg_names[i])
  form <- paste(feat, " ~ scode + poly(age, degree = 2) + euler_norm + (1|subID)") # fixed + random effects
  
  model = lmer(formula(form),
                    control = lmerControl(optimizer = "nloptwrap"),
                    data = data_hc)
  pval <- anova(model)$P[3] # if agesq not significant then exclude it:
  if(pval<0.05) {
    m = 1 # Full model
    form <- paste(feat, " ~ scode + poly(age, degree = 2) + euler_norm + (1|subID)")
    model1 = lmer(formula(form),
                 control = lmerControl(optimizer = "nloptwrap"),
                 data = data_hc)
    
    form <- paste(feat, " ~ scode + poly(age, degree = 2) + euler_norm")
    model2 = lme(formula(form), random = ~1 + poly(age, degree = 2)|subID, data = data_hc, 
                control=lmeControl(opt='optim'), 
                method='REML')
    pred1<-predict( model1, data_hc, na.action = na.omit)
    pred2<-predict( model2, data_hc, na.action = na.omit)
    correl = cor.test(pred1,pred2)$estimate
    
  } else {
    m = 0  # Reduced model
    form <- paste(feat, " ~ scode + age + euler_norm + (1|subID)")
    model1 = lmer(formula(form),
                  control = lmerControl(optimizer = "nloptwrap"),
                  data = data_hc)
    
    form <- paste(feat, " ~ scode + age + euler_norm")
    model2 = lme(formula(form), random = ~ 1 + age |subID, data = data_hc, 
                control=lmeControl(opt='optim'), 
                method='REML')
    pred1<-predict( model1, data_hc, na.action = na.omit)
    pred2<-predict( model2, data_hc, na.action = na.omit)
    correl = cor.test(pred1,pred2)$estimate
  }
  
  subIDs <- unique(data_pat$subID)
  for (sub in subIDs){
    data_sub <- subset(data_pat, subID == sub)
    # Dynamic predictions for new subject (patient) using model2 (lme)
    pred_sub <- IndvPred_lme(model2, newdata = data_sub, all_times = F,
                                   timeVar = "age", M = 500, return_data = T) # Gives predictions at ages after true ages
    n_tp <- nrow(data_sub)
    for (tp in 1:n_tp){
      Z <- (data_sub[tp,1] - pred_sub$pred[tp]) / model2$sigma
      df_Z = rbind(df_Z, data.frame("region"=reg_names[i], "model"= m, "corr pred"=correl, "subID"=sub, "timepoint"=tp, "z"=Z))
    }
  }
}

write.csv(df_Z,"/data_J/results/Zs_lmer.csv", row.names = FALSE)