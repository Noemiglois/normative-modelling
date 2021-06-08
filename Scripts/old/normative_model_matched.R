############################# DATASET  ########################################

df <- read.csv("/data_J/Data/lC_harmonized_thickness_matched.csv", header = T)

idx_pat <- which(df$dcode== 1)  # 1: patient
idx_hc  <- which(df$dcode== 0)  # 0: healthy control
data_pat <- df[idx_pat, ]
data_hc  <- df[idx_hc, ]

pat_subIDs <- unique(data_pat$subID)
ncores = detectCores()
reg_names <- colnames(df%>%select(ends_with("thickness"))) # Select thickness features names

######################### NORMATIVE MODEL ##################################
f_normative_model <- function(i){
  data_pat <- data_pat[, c(reg_names[i], "age", "scode", "subID", "euler", "timepoint")]
  data_hc  <- data_hc[ , c(reg_names[i], "age", "scode", "subID", "euler", "timepoint")]
  colnames(data_pat) <- c("reg", "age", "scode", "subID", "euler", "timepoint")
  colnames(data_hc)  <- c("reg", "age", "scode", "subID", "euler", "timepoint")
  
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

start_time <- Sys.time()
df_Z1 = NULL
df_Z2 = NULL
df_Z <- do.call(rbind, mclapply(1:length(reg_names), 
                                FUN= f_normative_model, 
                                mc.cores=ncores))
end_time <- Sys.time()
end_time - start_time

write.csv(df_Z,"/data_J/Results/Zs_matched.csv", row.names = FALSE)
rm(df_Z, df_Z1, df_Z2, start_time, end_time, data_hc, data_pat, pat_subIDs, idx_pat, idx_hc)
Zs <- read.csv("/data_J/results/Zs_matched.csv", header = T)







