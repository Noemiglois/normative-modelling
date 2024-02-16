FDR_correction_pval <- function(i, p_values){
  p_values_fdr <- p.adjust(p_values[,i], method='fdr') 
  return(p_values_fdr)
}

Apply_FDR_Correction <- function(p_val){
  p_val_FDR <- t(do.call(rbind, mclapply(1:ncol(p_val),
                                         FDR_correction_pval,
                                         p_val,
                                         mc.cores=ncores-1)))
  colnames(p_val_FDR) <- c("reg", "scode","age", "euler", "dcode", "dcode_age")
  return(p_val_FDR)
}


#### PRUEBAS PARA COHENS'D ######
Regression_NormativeModel_Cohensd <- function(i, df, measure){
  
  # Select data
  reg_names  <- colnames(df%>%select(ends_with(measure)))
  diag <- df[ ,"dcode"] 

  var_names <- c("age", "scode", "subID", "euler", "timepoint")
  data <- df[, c(reg_names[i], var_names)]
  colnames(data) <- c("reg", var_names)

  # NORMATIVE MODEL
  model <- lme(as.formula("reg ~ scode + age + euler"), random = ~1 + age|subID,
               data = data, control=lmeControl(opt='optim', msMaxIter=200),
               method='REML')
  
  res  <- model$residuals[,"fixed"]
  res  <- model$residuals[,"subID"]
  d <- cohens_d(res ~ diag)
  #print(d, append_CL = TRUE)                   # Easy interpretation
  #interpret_d(d$Cohens_d, rules = "cohen1988") # Easy interpretation
  d$interpret <- interpret_d(d$Cohens_d, rules = "cohen1988")[1]
  return(d)
}

run_NormativeModel_Cohensd <- function(df, measure, parc, match, harmonization){
  
  filepath = paste0("/data_J/Results/Files/",parc,"/")
  filename = paste0("CohensD_",harmonization,"_",match,"_",parc,".csv")
  
  if (file.exists(paste0(filepath,filename))){
    df_d <- read.csv(paste0(filepath, filename), header = T)
  }
  
  else if (!file.exists(paste0(filepath,filename))){
    n <- length(colnames(df%>%select(ends_with(measure))))
    
    df_d <- do.call(rbind, mclapply(1:n, Regression_NormativeModel_Cohensd,
                                    df, measure,
                                    mc.cores = ncores-1))
    
    write.csv(df_d, paste0(filepath, filename), row.names = FALSE)
  }
  return(df_d)
}

########### Pruebas lme con train-test #############

NormativeModel <- function(i, df, measure){

  # Select data
  split <- train_test_split_f(df, 
                              perc = 0.9)
  df_train <- split[[1]]
  df_test  <- split[[2]]
  df_pat   <- split[[3]]
  
  reg_name  <- colnames(df%>%dplyr::select(ends_with(measure)))
  # cat("Noemi", reg_name[i], i, "\n") # Borrar
  var_names <- c("age", "scode", "subID", "euler", "timepoint")
  df_train  <- df_train[, c(reg_name[i], var_names)]
  df_test   <- df_test[,  c(reg_name[i], var_names)]
  df_pat    <- df_pat[,   c(reg_name[i], var_names)]
  
  colnames(df_train) <- c("reg", var_names)
  colnames(df_test)  <- c("reg", var_names)
  colnames(df_pat)   <- c("reg", var_names)
  
  df_Z <- NULL
  
  # NORMATIVE MODEL
  model <- lme(as.formula("reg ~ scode + age + euler"), random = ~1 + age|subID,
               data = df_train, control=lmeControl(opt='optim', msMaxIter=300),
               method='REML')
  
  # PREDICTED VALUES FOR TRAIN CONTROLS
  predictions <- predict(model, df_train, na.action = na.omit, interval = "confidence")
  
  # COMPUTE Zs FOR TRAIN CONTROLS
  for (hc in 1:nrow(df_train)){
    # hc = 1
    Z <- (df_train$reg[hc] - predictions[hc]) / model$sigma
    df_Z <- rbind(df_Z, data.frame("region" = reg_name[i],
                                     "subID"  = df_train$subID[hc],
                                     "timepoint" = df_train$timepoint[hc],
                                     "z"     = Z,
                                     "group" = 'train_control',
                                     "sigma" = model$sigma))
  }
  
  # DYNAMIC PREDICTIONS FOR TEST CONTROLS
  test_subIDs <- unique(df_test$subID)   # Get unique subIDs from test controls
  test_subIDs <- test_subIDs[test_subIDs!=655] # This control gives problems when computing IndvPred (we delete it for now)
  test_subIDs <- test_subIDs[test_subIDs!=567] # This control gives problems when computing IndvPred (we delete it for now)
  
  for (sub in test_subIDs){
    # sub=14
    data_sub <- subset(df_test, subID == sub)
    prediction_sub <- IndvPred_lme(model, newdata = data_sub, 
                                   all_times = F, timeVar = "age", 
                                   M = 500, return_data = T)
    # COMPUTE Zs FOR TEST CONTROLS
    n_tp <- nrow(data_sub)
    for (tp in 1:n_tp){
      # tp=2
      Z  <- (data_sub[tp,1] - prediction_sub$pred[tp]) / model$sigma
      df_Z <- rbind(df_Z, data.frame("region"   = reg_name[i], 
                                       "subID"  = sub, 
                                       "timepoint" = tp, 
                                       "z"    = Z,
                                       "group"='test_control',
                                       "sigma" = model$sigma))
    }
  }
  
  # DYNAMIC PREDICTIONS FOR PATIENTS
  pat_subIDs <- unique(df_pat$subID)   # Get unique subIDs from patients
  pat_subIDs <- pat_subIDs[pat_subIDs!=949] # This patient gives problems when computing IndvPred (we delete it for now)
  for (sub in pat_subIDs){
    # cat(sub, "\n") # borrar
    data_sub <- subset(df_pat, subID == sub)
    prediction_sub <- IndvPred_lme(model, newdata = data_sub, 
                                   all_times = F, timeVar = "age", 
                                   M = 500, return_data = T)
    # COMPUTE Zs FOR PATIENTS
    n_tp <- nrow(data_sub)
    for (tp in 1:n_tp){
      Z  <- (data_sub[tp,1] - prediction_sub$pred[tp]) / model$sigma
      df_Z <- rbind(df_Z, data.frame("region" = reg_name[i], 
                                       "subID"  = sub, 
                                       "timepoint" = tp, 
                                       "z"    = Z,
                                       "group"='test_patient',
                                       "sigma" = model$sigma))
    }
  }
  return(df_Z)
  #return(list(df_Z, model$sigma))
}

run_NormativeModel <- function(df, measure, parc, match, harmonization){

  filepath = paste0("/data_J/Results/Files/",parc,"/")
  filename = paste0("Zs_",harmonization,"_",match,"_",parc,".csv")
  
  if (file.exists(paste0(filepath,filename))){
    df_Zs <- read.csv(paste0(filepath, filename), header = T)
  }
  
  else if (!file.exists(paste0(filepath,filename))){
    n <- length(colnames(df%>%select(ends_with(measure))))

    # df_Zs <- do.call(rbind, lapply(1:n, NormativeModel,
    #                                 df, measure))
    
    df_Zs <- do.call(rbind, mclapply(1:n, NormativeModel,
                                     df, measure,
                                     mc.cores = ncores-1))
    write.csv(df_Zs, paste0(filepath, filename), row.names = FALSE)
  }
  return(df_Zs)
}

##############prueba#################

RegressionModel_AgeDiagnosis <- function(i, df, Z, measure, exclude_deviants){
  
  pval <- NULL
  var_names  <- c("age", "scode", "subID", "euler", "dcode", "dcode_age")
  reg_names  <- colnames(df%>%select(ends_with(measure)))
  
  if (!exclude_deviants){ # NOT EXCLUDING DEVIANTS
    data_hc_pat <- df[, c(var_names, reg_names[i])]
    colnames(data_hc_pat) <- c(var_names, "reg")
    
    model <- lme(as.formula("reg ~ scode + age + euler + dcode + dcode_age"), 
                 random = ~ 1 + age|subID, data = data_hc_pat,
                 control= lmeControl(opt='optim', msMaxIter=150), 
                 method = 'REML')
    pval <- anova(model)$p 
  }
  else if (exclude_deviants){ # EXCLUDING DEVIANTS
    subIDs <- Z$subID[Z$region==reg_names[i]] # (x subIDs)/467
    data_nondev <- df %>% 
      select(reg_names[i], age, scode, subID, euler, dcode, dcode_age) %>% 
      filter(subID %in% subIDs)
    
    colnames(data_nondev) <- c("reg", var_names)
    
    model <- lme(as.formula("reg ~ scode + age + euler + dcode + dcode_age"),
                 random = ~1 + age|subID, data = data_nondev,
                 control=lmeControl(opt='optim', msMaxIter=150), method='REML')
    pval <- anova(model)$p
  }
  # Añadir un tercer modelo que sólo incluya dcode para aquellas regiones 
  # en las que la variable dcode_age no es significativa
  return(pval)
}

run_AgeDiagnosisModel <- function(df, measure, Z_dev, exclude_deviants){

  n <- length(colnames(df %>% select(ends_with(measure))))
  
  if (exclude_deviants){ # EXCLUDING DEVIANTS
    Z_nondev <- Z_dev %>%
      mutate(dev = ifelse(z < -1.96 | z > 1.96, 1, 0)) %>%
      filter(dev==0) # Filter the subjects that are not deviates
  }
  else if (!exclude_deviants){ # NOT EXCLUDING DEVIANTS
    Z_nondev <- Z_dev
  }
  
  p_val <- do.call(rbind, mclapply(1:n,
                                   RegressionModel_AgeDiagnosis,
                                   df = df,
                                   Z = Z_nondev,
                                   measure = measure,
                                   exclude_deviants = exclude_deviants,
                                   mc.cores = ncores-1))
  
  colnames(p_val) <- c("reg", "scode", "age", "euler", "dcode", "dcode_age")
  
  return(p_val)
}