Regression_NormativeModel <- function(i, df, measure){
  # df <- df_lC_NO_matched
  # measure = "CT_freesurfer"
  # parc = "parc68"
  # match = "NOmatch"
  # harmonization = "lC"
  # i=1
  # sub = 100
  # tp = 1
  
  df_Z <- df_Z1 <- df_Z2 <- NULL
  
  # Select data
  reg_names  <- colnames(df%>%select(ends_with(measure)))
  idx_pat    <- which(df$dcode== 1)  # 1: patient
  idx_hc     <- which(df$dcode== 0)  # 0: healthy control
  data_pat   <- df[idx_pat, ]        # Select rows corresponding to patients
  data_hc    <- df[idx_hc,  ]        # Select rows corresponding to healthy controls
  pat_subIDs <- unique(data_pat$subID)   # Get unique subIDs from patients

  var_names <- c("age", "scode", "subID", "euler", "timepoint")
  data_pat <- data_pat[, c(reg_names[i], var_names)]
  data_hc  <- data_hc[ , c(reg_names[i], var_names)]
  colnames(data_pat) <- c("reg", var_names)
  colnames(data_hc)  <- c("reg", var_names)
  
  # NORMATIVE MODEL
  model <- lme(as.formula("reg ~ scode + age + euler"), random = ~1 + age|subID,
              data = data_hc, control=lmeControl(opt='optim', msMaxIter=200),
              method='REML')
  
  # PREDICTED VALUES FOR CONTROLS
  predictions <- predict(model, data_hc, na.action = na.omit, interval = "confidence")
  
  # COMPUTE Zs FOR CONTROLS
  for (hc in 1:nrow(data_hc)){
    # hc = 1
    Z_hc <- (data_hc$reg[hc] - predictions[hc]) / model$sigma
    df_Z1 <- rbind(df_Z1, data.frame("region" = reg_names[i],
                                     "subID"  = data_hc$subID[hc],
                                     "timepoint" = data_hc$timepoint[hc],
                                     "z"     = Z_hc,
                                     "group" = 'control'))
  }

  # DYNAMIC PREDICTIONS FOR PATIENTS
  pat_subIDs <- pat_subIDs[pat_subIDs!=949] # This patient gives problems when computing IndvPred (we delete it for now)
  for (sub in pat_subIDs){
    data_sub <- subset(data_pat, subID == sub)
    prediction_sub <- IndvPred_lme(model, newdata = data_sub, 
                                   all_times = F, timeVar = "age", 
                                   M = 500, return_data = T)
    # COMPUTE Zs FOR PATIENTS
    n_tp <- nrow(data_sub)
    for (tp in 1:n_tp){
      Z  <- (data_sub[tp,1] - prediction_sub$pred[tp]) / model$sigma
      df_Z2 <- rbind(df_Z2, data.frame("region" = reg_names[i], 
                                       "subID"  = sub, 
                                       "timepoint" = tp, 
                                       "z"    = Z,
                                       "group"='patient'))
    }
  }
  df_Z <- rbind(df_Z1, df_Z2)
  return(df_Z)
}

run_NormativeModel <- function(df, measure, parc, match, harmonization){
  # df = df_lC_NO_matched
  # measure = "CT_freesurfer"
  # parc = "parc68"
  # match = "NOmatch"
  # harmonization = "lC"

  filepath = paste0("/data_J/Results/Files/",parc,"/")
  filename = paste0("Zs_",harmonization,"_",match,"_",parc,".csv")
  
  if (file.exists(paste0(filepath,filename))){
    df_Z <- read.csv(paste0(filepath, filename), header = T)
  }
  
  else if (!file.exists(paste0(filepath,filename))){
    n <- length(colnames(df%>%select(ends_with(measure))))
    
    df_Z <- do.call(rbind, mclapply(1:n, Regression_NormativeModel,
                                    df, measure,
                                    mc.cores = ncores-1))
    
    write.csv(df_Z, paste0(filepath, filename), row.names = FALSE)
  }
  return(df_Z)
}

RegressionModel_AgeDiagnosis <- function(i, df, Z, measure, exclude_deviants){
  
  pval <- NULL
  var_names  <- c("age", "scode", "subID", "euler", "dcode", "dcode_age")
  reg_names  <- colnames(df%>%select(ends_with(measure)))
  
  if (!exclude_deviants){ # NOT EXCLUDING DEVIANTS
    data_hc_pat <- df[, c(reg_names[i], var_names)]
    colnames(data_hc_pat) <- c("reg", var_names)
    
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
  
  return(pval)
}

run_AgeDiagnosisModel <- function(df, measure, Z_dev, exclude_deviants){

  # df <- df_lC_matched
  # measure <- "CT_freesurfer"
  # Z_dev <- Zs_match
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



# PRUEBAS PARA COHENS'D

Regression_NormativeModel_Cohensd <- function(i, df, measure){
  # df <- df_lC_matched_308
  # measure = "thickness"
  # parc = "parc308"
  # match = "match"
  # harmonization = "lC"
  # i=1
  # sub = 100
  # tp = 1
  
  df_Z <- df_Z1 <- df_Z2 <- NULL
  
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
  # df = df_lC_NO_matched
  # measure = "CT_freesurfer"
  # parc = "parc68"
  # match = "NOmatch"
  # harmonization = "lC"
  
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


