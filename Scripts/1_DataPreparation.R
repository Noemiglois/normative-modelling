getData <- function(parc){
  # This function will combine the original datasets into the one 
  # that will be using for analysis.
  
  # CREATION OF THE 68 REGIONS PARCELLATION DATASET
  if (parc == "68"){ 
    parc = "parc68"
    if (!file.exists(paste0("/data_J/Data/parc68/",parc,".csv"))){
      cat("File parc68 doesn't exist, creating...\n")
      
      # Measures from both hemispheres are combined. 
      # # IGI
      # df_lh_Igi   <- read.csv("/data_J/Data/metrics/lh.aparc.lgi.csv")
      # df_rh_Igi   <- read.csv("/data_J/Data/metrics/rh.aparc.lgi.csv")
      # df_Igi <- merge(df_lh_Igi, df_rh_Igi, by="ID", all.x=T)
      # rm(df_lh_Igi, df_rh_Igi)
      
      # CT
      df_lh_CT   <- read.csv("/data_J/Data/metrics/lh.aparc.thickness.csv")[,1:35]
      df_rh_CT   <- read.csv("/data_J/Data/metrics/rh.aparc.thickness.csv")[,1:35]
      df_CT <- merge(df_lh_CT, df_rh_CT, by="ID", all.x=T)
      rm(df_lh_CT, df_rh_CT)
      
      # Vol
      df_lh_Vol   <- read.csv("/data_J/Data/metrics/lh.aparc.volume.csv")
      df_rh_Vol   <- read.csv("/data_J/Data/metrics/rh.aparc.volume.csv")
      df_Vol <- merge(df_lh_Vol, df_rh_Vol, by="ID", all.x=T)
      rm(df_lh_Vol, df_rh_Vol)
      
      df <- merge(df_CT, df_Vol, by="ID", all.x=T)
      rm(df_Vol, df_CT)
      
      # Covariates dataframe is also created by selecting those relevant ones: 
      df_morphosim <- read.csv("/data_J/Data/planU_morphosim_rawdatabase.csv")
      df_euler     <- read.csv("/data_J/Data/planU_euler.csv")
      temp   <- as.vector(df_euler[, 3])
      df_euler[,3] <- impute(temp, fun=median)
      df_raw <- merge(df_morphosim, df_euler, by="ID", all.x=T)
      df_raw <- df_raw[, c("age", "scode", "dcode", "acode", "timepoint",
                           "subID", "euler_lh", "euler_rh","ID")]
      
      # Final dataset combining <CT measures> and <covariates>:
      df <- merge(df_raw, df, by="ID", all.x=T)
      
      # Columns creation: euler and dcode_age
      df <- transform(df, euler=(euler_lh+euler_rh)/2) # Create euler variable
      df <- transform(df, dcode_age=dcode*age)         # Create dcode_age variable
      df <- transform(df, age2=age**2)                 # Create age squared variable
      df <- transform(df, age3=age**3)                 # Create age cubed variable
      df <- df[order(df$age), ]                        # Sort by age (not sure why)
      df[ ,c('euler_rh', 'euler_lh')] <- list(NULL)
      
      # Delete subjects with tp 1 and tp 3, but not tp 2 (project criteria)
      # Alternative: change tp3 by tp2
      # Alternative: not delete those subjects --> ¿what would happen?
      df <- df[df$subID != 133, ] 
      df <- df[df$subID != 267, ] 
      df <- df[df$subID != 804, ] 
      df <- df[df$subID != 822, ] 
      df <- df[df$subID != 879, ] 
      df <- df[df$subID != 931, ] 
      
      rm(df_euler, df_morphosim, df_raw, temp)
      
      metric_reg_names <- sub("CT_freesurfer", "CT",colnames(df)) # Delete termination to get a cleaner name
      colnames(df) <- metric_reg_names
      metric_reg_names <- sub("VOL_freesurfer", "Vol",colnames(df)) # Delete termination to get a cleaner name
      colnames(df) <- metric_reg_names
      rm(metric_reg_names)
      
      feat_names  <- c(colnames(df%>%dplyr::select(ends_with("CT"))), 
                       colnames(df%>%dplyr::select(ends_with("Vol"))))
      cov_names   <- setdiff(colnames(df), feat_names)
      df <- df[,c(sort(cov_names), feat_names)]
      rm(cov_names, feat_names)
      
      write.csv(df, paste0("/data_J/Data/parc68/", parc,".csv"), row.names = FALSE)
      cat("Done!\n")}
    
    else if (file.exists(paste0("/data_J/Data/parc68/",parc,".csv"))){
      cat("File parc68 already exists, reading file...\n")
      df <- read.csv(paste0("/data_J/Data/parc68/",parc,".csv"))
    }
  }
  
  # CREATION OF THE 308 REGIONS PARCELLATION DATASET
  if (parc == "308"){ 
    parc = "parc308"
    if (!file.exists(paste0("/data_J/Data/parc308/",parc,".csv"))){ 
      cat("File parc308 doesn't exist, creating....\n")
      
      # Measures (CT) and covariates (those relevant) are combined.
      df_morphosim <- read.csv("/data_J/Data/planU_morphosim_rawdatabase.csv")
      df_euler <- read.csv("/data_J/Data/planU_euler.csv")
      temp <- as.vector(df_euler[,3])
      df_euler[,3] <- impute(temp, fun=median)
      df_raw <- merge(df_morphosim, df_euler, by="ID", all.x=T)
      
      cov_names <- c("age", "scode", "dcode", "acode", "timepoint",
                     "subID", "euler_lh", "euler_rh", "ID")
      feat_names_CT  <- colnames(df_raw%>%dplyr::select(ends_with("thickness")))
      feat_names_Vol <- colnames(df_raw%>%dplyr::select(ends_with("volume")))
      
      df <- df_raw[, c(cov_names, feat_names_CT, feat_names_Vol)]
      
      # Columns creation: euler and dcode_age
      df <- transform(df, euler=(euler_lh+euler_rh)/2) # Create euler variable
      df <- transform(df, dcode_age=dcode*age)         # Create dcode_age variable
      df <- transform(df, age2=age**2)                 # Create age squared variable
      df <- transform(df, age3=age**3)                 # Create age cubed variable
      df <- df[order(df$age), ]                        # Sort by age (not sure why)
      df[ ,c('euler_rh', 'euler_lh')] <- list(NULL)
      
      # Delete subjects with tp 1 and tp 3, but not tp 2 (project criteria)
      # Alternative: change tp3 by tp2
      # Alternative: not delete those subjects --> ¿what would happen?
      df <- df[df$subID != 133, ] 
      df <- df[df$subID != 267, ] 
      df <- df[df$subID != 804, ] 
      df <- df[df$subID != 822, ] 
      df <- df[df$subID != 879, ] 
      df <- df[df$subID != 931, ] 
      
      rm(df_euler, df_morphosim, df_raw, cov_names, feat_names_CT, feat_names_Vol, temp)
      
      metric_reg_names <- sub("thickness", "CT",colnames(df)) # Delete termination to get a cleaner name
      colnames(df) <- metric_reg_names
      metric_reg_names <- sub("volume", "Vol",colnames(df)) # Delete termination to get a cleaner name
      colnames(df) <- metric_reg_names
      rm(metric_reg_names)
      
      feat_names  <- c(colnames(df%>%dplyr::select(ends_with("CT"))), 
                       colnames(df%>%dplyr::select(ends_with("Vol"))))
      cov_names   <- setdiff(colnames(df), feat_names)
      df <- df[,c(sort(cov_names), feat_names)]
      rm(cov_names, feat_names)
            
      write.csv(df,"/data_J/Data/parc308/parc308.csv", row.names = FALSE)
      cat("Done!\n")}
    
    else if (file.exists(paste0("/data_J/Data/parc308/",parc,".csv"))){ 
      cat("File parc308 already exists, reading file...\n") 
      df <- read.csv(paste0("/data_J/Data/parc308/", parc, ".csv"))
    }
  }
  return(df)
}

DataPreprocessing <- function(df, parc, measure, harmonization, match){
  # This function will prepare the df by harmonizing or matching.
  
  # HARMONIZATION WITH LONGITUDINAL COMBAT
  if (harmonization){ 
      cat("File with longCombat harmonization doesn't exist, creating...\n")
      measures <- c("CT", "Vol") 
      measure2 <- setdiff(measures, measure)
      
      feat_names <- colnames(df%>%select(ends_with(measure)))
      
      lC_measure <- longCombat(idvar    = 'subID', 
                               timevar  = 'timepoint',
                               batchvar = 'acode', 
                               features = feat_names, 
                               formula  = "age + dcode*timepoint",
                               ranef = "(1|subID)",
                               data  = df)
            
      data_lC_measure <- lC_measure$data_combat
      colnames(data_lC_measure)[4:ncol(data_lC_measure)] <- feat_names
      
      all_colnames <- colnames(df)
      data_lC_measure_colnames <- colnames(data_lC_measure)
      diff_colnames <- c(setdiff(all_colnames, data_lC_measure_colnames), 
                         setdiff(data_lC_measure_colnames, all_colnames))
      
      df_to_merge <- df[c(diff_colnames,"subID","timepoint")] # Subset df
      df <- merge(df_to_merge, data_lC_measure, by=c("subID","timepoint"))
      
      feat_names  <- c(colnames(df%>%dplyr::select(ends_with("CT"))), 
                       colnames(df%>%dplyr::select(ends_with("Vol"))))
      cov_names   <- setdiff(colnames(df), feat_names)
      feat_names <- colnames(df%>%select(ends_with(measure)))
      df <- df[,c(sort(cov_names), feat_names)]
      write.csv(df,paste0("/data_J/Data/", paste0("parc", parc),"/",paste0("parc", parc),paste0("_", measure),"_lC_harmonized.csv"), row.names = FALSE)  }
  
  # MATCH-IT ALGORITHM
  if (match){
      cat("File with MATCH-IT doesn't exist, creating...\n")
      set.seed(99) # Otherwise every match-it would lead to different samples
      
      # Baseline
      df_bl  <- df[which(df$timepoint==1), ] 
      temp1  <- df_bl[ , c("ID","age","scode","euler","dcode")] # sample to match
      temp2 <- matchit(dcode ~ age+scode, data=temp1) # matching method by default: nearest neighbor
      # formula: two-sided formula() object containing the treatment and 
      # the covariates to be used in creating the distance measure used in the matching. 
      # This formula will be supplied to the functions that estimate the distance measure. 
      # The formula should be specified as A ~ X1 + X2 + ... where 
      # A represents the treatment variable and X1 and X2 are covariates.
      temp3 <- match.data(temp2)[1:ncol(temp1)] # matched sample
      
      df_bl_matched  <- subset(df_bl,  ID %in% temp3$ID) # extracting matched sample for baseline
      
      # Follow up 1
      df_fu1  <- df[which(df$timepoint==2), ]
      df_fu1  <- subset(df_fu1,  subID %in% df_bl_matched$subID)  # Get only the subjects that were in bl matched
      
      temp1 <- df_fu1[ , c("ID","age","scode","euler","dcode")]   # samples to match
      temp2 <- matchit(dcode ~ age+scode, data=temp1 )   # matching
      temp3 <- match.data(temp2)[1:ncol(temp1)]          # matched sample
      
      df_fu1_matched  <- subset(df_fu1,  ID %in% temp3$ID) # extracting matched sample for fu1
      
      # Follow up 2
      df_fu2  <- df[which(df$timepoint==3), ] 
      df_fu2  <- subset(df_fu2,  subID %in% df_fu1_matched$subID)  # Get only the subjects that were in fu1 matched
      
      temp1 <- df_fu2[ , c("ID","age","scode","euler","dcode")]   # samples to match
      temp2 <- matchit(dcode ~ age+scode, data=temp1)      # matching
      temp3 <- match.data(temp2)[1:ncol(temp1)]            # matched sample
      
      df_fu2_matched  <- subset(df_fu2,  ID %in% temp3$ID) # extracting matched sample for fu2
      
      df <- rbind(df_bl_matched, df_fu1_matched, df_fu2_matched)
      
      write.csv(df, paste0("/data_J/Data/", paste0("parc", parc),"/",paste0("parc", parc),paste0("_", measure) , "_harmonized_matched.csv"), row.names = FALSE)
  }
  cat("Done!\n")
  return(df)
}

DataPreprocessing_males <- function(df, parc, measure, harmonization, match){
  # This function will prepare the df by harmonizing or matching.
  
  # HARMONIZATION WITH LONGITUDINAL COMBAT
  if (harmonization){ 
      cat("File with longCombat harmonization doesn't exist, creating...\n")
      measures <- c("CT", "Vol") 
      measure2 <- setdiff(measures, measure)

      feat_names <- colnames(df%>%select(ends_with(measure)))
      lC_measure <- longCombat(idvar    = 'subID', 
                               timevar  = 'timepoint',
                               batchvar = 'acode', 
                               features = feat_names, 
                               formula  = "age + dcode*timepoint",
                               ranef = "(1|subID)",
                               data  = df)
      
      data_lC_measure <- lC_measure$data_combat
      colnames(data_lC_measure)[4:ncol(data_lC_measure)] <- feat_names
      
      all_colnames <- colnames(df)
      data_lC_measure_colnames <- colnames(data_lC_measure)
      diff_colnames <- c(setdiff(all_colnames, data_lC_measure_colnames), 
                         setdiff(data_lC_measure_colnames, all_colnames))
      
      df_to_merge <- df[c(diff_colnames,"subID","timepoint")] # Subset df
      df <- merge(df_to_merge, data_lC_measure, by=c("subID","timepoint"))
      
      feat_names  <- c(colnames(df%>%dplyr::select(ends_with("CT"))), 
                       colnames(df%>%dplyr::select(ends_with("Vol"))))
      cov_names   <- setdiff(colnames(df), feat_names)
      feat_names <- colnames(df%>%select(ends_with(measure)))
      df <- df[,c(sort(cov_names), feat_names)]
      
      write.csv(df,paste0("/data_J/Data/", paste0("parc", parc),"/",paste0("parc", parc),paste0("_", measure),"_lC_harmonized_MALESONLY.csv"), row.names = FALSE)
}
  
  # MATCH-IT ALGORITHM
  if (match){
      cat("File with MATCH-IT doesn't exist, creating...\n")
      set.seed(99) # Otherwise every match-it would lead to different samples
      
      # Baseline
      df_bl  <- df[which(df$timepoint==1), ] 
      temp1  <- df_bl[ , c("ID","age","euler","dcode")] # sample to match
      temp2 <- matchit(dcode ~ age, data=temp1) # matching method by default: nearest neighbor
      # formula: two-sided formula() object containing the treatment and 
      # the covariates to be used in creating the distance measure used in the matching. 
      # This formula will be supplied to the functions that estimate the distance measure. 
      # The formula should be specified as A ~ X1 + X2 + ... where 
      # A represents the treatment variable and X1 and X2 are covariates.
      temp3 <- match.data(temp2)[1:ncol(temp1)] # matched sample
      
      df_bl_matched  <- subset(df_bl,  ID %in% temp3$ID) # extracting matched sample for baseline
      
      # Follow up 1
      df_fu1  <- df[which(df$timepoint==2), ]
      df_fu1  <- subset(df_fu1,  subID %in% df_bl_matched$subID)  # Get only the subjects that were in bl matched
      
      temp1 <- df_fu1[ , c("ID","age","euler","dcode")]   # samples to match
      temp2 <- matchit(dcode ~ age, data=temp1 )   # matching
      temp3 <- match.data(temp2)[1:ncol(temp1)]          # matched sample
      
      df_fu1_matched  <- subset(df_fu1,  ID %in% temp3$ID) # extracting matched sample for fu1
      
      # Follow up 2
      df_fu2  <- df[which(df$timepoint==3), ] 
      df_fu2  <- subset(df_fu2,  subID %in% df_fu1_matched$subID)  # Get only the subjects that were in fu1 matched
      
      temp1 <- df_fu2[ , c("ID","age","euler","dcode")]   # samples to match
      temp2 <- matchit(dcode ~ age, data=temp1)      # matching
      temp3 <- match.data(temp2)[1:ncol(temp1)]            # matched sample
      
      df_fu2_matched  <- subset(df_fu2,  ID %in% temp3$ID) # extracting matched sample for fu2
      
      df <- rbind(df_bl_matched, df_fu1_matched, df_fu2_matched)
      
      write.csv(df, paste0("/data_J/Data/", paste0("parc", parc),"/",paste0("parc", parc),paste0("_", measure) , "_harmonized_matched_MALESONLY.csv"), row.names = FALSE)
  }
  cat("Done!\n")
  return(df)
}