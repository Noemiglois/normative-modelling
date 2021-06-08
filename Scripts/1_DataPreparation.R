DataPreparation <- function(parc, harmonization, match){
  
  if (parc == "parc35"){ 
    measure = "CT_freesurfer"
    
    if (!file.exists(paste0("/data_J/Data/parc35/",parc,".csv"))){
      cat("File parc35 doesn't exist, creating...\n") # Borrar
      df_lh   <- read.csv("/data_J/Data/lh.aparc.thickness.csv")[,1:35]
      df_rh   <- read.csv("/data_J/Data/rh.aparc.thickness.csv")[,1:35]
      df_parc <- merge(df_lh, df_rh, by="ID", all.x=T)
      
      df_morphosim <- read.csv("/data_J/Data/planU_morphosim_rawdatabase.csv")
      df_euler <- read.csv("/data_J/Data/planU_euler.csv")
      temp <- as.vector(df_euler[,3])
      df_euler[,3] <- impute(temp, fun=median)
      df_raw <- merge(df_morphosim, df_euler, by="ID", all.x=T)
      df_raw <- df_raw[, c("age", "scode", "dcode", "acode", "timepoint",
                           "subID", "euler_lh", "euler_rh","ID")]
      
      df <- merge(df_raw, df_parc, by="ID", all.x=T)
      
      # Columns creation: euler and dcode_age
      df <- transform(df, euler=(euler_lh+euler_rh)/2) # Create euler variable
      df <- transform(df, dcode_age=dcode*age)         # Create dcode_age variable
      df <- df[order(df$age), ]                        # Sort by age (not sure why)
      
      write.csv(df,paste0("/data_J/Data/parc35/",parc,".csv"), row.names = FALSE)
    }
    else if (file.exists(paste0("/data_J/Data/parc35/",parc,".csv"))){
      cat("File parc35 already exists, reading file...\n") # Borrar
      df <- read.csv(paste0("/data_J/Data/parc35/",parc,".csv"))
    }
  }
  
  else if (parc == "parc308"){
    measure = "thickness"
    if (!file.exists(paste0("/data_J/Data/parc308/",parc,".csv"))){ 
      cat("File parc308 doesn't exist, creating....\n") # Borrar
      df_morphosim <- read.csv("/data_J/Data/planU_morphosim_rawdatabase.csv")
      df_euler <- read.csv("/data_J/Data/planU_euler.csv")
      temp <- as.vector(df_euler[,3])
      df_euler[,3] <- impute(temp, fun=median)
      df_raw <- merge(df_morphosim, df_euler, by="ID", all.x=T)
      
      cov_names <- c("age", "scode", "dcode", "acode", "timepoint",
                     "subID", "euler_lh", "euler_rh", "ID")
      feat_names_thickness <- colnames(df_raw%>%select(ends_with("thickness")))
      
      df <- df_raw[, c(cov_names, feat_names_thickness)]
      
      # Columns creation: euler and dcode_age
      df <- transform(df, euler=(euler_lh+euler_rh)/2) # Create euler variable
      df <- transform(df, dcode_age=dcode*age)         # Create dcode_age variable
      df <- df[order(df$age), ]                        # Sort by age (not sure why)
      
      write.csv(df,"/data_J/Data/parc308/parc308.csv", row.names = FALSE)}
    else if (file.exists(paste0("/data_J/Data/parc308/",parc,".csv"))){ 
      cat("File parc308 already exists, reading file...\n") # Borrar
      df <- read.csv(paste0("/data_J/Data/parc308/",parc,".csv"))
      
      }
  }

  # Harmonization
  if (harmonization=="lC"){ 
    if (!file.exists(paste0("/data_J/Data/", parc,"/",parc,"_lC_harmonized.csv"))){
      cat("File with longCombat harmonization doesn't exist, creating...\n") # Borrar
      
      formula = 'age + scode + dcode*timepoint'
      ranef = '(1|subID)'
      feat_names <- colnames(df%>%select(ends_with(measure)))
      
      lC_measure <- longCombat(idvar = 'subID', 
                               timevar = 'timepoint',
                               batchvar = 'acode', 
                               features = feat_names, 
                               formula = formula,
                               ranef = ranef,
                               data = df)
      
      data_lC_measure <- lC_measure$data_combat
      colnames(data_lC_measure)[4:ncol(data_lC_measure)] <- feat_names
      
      all_colnames <- colnames(df)
      data_lC_measure_colnames <- colnames(data_lC_measure)
      diff_colnames <- c(setdiff(all_colnames, data_lC_measure_colnames), 
                         setdiff(data_lC_measure_colnames, all_colnames))
      
      df_to_merge <- df[c(diff_colnames,"subID","timepoint")] # Subset df
      df <- merge(df_to_merge, data_lC_measure, by=c("subID","timepoint"))

      write.csv(df,paste0("/data_J/Data/", parc,"/",parc,"_lC_harmonized.csv"), row.names = FALSE)
      }
    else if (file.exists(paste0("/data_J/Data/", parc,"/",parc,"_lC_harmonized.csv"))){
      cat("File with longCombat harmonization already exists, reading file...\n") # Borrar
      df <- read.csv(paste0("/data_J/Data/", parc,"/",parc,"_lC_harmonized.csv"))
    }
  } 
  
  else if (harmonization=="nC"){ 
    if (!file.exists(paste0("/data_J/Data/", parc,"/",parc,"_nC_harmonized.csv"))){
      cat("File with neuroCombat harmonization doesn't exist, creating...\n") # Borrar
      
      batch <- as.numeric(df$acode)
      age <- as.numeric(df$age)   
      scode <- as.factor(df$scode) # Categorical variable
      dcode <- as.factor(df$dcode) # Categorical variable
      timepoint <- as.factor(df$timepoint) # Categorical variable
      
      feat_names <- colnames(df%>%select(ends_with(measure)))
      
      mod <- model.matrix(~age + scode + dcode + timepoint)
      
      df_measure <- df[feat_names] # Subset df
      
      dat <- t(as.matrix(df_measure))
      nC_measure <- neuroCombat(dat=dat, batch=batch, mod=mod)
      data_nC_measure <-  as.data.frame(t(nC_measure$dat.combat))
      
      data_nC_measure$subID = df$subID
      data_nC_measure$timepoint = df$timepoint
      
      data_nC_measure_colnames <- colnames(data_nC_measure)
      all_colnames <- colnames(df)
      diff_colnames<-c(setdiff(all_colnames, data_nC_measure_colnames), 
                       setdiff(data_nC_measure_colnames, all_colnames))
      df_to_merge <- df[c(diff_colnames, "subID", "timepoint")]
      df <- merge(df_to_merge, data_nC_measure, by=c("subID","timepoint"))
      
      write.csv(df,paste0("/data_J/Data/", parc,"/",parc,"_nC_harmonized.csv"), row.names = FALSE)
    }
    else if (file.exists(paste0("/data_J/Data/", parc,"/",parc,"_nC_harmonized.csv"))){
      df <- read.csv(paste0("/data_J/Data/", parc,"/",parc,"_nC_harmonized.csv"))
      cat("File with neuroCombat harmonization already exists, reading file...\n") # Borrar
    }
  }

  # Matching
  if (match){
    if (!file.exists(paste0("/data_J/Data/", parc,"/",parc,"_",harmonization,"_harmonized_MATCHED.csv"))){
      cat("File with MATCH-IT doesn't exist, creating...\n") # Borrar

      set.seed(99) # Every matchit would lead to different samples
      
      # Baseline
      df_bl  <- df[which(df$timepoint==1), ] 
      
      temp1  <- df_bl[ , c("ID","age","scode","euler","dcode")] # samples to match
      temp2 <- matchit(dcode ~ age+scode, data=temp1)    # matching
      temp3 <- match.data(temp2)[1:ncol(temp1)]          # matched sample

      df_bl_matched  <- subset(df_bl,  ID %in% temp3$ID) # extracting matched sample for baseline
      
      # Follow up 1
      df_fu1  <- df[which(df$timepoint==2), ]
      df_fu1  <- subset(df_fu1,  subID %in% df_bl_matched$subID)  # Get only the subjects that were in bl matched
      
      temp1 <- df_fu1[ , c("ID","age","scode","euler","dcode")]   # samples to match
      temp2 <- matchit(dcode ~ age+scode, data=temp1 )   # matching
      temp3 <- match.data(temp2)[1:ncol(temp1)]          # matched sample
      
      df_fu1_matched  <- subset(df_fu1,  ID %in% temp3$ID) # extracting matched sample for baseline
      
      # Follow up 2
      df_fu2  <- df[which(df$timepoint==3), ] 
      df_fu2  <- subset(df_fu2,  subID %in% df_fu1_matched$subID)  # Get only the subjects that were in fu1 matched
      
      temp1 <- df_fu2[ , c("ID","age","scode","euler","dcode")]   # samples to match
      temp2 <- matchit(dcode ~ age+scode, data=temp1)      # matching
      temp3 <- match.data(temp2)[1:ncol(temp1)]            # matched sample
      
      df_fu2_matched  <- subset(df_fu2,  ID %in% temp3$ID) # extracting matched sample for baseline
      
      df <- rbind(df_bl_matched, df_fu1_matched, df_fu2_matched)
      
      write.csv(df,paste0("/data_J/Data/", parc,"/",parc,"_",harmonization,"_harmonized_MATCHED.csv"), row.names = FALSE)
    }
    else if (file.exists(paste0("/data_J/Data/", parc,"/",parc,"_",harmonization,"_harmonized_MATCHED.csv"))){
      df <- read.csv(paste0("/data_J/Data/", parc,"/",parc,"_",harmonization,"_harmonized_MATCHED.csv"))
      cat("File with MATCH-IT already exists, reading file...\n") # Borrar
    }
  }
  cat("Done!\n")
  return(df)
}