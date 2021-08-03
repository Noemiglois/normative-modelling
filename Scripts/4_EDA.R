EDA_match_ages <- function(df, match){
  # df <- df_lC_NO_matched
  idx_pat    <- which(df$dcode== 1)  # Index of 1: patients
  idx_hc     <- which(df$dcode== 0)  # Index of 0: healthy controls
  data_pat   <- df[idx_pat, ]        # Select rows corresponding to patients
  data_hc    <- df[idx_hc,  ]        # Select rows corresponding to healthy controls
  
  age_pat <- sort(unique(floor(data_pat$age)), decreasing = T) # Floor = truncate
  age_hc  <- sort(unique(floor(data_hc$age)) , decreasing = T) # and sort by decreasing order

  # setdiff(age_hc, age_pat) # Ages that are in controls but not in patients
  # setdiff(age_pat,age_hc)  # Ages that are in patients but not in controls
  ages_no_match <- sort(c(setdiff(age_hc, age_pat),setdiff(age_pat,age_hc))) # Ages that doesn't match
  return(cat(match, "dataset: \n Ages that doesn't match in patients vs controls: ",
             ages_no_match, "\n Min age: ", min(df$age), "\n Max age: ", max(df$age)))
}


variancePartition_lm <- function(df, measure, lab, par){
  ## PREPARE DATA ##
  
  reg_names <- colnames(df%>%dplyr::select(ends_with(measure)))
  
  # Rename predictors
  names(df)[names(df) == "age"]   <- "Age"
  names(df)[names(df) == "scode"] <- "Sex"
  names(df)[names(df) == "euler"] <- "Euler_number"
  names(df)[names(df) == "dcode"] <- "Diagnosis"
  names(df)[names(df) == "subID"] <- "Individual"
  
  # Type coercion
  df$Sex <- as.factor(df$Sex)
  df$Individual <- as.factor(df$Individual)
  df$Diagnosis  <- as.factor(df$Diagnosis)
  
  ## VIOLIN PLOTS ##
  
  # Timepoints
  n_tp <- c(1, 2, 3)
  for (tp in n_tp){
    # Timepoint selection
    df_tp <- df %>%
      filter(timepoint==tp)
    
    # Select regions, transpose
    df_base_ima <- as.data.frame(t(df_tp[,reg_names]))
    
    # Create vectors for each predictor
    Age <- df_tp$Age
    Sex <- df_tp$Sex
    Diagnosis    <- df_tp$Diagnosis
    Euler_number <- df_tp$Euler_number
    Individual   <- df_tp$Individual
    
    # Create formula
    # form <- ~Age+(1|Sex)+(1|Diagnosis)+Euler_numer+(1|Individual)
    form <- ~ Age+Sex+Diagnosis+Euler_number # General/classic linear model
    
    # Create dataframe with predictors
    info <- df_tp[,c("Age", "Sex", "Diagnosis","Euler_number")]
    
    # Run variance partition
    varPart <- fitExtractVarPartModel(df_base_ima, form, info)
    
    # Sort columns
    vp <- sortCols(varPart)
    
    # Plot sorted columns as violins, setting colors
    #theme_set(theme_gray(base_size = 40))
    violins <- plotVarPart(vp, col = c("#5DC863FF", "#3B528BFF", "#21908CFF",
                                      "#440154FF", "lightgrey"))
    # Save as png
    png(file=paste0("/data_J/Results/Plots/",par,"/violin_tp",tp,"_",lab,".png"),
        width=5, height=2.5, units="in",res=600)
    print(violins)
    dev.off()
  }
}


variancePartition_lme <- function(df, measure, lab, par){
  ## PREPARE DATA ##
  reg_names <- colnames(df%>%dplyr::select(ends_with(measure)))
  
  # Rename predictors
  names(df)[names(df) == "age"]   <- "Age"
  names(df)[names(df) == "scode"] <- "Sex"
  names(df)[names(df) == "euler"] <- "Euler_number"
  names(df)[names(df) == "dcode"] <- "Diagnosis"
  names(df)[names(df) == "subID"] <- "Individual"
  
  # Type coercion
  df$Sex        <- as.factor(df$Sex)
  df$Individual <- as.factor(df$Individual)
  df$Diagnosis  <- as.factor(df$Diagnosis)
  
  ## VIOLIN PLOTS ##
  # Select regions, transpose
  df_base_ima <- as.data.frame(t(df[,reg_names]))
  
  # Create vectors for each predictor
  Age <- df$Age
  Sex <- df$Sex
  Diagnosis    <- df$Diagnosis
  Euler_number <- df$Euler_number
  Individual   <- df$Individual
  
  # Create formula
  # form <- ~ Age+Sex+Diagnosis+Euler_number # General/classic linear model
  # REGRESSION: FIXED: ~ scode + age + euler, RANDOM = ~1 + age|subID
  form <- ~ Age + (1|Sex) + (1|Diagnosis) + Euler_number + (1|Individual) # Linear mixed effects model
  
  # Create dataframe with predictors
  info <- df[,c("Age", "Sex", "Diagnosis", "Euler_number", "Individual")]
  
  # Run variance partition
  varPart <- fitExtractVarPartModel(df_base_ima, form, info, BPPARAM=SnowParam(25))
  
  # Sort columns
  vp <- sortCols(varPart)
  
  # Plot sorted columns as violins, setting colors
  theme_set(theme_gray(base_size = 40))
  violins <- plotVarPart(vp, col = c("#440154FF", "#5DC863FF", "#3B528BFF",
                                     "#21908CFF", "#440154FF", "lightgrey"))

  # Save as png
  png(file=paste0("/data_J/Results/Plots/", par, "/violin_", lab, ".png"),
      width=4, height=2.5, units="in", res=600)
  
  print(violins)
  
  dev.off()
}


pdf_plot <- function(Z, lab, par){
  # Histogram overlaid with kernel density curve
  for (g in unique(Z$group)){
    Zg <- subset(Z, group==g)
    pdf <- ggplot(Zg, aes(x=z)) +
      geom_histogram(aes(y=..density..),        # Histogram with density
                     binwidth=.3, colour="black", fill="white") +
      geom_density(alpha=.2, fill="#FF6666")  + # Transparent density plot
      geom_vline(aes(xintercept=mean(z, na.rm=T)),   # Ignore NA values for mean
                 color="red", linetype="dashed", size=0.5) +
      labs(title = paste0(g),
           x = "z-score",
           y = "density") +
      theme(plot.title = element_text(size=22)) +
      xlim(-5, 5)
    
    
    # Save as png
    png(file=paste0("/data_J/Results/Plots/",par,"/pdf_",g,"_",lab,".png"),
        width=10, height=5, units="in", res=300)
    print(pdf)
    dev.off()
  }
}

mean_Z_across_regions_plot <- function(Z, par, lab){
  for (g in unique(Z$group)){
    Zg <- subset(Z, group==g)
    z_mean_reg_match <- Zg %>%
      group_by(region) %>%
      dplyr::summarise(media = mean(z))
    x <- 1:length(z_mean_reg_match$media)
    means_plot <- plot(x, z_mean_reg_match$media, type = "l", lwd=1, 
                       xlab="region", ylab="mean z-score", 
                       main="Mean z-score per region (NO MATCH)")
    
    # Save as png
    png(file=paste0("/data_J/Results/Plots/",par,"/mean_Z_",g,"_",lab,".png"), 
        width=5, height=10, units="in", res=300)
    print(means_plot)
    dev.off()
  }
}


Z_across_regions_plot <- function(Z, par, lab){
  fig <- ggplot(Z, aes(y=z, x=group)) +
    geom_boxplot()
  # Save as png
  png(file=paste0("/data_J/Results/Plots/", par,"/boxplot_Z_", lab,".png"), width=10,
      height=5, units="in", res=300)
  print(fig)
  dev.off()
}




