df <- df_lC_matched[1:40]
measure = "CT_freesurfer"
i=1

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
             data = data_hc, control=lmeControl(opt='optim', msMaxIter=150),
             method='REML')
model$sigma
  
methods(sigma)
  

```{r eval=FALSE, include=FALSE}
# PATIENTS
par(mfrow = c(1, 3))
for(tp in 1:3){
  Z <- Zs_match %>%
    filter(group=="patient") %>% 
    filter(timepoint==tp) %>% 
    mutate(deviant = ifelse(z > -1.96 & z < 1.96, 0, z)) %>%
    mutate(deviant = ifelse(z < -1.96, -1, deviant)) %>%
    mutate(deviant = ifelse(z >  1.96, 1 , deviant)) %>%
    filter(deviant==1)
  breakpoints <- seq(from = 0, to = 1000, by = 1)
  hist(Z$subID, breaks = breakpoints, main = paste0("Patients, timepoint ",tp), ylab = "Frecuencia", xlab = "patient")
}

# CONTROLS
par(mfrow = c(1, 3))
for(tp in 1:3){
  Z <- Zs_match %>%
    filter(group=="control") %>% 
    filter(timepoint==tp) %>% 
    mutate(deviant = ifelse(z > -1.96 & z < 1.96, 0, z)) %>%
    mutate(deviant = ifelse(z < -1.96, -1, deviant)) %>%
    mutate(deviant = ifelse(z >  1.96, 1 , deviant)) %>%
    filter(deviant==1)
  breakpoints <- seq(from = 0, to = 1000, by = 1)
  hist(Z$subID, breaks = breakpoints, main = paste0("Control, timepoint ",tp), ylab = "Frecuencia", xlab = "patient")
}


# controles y sujetos en la misma gráfica, diferenciar por
for (g in c("patient", "control")){
  for(tp in 1:3){
    data <- Zs_match %>%
      filter(group==g) %>% 
      mutate(deviant = ifelse(z > 1.96 | z < -1.96, 1, 0)) %>%
      filter(deviant!=0) %>% 
      mutate(timepoint = factor(timepoint, labels = c("bl", "fu1", "fu2")),
             subID = factor(subID))
    
    # Bar chart side by side
    plot <- ggplot(data, aes(x = subID, fill = timepoint)) +
      geom_bar(position = position_dodge()) + 
      
      labs(title=paste0(g,", timepoint ", tp), x= g , y = "n regions") + 
      scale_fill_manual(values=c("#999999", "#E69F00")) +
      theme_minimal() + 
      theme(axis.text.x = element_text(size = 5, angle = 90))   
    print(plot)
  }
}


data <- Zs_match   %>%
  mutate(deviant = ifelse(z > -1.96 & z < 1.96, 0, z)) %>%
  mutate(deviant = ifelse(z < -1.96, -1, deviant)) %>%
  mutate(deviant = ifelse(z >  1.96,  1, deviant)) %>%
  filter(deviant!=0) %>% 
  group_by(subID, group, deviant, timepoint) %>%
  dplyr::summarise(n=n()) %>%
  mutate(deviant = factor(deviant, labels = c("infra-dev","supra-dev")),
         subID = factor(subID),
         group = factor(group),
         timepoint = factor(timepoint, labels = c("Baseline","FU1","FU2")),
  )

# Bar chart side by side
plot <- ggplot(data, aes(x = subID, y = n, color = deviant)) +
  geom_bar(stat = "identity", fill="white") # , position = "fill" 
#labs(title=paste0(g,", timepoint ", tp), x="Dose (mg)", y = "Length") + 
#scale_fill_manual(values=c("#999999", "#E69F00")) +
#theme_minimal()
print(plot)
```
