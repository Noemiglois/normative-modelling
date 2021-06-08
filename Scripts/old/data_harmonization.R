# PACKAGES ----------------------------------------------------------------
library(devtools)
library('longCombat')
library('neuroCombat')
library(dplyr)
library(Hmisc)
# devtools::install_version("latticeExtra", version="0.6-28")
# install.packages("Hmisc")
# sudo apt install r-cran-Hmisc

# EXPLORATORY DATA ANALYSIS --------------------------------------------------------------------
df_morphosim <- read.csv("/home/jjanssen/noemi_r/morphosim/data/planU_morphosim_rawdatabase.csv")
df_euler <- read.csv("/home/jjanssen/noemi_r/morphosim/data/planU_euler.csv")

"
Euler dataframe has 1 NaN value. As it is only one, and we will use this 
covariate to perform the harmonization, we decided to impute it by the median.
"
x<-as.vector(df_euler[,3])
df_euler[,3]<-impute(x,fun=median)

df <- merge(df_morphosim, df_euler, by="ID", all.x=T) # by default = merge(df_morphosim, df_euler)

"
We check whether each subject has the different images taken in the same scanner.
To do so, we count the number of unique IDs and compare that amount with the 
number of 0 in table (subID, acode), which adds up how many images were taken in
each scanner for each subject.
"
dim(df[4]) # 1087 images
n_sub <- dim(unique(df[4]))[1] # 467 different subjects
sum(table(df$subID, df$acode)==0)== n_sub 

"Now we start by building the subsets for each feature, cointaining 308 measures 
each, plus the 8 covariates. Also, we check whether there are any NA values.
"

cov_names <- c("age", "scode", "dcode", "acode", "timepoint", "subID", "euler_lh", "euler_rh")

# Area
feat_names_area <- colnames(df%>%select(ends_with("area")))
df_area <- df[c(cov_names, feat_names_area)] # Subset df
any(is.na(df_area)) # Check NA

# Curvind
feat_names_curvind <- colnames(df%>%select(ends_with("curvind")))
df_curvind <- df[c(cov_names, feat_names_curvind)] # Subset df
any(is.na(df_curvind)) # Check NA

# Foldind
feat_names_foldind <- colnames(df%>%select(ends_with("foldind")))
df_foldind <- df[c(cov_names, feat_names_foldind)] # Subset df
any(is.na(df_foldind)) # Check NA

# Gauscurv
feat_names_gauscurv <- colnames(df%>%select(ends_with("gauscurv")))
df_gauscurv <- df[c(cov_names, feat_names_gauscurv)] # Subset df
any(is.na(df_gauscurv)) # Check NA

# Meancurv
feat_names_meancurv <- colnames(df%>%select(ends_with("meancurv")))
df_meancurv <- df[c(cov_names, feat_names_meancurv)] # Subset df
any(is.na(df_meancurv)) # Check NA

# Thickness
feat_names_thickness <- colnames(df%>%select(ends_with("thickness")))
df_thickness <- df[c(cov_names, feat_names_thickness)] # Subset df
any(is.na(df_thickness)) # Check NA
all_colnames<-c(colnames(df[1:97]),colnames(df[feat_names_thickness]))

# Volume
feat_names_volume <- colnames(df%>%select(ends_with("volume")))
df_volume <- df[c(cov_names, feat_names_volume)] # Subset df
any(is.na(df_volume)) # Check NA

# -------------- LONG COMBAT HARMONIZATION: -----------------------------------

"
# 467 subjects with 1, 2 or 3 time points each
# 2 scanners / batches
# 8 covariates + 7x308 features
"

formula='age + scode + dcode*timepoint'
ranef='(1|subID)'

# 1. AREA -----------------------------------------------------------------
lC_area <- longCombat(idvar='subID', 
                      timevar='timepoint',
                      batchvar='acode', 
                      features=feat_names_area, 
                      formula=formula,
                      ranef=ranef,
                      data=df_area)

data_lC_area <- lC_area$data_combat
colnames(data_lC_area)[4:311] <- feat_names_area

# 2. CURVIND --------------------------------------------------------------
lC_curvind <- longCombat(idvar='subID', 
                         timevar='timepoint',
                         batchvar='acode', 
                         features=feat_names_curvind, 
                         formula=formula,
                         ranef=ranef,
                         data=df_curvind)

data_lC_curvind <- lC_curvind$data_combat
colnames(data_lC_curvind)[4:311] <- feat_names_curvind

# 3. FOLDIND --------------------------------------------------------------
lC_foldind <- longCombat(idvar='subID', 
                         timevar='timepoint',
                         batchvar='acode', 
                         features=feat_names_foldind, 
                         formula=formula,
                         ranef=ranef,
                         data=df_foldind)

data_lC_foldind <- lC_foldind$data_combat
colnames(data_lC_foldind)[4:311] <- feat_names_foldind

# 4. GAUSCURV ---------------------------------------------------------------
lC_gauscurv <- longCombat(idvar='subID', 
                          timevar='timepoint',
                          batchvar='acode', 
                          features=feat_names_gauscurv, 
                          formula=formula,
                          ranef=ranef,
                          data=df_gauscurv)

data_lC_gauscurv <- lC_gauscurv$data_combat
colnames(data_lC_gauscurv)[4:311] <- feat_names_gauscurv

# 5. MEANCURV -------------------------------------------------------------
lC_meancurv <- longCombat(idvar='subID', 
                          timevar='timepoint',
                          batchvar='acode', 
                          features=feat_names_meancurv, 
                          formula=formula,
                          ranef=ranef,
                          data=df_meancurv)

data_lC_meancurv <- lC_meancurv$data_combat
colnames(data_lC_meancurv)[4:311] <- feat_names_meancurv

# 6. THICKNESS ------------------------------------------------------------
lC_thickness <- longCombat(idvar='subID', 
                           timevar='timepoint',
                           batchvar='acode', 
                           features=feat_names_thickness, 
                           formula=formula,
                           ranef=ranef,
                           data=df_thickness)

data_lC_thickness <- lC_thickness$data_combat
colnames(data_lC_thickness)[4:311] <- feat_names_thickness

# 7. VOLUME --------------------------------------------------------------
lC_volume <- longCombat(idvar='subID', 
                        timevar='timepoint',
                        batchvar='acode', 
                        features=feat_names_volume, 
                        formula=formula,
                        ranef=ranef,
                        data=df_volume)

data_lC_volume <- lC_volume$data_combat
colnames(data_lC_volume)[4:311] <- feat_names_volume

# Save thickness lC harmonized dataset-----
data_lC_thickness_colnames<-colnames(data_lC_thickness)
diff_colnames<-c(setdiff(all_colnames, data_lC_thickness_colnames), setdiff(data_lC_thickness_colnames, all_colnames))
df_to_merge <- df[c(diff_colnames,"subID","timepoint")] # Subset df
thickness_lC<-merge(df_to_merge, data_lC_thickness, by=c("subID","timepoint"))

write.csv(thickness_lC,"/home/jjanssen/noemi_r/morphosim/harmonized_data/longCombat/morphosim_lC_harmonized_thickness.csv", row.names = FALSE)

# Save whole lC harmonized dataset-----
df_lC1 <- merge(data_lC_area, data_lC_curvind)
df_lC2 <- merge(df_lC1, data_lC_foldind)
df_lC3 <- merge(df_lC2, data_lC_gauscurv)
df_lC4 <- merge(df_lC3, data_lC_meancurv)
df_lC5 <- merge(df_lC4, data_lC_thickness)
df_lC <- merge(df_lC5, data_lC_volume)
rm(df_lC1,df_lC2,df_lC3,df_lC4,df_lC5)

df_lC_colnames<-colnames(df_lC)
diff_colnames<-c(setdiff(colnames(df), df_lC_colnames), setdiff(df_lC_colnames, colnames(df)))
df_to_merge <- df[c(diff_colnames,"subID","timepoint")] # Subset df
df_lC<-merge(df_to_merge, df_lC, by=c("subID","timepoint"))

write.csv(df_lC,"/home/jjanssen/noemi_r/morphosim/harmonized_data/longCombat/morphosim_lC_harmonized_data.csv", row.names = FALSE)


# -------------- NEURO COMBAT HARMONIZATION: -----------------------------------

"
# Columns (n) are participants (n=1.087)
# Rows (p) are features (p=314)
# p<<n --> without empirical Bayes
# Data matrix should be (pxn), so we have to transpose the previous ones.
# Batch is a numeric vector of length n indicating the scanner
"
batch <- as.numeric(df$acode)

age <- as.numeric(df$age) # Continuous variable
scode <- as.factor(df$scode) # Categorical variable
dcode <- as.factor(df$dcode) # Categorical variable
timepoint <- as.factor(df$timepoint) # Categorical variable

mod <- model.matrix(~age + scode + dcode + timepoint)

# 1. AREA -----------------------------------------------------------------
feat_names_area <- colnames(df%>%select(ends_with("area")))
df_area <- df[feat_names_area] # Subset df
any(is.na(df_area)) # Check NA

dat <- t(as.matrix(df_area))
nC_area <- neuroCombat(dat=dat, batch=batch, mod=mod)
data_nC_area <-  as.data.frame(t(nC_area$dat.combat))

# 2. CURVIND -----------------------------------------------------------------
feat_names_curvind <- colnames(df%>%select(ends_with("curvind")))
df_curvind <- df[feat_names_curvind] # Subset df
any(is.na(df_curvind)) # Check NA

dat <- t(as.matrix(df_curvind))
nC_curvind <- neuroCombat(dat=dat, batch=batch, mod=mod)
data_nC_curvind <-  as.data.frame(t(nC_curvind$dat.combat))

# 3. FOLDIND -----------------------------------------------------------------
feat_names_foldind <- colnames(df%>%select(ends_with("foldind")))
df_foldind <- df[feat_names_foldind] # Subset df
any(is.na(df_foldind)) # Check NA

dat <- t(as.matrix(df_foldind))
nC_foldind <- neuroCombat(dat=dat, batch=batch, mod=mod)
data_nC_foldind <-  as.data.frame(t(nC_foldind$dat.combat))

# 4. GAUSCURV -----------------------------------------------------------------
feat_names_gauscurv <- colnames(df%>%select(ends_with("gauscurv")))
df_gauscurv <- df[feat_names_gauscurv] # Subset df
any(is.na(df_gauscurv)) # Check NA

dat <- t(as.matrix(df_gauscurv))
nC_gauscurv <- neuroCombat(dat=dat, batch=batch, mod=mod)
data_nC_gauscurv <-  as.data.frame(t(nC_gauscurv$dat.combat))

# 5. MEANCURV -----------------------------------------------------------------
feat_names_meancurv <- colnames(df%>%select(ends_with("meancurv")))
df_meancurv <- df[feat_names_meancurv] # Subset df
any(is.na(df_meancurv)) # Check NA

dat <- t(as.matrix(df_meancurv))
nC_meancurv <- neuroCombat(dat=dat, batch=batch, mod=mod)
data_nC_meancurv <-  as.data.frame(t(nC_meancurv$dat.combat))

# 6. THICKNESS -----------------------------------------------------------------
feat_names_thickness <- colnames(df%>%select(ends_with("thickness")))
df_thickness <- df[feat_names_thickness] # Subset df
any(is.na(df_thickness)) # Check NA

dat <- t(as.matrix(df_thickness))
nC_thickness <- neuroCombat(dat=dat, batch=batch, mod=mod)
data_nC_thickness <-  as.data.frame(t(nC_thickness$dat.combat))

# 7. VOLUME -----------------------------------------------------------------
feat_names_volume <- colnames(df%>%select(ends_with("volume")))
df_volume <- df[feat_names_volume] # Subset df
any(is.na(df_volume)) # Check NA

dat <- t(as.matrix(df_volume))
nC_volume <- neuroCombat(dat=dat, batch=batch, mod=mod)
data_nC_volume <-  as.data.frame(t(nC_volume$dat.combat))

# Save thickness nC harmonized dataset -----
data_nC_thickness$subID = df$subID
data_nC_thickness$timepoint = df$timepoint

data_nC_thickness_colnames<-colnames(data_nC_thickness)
all_colnames<-c(colnames(df[1:97]), colnames(df[feat_names_thickness]))
diff_colnames<-c(setdiff(all_colnames, data_nC_thickness_colnames), setdiff(data_nC_thickness_colnames, all_colnames))
df_to_merge <- df[c(diff_colnames, "subID", "timepoint")] # Subset df
thickness_nC<-merge(df_to_merge, data_nC_thickness, by=c("subID","timepoint"))

write.csv(thickness_nC,"/home/jjanssen/noemi_r/morphosim/harmonized_data/neuroCombat/morphosim_nC_harmonized_thickness.csv", row.names = FALSE)

# Save whole nC harmonized dataset -----
data_nC_area$subID = df$subID
data_nC_area$timepoint = df$timepoint

data_nC_curvind$subID = df$subID
data_nC_curvind$timepoint = df$timepoint

data_nC_foldind$subID = df$subID
data_nC_foldind$timepoint = df$timepoint

data_nC_gauscurv$subID = df$subID
data_nC_gauscurv$timepoint = df$timepoint

data_nC_meancurv$subID = df$subID
data_nC_meancurv$timepoint = df$timepoint

data_nC_volume$subID = df$subID
data_nC_volume$timepoint = df$timepoint

df_nC1 <- merge(data_nC_area, data_nC_curvind)
df_nC2 <- merge(df_nC1, data_nC_foldind)
df_nC3 <- merge(df_nC2, data_nC_gauscurv)
df_nC4 <- merge(df_nC3, data_nC_meancurv)
df_nC5 <- merge(df_nC4, data_nC_thickness)
df_nC <- merge(df_nC5, data_nC_volume)
rm(df_nC1, df_nC2, df_nC3, df_nC4, df_nC5)

df_nC_colnames<-colnames(df_nC)
diff_colnames<-c(setdiff(colnames(df), df_nC_colnames), setdiff(df_nC_colnames, colnames(df)))
df_to_merge <- df[c(diff_colnames,"subID","timepoint")] # Subset df
df_nC<-merge(df_to_merge, df_nC, by=c("subID","timepoint"))

write.csv(df_nC,"/home/jjanssen/noemi_r/morphosim/harmonized_data/neuroCombat/morphosim_nC_harmonized_data.csv", row.names = FALSE)

