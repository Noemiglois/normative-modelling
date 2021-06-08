# PACKAGES ----------------------------------------------------------------
library(devtools)
library('longCombat')
library('neuroCombat')
library(dplyr)

# EXPLORATORY DATA ANALYSIS --------------------------------------------------------------------
df <- read.csv("/Users/noemi/OneDrive/Escritorio/planU.peps.longit.aparc.thickness.NOcombat.csv")
all_colnames<-colnames(df)
cov_names <- c("age", "scode", "dcode", "acode", "timepoint", "subID")

# Freesurfer
feat_names_freesurfer <- colnames(df%>%select(ends_with("_freesurfer")))
df_freesurfer <- df[c(cov_names, feat_names_freesurfer)] # Subset df

# -------------- LONG COMBAT HARMONIZATION: -----------------------------------
formula='age + scode + dcode*timepoint'
ranef='(1|subID)'

lC_freesurfer <- longCombat(idvar='subID', 
                      timevar='timepoint',
                      batchvar='acode', 
                      features=feat_names_freesurfer, 
                      formula=formula,
                      ranef=ranef,
                      data=df_freesurfer)

data_lC_freesurfer <- lC_freesurfer$data_combat
colnames(data_lC_freesurfer)[4:73] <- feat_names_freesurfer

data_lC_freesurfer_colnames<-colnames(data_lC_freesurfer)
diff_colnames<-c(setdiff(all_colnames, data_lC_freesurfer_colnames), setdiff(data_lC_freesurfer_colnames, all_colnames))
df_to_merge <- df[c(diff_colnames,"subID","timepoint")] # Subset df
df_lC<-merge(df_to_merge, data_lC_freesurfer, by=c("subID","timepoint"))
write.csv(df_lC,"/Users/noemi/OneDrive/Escritorio/planU.peps.longit.aparc.thickness.NOcombat.LongCombat.csv", row.names = FALSE)

# -------------- NEURO COMBAT HARMONIZATION: -----------------------------------

dat <- t(as.matrix(df[feat_names_freesurfer]))

batch <- as.numeric(df$acode)

age <- as.numeric(df$age) # Continuous variable
scode <- as.factor(df$scode) # Categorical variable
dcode <- as.factor(df$dcode) # Categorical variable
timepoint <- as.factor(df$timepoint) # Categorical variable

mod <- model.matrix(~age + scode + dcode + timepoint)

combat.harmonized <- neuroCombat(dat=dat, batch=batch, mod=mod)

data_nC_freesurfer <-  as.data.frame(t(combat.harmonized$dat.combat))

data_nC_freesurfer$subID=df$subID
data_nC_freesurfer$timepoint=df$timepoint

data_nC_freesurfer_colnames<-colnames(data_nC_freesurfer)

cols<-c(setdiff(all_colnames, data_nC_freesurfer_colnames), setdiff(data_nC_freesurfer_colnames, all_colnames))
df_prev <- df[c(cols,"subID","timepoint")] # Subset df
df_nC<-merge(df_prev, data_nC_freesurfer, by=c("subID","timepoint"))

write.csv(df_nC,"/Users/noemi/OneDrive/Escritorio/planU.peps.longit.aparc.thickness.NOcombat.NeuroCombat.csv", row.names = FALSE)


