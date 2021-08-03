list.of.packages <- c("viridis", "tidyverse", "MatchIt", "grid", "png", 
                      "gridExtra", "parallel", "nlme", "JMbayes", "plyr",
                      "BiocManager", "Biostrings", "lme4", "ggplot2", "dplyr",
                      "Hmisc", "devtools", "longCombat", "neuroCombat", 
                      "tinytex", "knitr", "ggExtra", "variancePartition",
                      "Rmisc", "doParallel", "hrbrthemes", "cowplot")

ip = as.data.frame(installed.packages()[,c(1,3:4)])
ip = ip[is.na(ip$Priority),1:2,drop=FALSE]
list.of.versions <- ip[list.of.packages,"Version"]


library(versions) # install.packages("versions")
install.versions(list.of.packages, list.of.versions)