# library(MatchIt)

df <- read.csv("/data_J/Data/lC_harmonized_thickness_reduced.csv")

df_bl  <- df[which(df$timepoint==1), ] # Baseline
df_fu1 <- df[which(df$timepoint==2), ] # Follow up 1
df_fu2 <- df[which(df$timepoint==3), ] # Follow up 2

temp1  <- df_bl[ ,c(4,1,2,5,7)]

set.seed(99) # Every matchit would lead to different samples
temp2 <- matchit(dcode ~ age+scode, data=temp1)    # matching
temp3 <- match.data(temp2)[1:ncol(temp1)]          # matched sample

df_bl_matched  <- subset(df_bl,  ID %in% temp3$ID)   # extracting matched sample for baseline
subIDs <- df_bl_matched$subID
rm(temp1, temp2, temp3)

df_fu1_matched <- subset(df_fu1, subID %in% df_bl_matched$subID)   # extracting matched sample for follow up 1
df_fu2_matched <- subset(df_fu2, subID %in% df_bl_matched$subID)   # extracting matched sample for follow up 2
rm(df_bl, df_fu1, df_fu2)

# table(df_bl_matched$dcode, df_bl_matched$scode)

df_matched <- rbind(df_bl_matched, df_fu1_matched, df_fu2_matched)
rm(df_bl_matched,df_fu1_matched,df_fu2_matched)

write.csv(df_matched,"/data_J/data/lC_harmonized_thickness_matched.csv", row.names = FALSE)

