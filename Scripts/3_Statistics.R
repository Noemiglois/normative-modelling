getStatistics <- function(Zs, lab){
  pcol <- viridis(1)
  n_tp <- unique(Zs$timepoint)

  for (tp in n_tp){
    dataf <- Zs %>% 
      filter(timepoint==tp)%>%
      group_by(region, timepoint)%>%
      summarise(num = sum(abs(z)>1.96), den = sum(abs(z)<1.96)) %>%
      summarise(region=region, ratio = num/den)
    
    prev <- ggplot(dataf, aes(ratio, fill = pcol)) + 
      geom_density(alpha = 0.8) + theme_minimal() +
      xlab("Percentage of subjects with atypical w-score") + 
      ylab("number of \n brain regions") +
      ggtitle(paste0("Median prevalence: ", round(median(dataf$ratio),3), ' (Timepoint ',tp,')')) +
      theme(legend.position = 'none')
    
    # Save as png
    png(file=paste0("/data_J/Results/Plots/prevalence_wscores_",tp,"_", lab,".png"), width=10, height=2.5, units="in",res=600)
    print(prev)
    dev.off()
    
    # Create another function to call when we want to order and save
    # in txt to use in docker visualizations
    ####### Region names ordered ################
    reg_order <- c("lh_bankssts_part1_thickness","lh_bankssts_part2_thickness","lh_caudalanteriorcingulate_part1_thickness","lh_caudalmiddlefrontal_part1_thickness","lh_caudalmiddlefrontal_part2_thickness","lh_caudalmiddlefrontal_part3_thickness","lh_caudalmiddlefrontal_part4_thickness","lh_cuneus_part1_thickness","lh_cuneus_part2_thickness","lh_entorhinal_part1_thickness","lh_fusiform_part1_thickness",
                   "lh_fusiform_part2_thickness","lh_fusiform_part3_thickness","lh_fusiform_part4_thickness","lh_fusiform_part5_thickness","lh_inferiorparietal_part1_thickness","lh_inferiorparietal_part2_thickness","lh_inferiorparietal_part3_thickness","lh_inferiorparietal_part4_thickness","lh_inferiorparietal_part5_thickness","lh_inferiorparietal_part6_thickness","lh_inferiorparietal_part7_thickness",
                   "lh_inferiorparietal_part8_thickness","lh_inferiortemporal_part1_thickness","lh_inferiortemporal_part2_thickness","lh_inferiortemporal_part3_thickness","lh_inferiortemporal_part4_thickness","lh_inferiortemporal_part5_thickness","lh_inferiortemporal_part6_thickness","lh_isthmuscingulate_part1_thickness","lh_isthmuscingulate_part2_thickness","lh_lateraloccipital_part1_thickness",
                   "lh_lateraloccipital_part2_thickness","lh_lateraloccipital_part3_thickness","lh_lateraloccipital_part4_thickness","lh_lateraloccipital_part5_thickness","lh_lateraloccipital_part6_thickness","lh_lateraloccipital_part7_thickness","lh_lateraloccipital_part8_thickness","lh_lateraloccipital_part9_thickness","lh_lateralorbitofrontal_part1_thickness","lh_lateralorbitofrontal_part2_thickness",
                   "lh_lateralorbitofrontal_part3_thickness","lh_lateralorbitofrontal_part4_thickness","lh_lingual_part1_thickness","lh_lingual_part2_thickness","lh_lingual_part3_thickness","lh_lingual_part4_thickness","lh_lingual_part5_thickness","lh_lingual_part6_thickness","lh_medialorbitofrontal_part1_thickness","lh_medialorbitofrontal_part2_thickness","lh_medialorbitofrontal_part3_thickness",
                   "lh_middletemporal_part1_thickness","lh_middletemporal_part2_thickness","lh_middletemporal_part3_thickness","lh_middletemporal_part4_thickness","lh_middletemporal_part5_thickness","lh_parahippocampal_part1_thickness","lh_parahippocampal_part2_thickness","lh_paracentral_part1_thickness","lh_paracentral_part2_thickness","lh_paracentral_part3_thickness","lh_parsopercularis_part1_thickness",
                   "lh_parsopercularis_part2_thickness","lh_parsopercularis_part3_thickness","lh_parsorbitalis_part1_thickness","lh_parstriangularis_part1_thickness","lh_parstriangularis_part2_thickness","lh_pericalcarine_part1_thickness","lh_pericalcarine_part2_thickness","lh_postcentral_part1_thickness","lh_postcentral_part2_thickness","lh_postcentral_part3_thickness","lh_postcentral_part4_thickness",
                   "lh_postcentral_part5_thickness","lh_postcentral_part6_thickness","lh_postcentral_part7_thickness","lh_postcentral_part8_thickness","lh_posteriorcingulate_part1_thickness","lh_posteriorcingulate_part2_thickness","lh_precentral_part1_thickness","lh_precentral_part2_thickness","lh_precentral_part3_thickness","lh_precentral_part4_thickness","lh_precentral_part5_thickness","lh_precentral_part6_thickness",
                   "lh_precentral_part7_thickness","lh_precentral_part8_thickness","lh_precentral_part9_thickness","lh_precuneus_part1_thickness","lh_precuneus_part2_thickness","lh_precuneus_part3_thickness","lh_precuneus_part4_thickness","lh_precuneus_part5_thickness","lh_precuneus_part6_thickness","lh_precuneus_part7_thickness","lh_rostralanteriorcingulate_part1_thickness","lh_rostralmiddlefrontal_part1_thickness",
                   "lh_rostralmiddlefrontal_part2_thickness","lh_rostralmiddlefrontal_part3_thickness","lh_rostralmiddlefrontal_part4_thickness","lh_rostralmiddlefrontal_part5_thickness","lh_rostralmiddlefrontal_part6_thickness","lh_rostralmiddlefrontal_part7_thickness","lh_rostralmiddlefrontal_part8_thickness","lh_rostralmiddlefrontal_part9_thickness","lh_rostralmiddlefrontal_part10_thickness","lh_superiorfrontal_part1_thickness",
                   "lh_superiorfrontal_part2_thickness","lh_superiorfrontal_part3_thickness","lh_superiorfrontal_part4_thickness","lh_superiorfrontal_part5_thickness","lh_superiorfrontal_part6_thickness","lh_superiorfrontal_part7_thickness","lh_superiorfrontal_part8_thickness","lh_superiorfrontal_part9_thickness","lh_superiorfrontal_part10_thickness","lh_superiorfrontal_part11_thickness","lh_superiorfrontal_part12_thickness",
                   "lh_superiorfrontal_part13_thickness","lh_superiorparietal_part1_thickness","lh_superiorparietal_part2_thickness","lh_superiorparietal_part3_thickness","lh_superiorparietal_part4_thickness","lh_superiorparietal_part5_thickness","lh_superiorparietal_part6_thickness","lh_superiorparietal_part7_thickness","lh_superiorparietal_part8_thickness","lh_superiorparietal_part9_thickness","lh_superiorparietal_part10_thickness",
                   "lh_superiortemporal_part1_thickness","lh_superiortemporal_part2_thickness","lh_superiortemporal_part3_thickness","lh_superiortemporal_part4_thickness","lh_superiortemporal_part5_thickness","lh_superiortemporal_part6_thickness","lh_superiortemporal_part7_thickness","lh_supramarginal_part1_thickness","lh_supramarginal_part2_thickness","lh_supramarginal_part3_thickness","lh_supramarginal_part4_thickness",
                   "lh_supramarginal_part5_thickness","lh_supramarginal_part6_thickness","lh_supramarginal_part7_thickness","lh_frontalpole_part1_thickness","lh_temporalpole_part1_thickness","lh_transversetemporal_part1_thickness","lh_insula_part1_thickness","lh_insula_part2_thickness","lh_insula_part3_thickness","lh_insula_part4_thickness","rh_bankssts_part1_thickness","rh_bankssts_part2_thickness","rh_caudalanteriorcingulate_part1_thickness",
                   "rh_caudalmiddlefrontal_part1_thickness","rh_caudalmiddlefrontal_part2_thickness","rh_caudalmiddlefrontal_part3_thickness","rh_caudalmiddlefrontal_part4_thickness","rh_cuneus_part1_thickness","rh_cuneus_part2_thickness","rh_cuneus_part3_thickness","rh_entorhinal_part1_thickness","rh_fusiform_part1_thickness","rh_fusiform_part2_thickness","rh_fusiform_part3_thickness","rh_fusiform_part4_thickness",
                   "rh_fusiform_part5_thickness","rh_inferiorparietal_part1_thickness","rh_inferiorparietal_part2_thickness","rh_inferiorparietal_part3_thickness","rh_inferiorparietal_part4_thickness","rh_inferiorparietal_part5_thickness","rh_inferiorparietal_part6_thickness","rh_inferiorparietal_part7_thickness","rh_inferiorparietal_part8_thickness","rh_inferiorparietal_part9_thickness","rh_inferiorparietal_part10_thickness",
                   "rh_inferiortemporal_part1_thickness","rh_inferiortemporal_part2_thickness","rh_inferiortemporal_part3_thickness","rh_inferiortemporal_part4_thickness","rh_inferiortemporal_part5_thickness","rh_isthmuscingulate_part1_thickness","rh_isthmuscingulate_part2_thickness","rh_lateraloccipital_part1_thickness","rh_lateraloccipital_part2_thickness","rh_lateraloccipital_part3_thickness","rh_lateraloccipital_part4_thickness",
                   "rh_lateraloccipital_part5_thickness","rh_lateraloccipital_part6_thickness","rh_lateraloccipital_part7_thickness","rh_lateraloccipital_part8_thickness","rh_lateraloccipital_part9_thickness","rh_lateralorbitofrontal_part1_thickness","rh_lateralorbitofrontal_part2_thickness","rh_lateralorbitofrontal_part3_thickness","rh_lateralorbitofrontal_part4_thickness","rh_lingual_part1_thickness","rh_lingual_part2_thickness",
                   "rh_lingual_part3_thickness","rh_lingual_part4_thickness","rh_lingual_part5_thickness","rh_lingual_part6_thickness","rh_medialorbitofrontal_part1_thickness","rh_medialorbitofrontal_part2_thickness","rh_medialorbitofrontal_part3_thickness","rh_middletemporal_part1_thickness","rh_middletemporal_part2_thickness","rh_middletemporal_part3_thickness","rh_middletemporal_part4_thickness","rh_middletemporal_part5_thickness",
                   "rh_middletemporal_part6_thickness","rh_parahippocampal_part1_thickness","rh_parahippocampal_part2_thickness","rh_paracentral_part1_thickness","rh_paracentral_part2_thickness","rh_paracentral_part3_thickness","rh_parsopercularis_part1_thickness","rh_parsopercularis_part2_thickness","rh_parsopercularis_part3_thickness","rh_parsorbitalis_part1_thickness","rh_parstriangularis_part1_thickness","rh_parstriangularis_part2_thickness",
                   "rh_parstriangularis_part3_thickness","rh_pericalcarine_part1_thickness","rh_pericalcarine_part2_thickness","rh_pericalcarine_part3_thickness","rh_postcentral_part1_thickness","rh_postcentral_part2_thickness","rh_postcentral_part3_thickness","rh_postcentral_part4_thickness","rh_postcentral_part5_thickness","rh_postcentral_part6_thickness","rh_postcentral_part7_thickness","rh_postcentral_part8_thickness","rh_posteriorcingulate_part1_thickness",
                   "rh_posteriorcingulate_part2_thickness","rh_precentral_part1_thickness","rh_precentral_part2_thickness","rh_precentral_part3_thickness","rh_precentral_part4_thickness","rh_precentral_part5_thickness","rh_precentral_part6_thickness","rh_precentral_part7_thickness","rh_precentral_part8_thickness","rh_precentral_part9_thickness","rh_precuneus_part1_thickness","rh_precuneus_part2_thickness","rh_precuneus_part3_thickness","rh_precuneus_part4_thickness",
                   "rh_precuneus_part5_thickness","rh_precuneus_part6_thickness","rh_precuneus_part7_thickness","rh_rostralanteriorcingulate_part1_thickness","rh_rostralmiddlefrontal_part1_thickness","rh_rostralmiddlefrontal_part2_thickness","rh_rostralmiddlefrontal_part3_thickness","rh_rostralmiddlefrontal_part4_thickness","rh_rostralmiddlefrontal_part5_thickness","rh_rostralmiddlefrontal_part6_thickness","rh_rostralmiddlefrontal_part7_thickness","rh_rostralmiddlefrontal_part8_thickness",
                   "rh_rostralmiddlefrontal_part9_thickness","rh_rostralmiddlefrontal_part10_thickness","rh_superiorfrontal_part1_thickness","rh_superiorfrontal_part2_thickness","rh_superiorfrontal_part3_thickness","rh_superiorfrontal_part4_thickness","rh_superiorfrontal_part5_thickness","rh_superiorfrontal_part6_thickness","rh_superiorfrontal_part7_thickness","rh_superiorfrontal_part8_thickness","rh_superiorfrontal_part9_thickness","rh_superiorfrontal_part10_thickness","rh_superiorfrontal_part11_thickness",
                   "rh_superiorfrontal_part12_thickness","rh_superiorfrontal_part13_thickness","rh_superiorparietal_part1_thickness","rh_superiorparietal_part2_thickness","rh_superiorparietal_part3_thickness","rh_superiorparietal_part4_thickness","rh_superiorparietal_part5_thickness","rh_superiorparietal_part6_thickness","rh_superiorparietal_part7_thickness","rh_superiorparietal_part8_thickness","rh_superiorparietal_part9_thickness",
                   "rh_superiorparietal_part10_thickness","rh_superiortemporal_part1_thickness","rh_superiortemporal_part2_thickness","rh_superiortemporal_part3_thickness","rh_superiortemporal_part4_thickness","rh_superiortemporal_part5_thickness","rh_superiortemporal_part6_thickness","rh_supramarginal_part1_thickness","rh_supramarginal_part2_thickness","rh_supramarginal_part3_thickness","rh_supramarginal_part4_thickness","rh_supramarginal_part5_thickness","rh_supramarginal_part6_thickness",
                   "rh_supramarginal_part7_thickness","rh_frontalpole_part1_thickness","rh_temporalpole_part1_thickness","rh_transversetemporal_part1_thickness","rh_insula_part1_thickness","rh_insula_part2_thickness","rh_insula_part3_thickness","rh_insula_part4_thickness")
    
    
    #######
    #  Ordering the 308 regions to save in txt file
    ratio_ordered <- left_join(data.frame(region = reg_order),
                               dataf, 
                               by = "region")
    
    ratio_ordered[,"ratio"] <- -log(ratio_ordered[,"ratio"]) 
    # - log(0) = inf --> replace inf by 0.0
    ratio_ordered <- do.call(data.frame,lapply(ratio_ordered, function(x) replace(x, is.infinite(x),0)))
    write.table(ratio_ordered[,"ratio"], 
                paste0("/data_J/Results/Files/GlobalRatioOrdered_",tp,"_", lab,".txt"), 
                sep = "\t", 
                row.names = FALSE,
                col.names = FALSE)
  }
}
