## Data

This project utilizes the Utrecht Schizophrenia Project database (293 healthy controls, 168 schizophrenia patients). -it is a longitudinal database (baseline, 4-year follow-up, sometimes second follow-up image at 7.5 years). The database comprises 1075 images processed using FreeSurfer, providing cortical thickness, brain volume, and Euler number data for 68 or 308 regions. Additionally, demographic variables including age, sex, diagnosis, IQ, and symptom presence are collected.

![](images/atlas.png)



## Project structure
This repository contains a folder with four R functions:

- 1_DataPreparation.R: functions for automating data preprocessing tasks to ensure data consistency and quality before analysis.
- 2_RegressionModel.R: functions for performing normative modeling and regression analysis on neuroimaging data, generating Z-scores and Cohen's d effect sizes for specified measures and regions of interest, and includes functions for model evaluation and correction for false discovery rate (FDR) in p-values.
- 3_Statistics.R: functions for computing statistics and generating density plots, histograms, boxplots, and mean z-scores of brain regions, saving the results as PNG files.
- 4_EDA.R: functions for performing exploratory data analysis (EDA) and variance partitioning using linear mixed-effects and fixed-effects models, generating violin plots for each data group and saving them as PNG files.

and the rest of the scripts that are .Rmd formatted, which use those functions.

## Results
![](images/models.png)
![](images/zs_s1_s2.png)
![](images/deviance-over-time.png)
![](images/deviance-over-time_1.png)












