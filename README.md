# UC3M Multivariate Analysis Projects

Repository with all the `R` Code of the project of Multivariate Analysis

# Details of folders:

--------
## data
--------

* Raw data with information of clinical trials.

--------
## rda
--------

* After wranglinf, this is the oficial dataset used in the project `rda/clinical_trial_complete.rda`

--------
## figure_output
--------

`.png` with some plot results. Avoid big size plot taking a print-screen of the graph and using as image.

--------
## Scripts
--------
* 1-data-wrangling.R: Wrangling of raw data set (raw data is in file "data")

* 2-visual_Analysis: file with the exploration of the data set
  
* 3-Sample_Mean_Analysis: file with the sample estimator and outliers study
  
* 4-PCA: file with the PCA study by the three categories: Sex, Pain and death

--------
## Report
--------

Details of the different reports for this projects. All folders contains the `.rmd` and `.pdf` format


1. **Report1-Data_Presentation**

First exploratory analysis of the data set.

2. **Report2-PCA**

3.1) File "Project_Report.Rmd" is the first exploration of the data set

3.2) File "Step1_Project.Rmd" is the first study of the data set:
  
a) Visual Analysis
  
b) Sample Estimators
  
c) PCA
    

3. **Report3-Cluster_Factor_Analysis**

Step 2 of the proyect!!
 
------------ 
# Life Lesson:
------------

First time I don't use the root directory of the project in Markdown. Consider two aspects:

* Change th e configuration of the root chunk in the section of global option of Markdown. Reference: [rmardown-cookbook](https://bookdown.org/yihui/rmarkdown-cookbook/working-directory.html)

* For the images, moves ones directory back. Reference: [stackoverflow](https://stackoverflow.com/questions/48813596/moving-one-directory-backward-in-r-path)
 
 