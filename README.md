# STAT540 Group Project

Welcome to the repository for the Zinc group! For this group project, we initially chose this [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5069352/). Data can be found on ArrayExpress with accession numbers: [E-MTAB-5060](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-5060) and [E-MTAB-5061](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-5061?accession=E-MTAB-5061).

We changed our dataset to this one: [GSE41762](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41762). This dataset was used in two papers, one by [Mahdi et al.,(2012)](https://pubmed.ncbi.nlm.nih.gov/23140642/), and the other by [Tang et al., (2014)](https://pubmed.ncbi.nlm.nih.gov/25298321/).

## Team members:

-   Cindy (XiaoYu) Zhang
-   Lucy (Shuxin) Chi
-   Janet (Jia He) Zhang

## Table of Content: 

1. [Microarray Linear Regression](https://github.com/STAT540-UBC-2023/project-zinc/tree/main/MicroarrayLinearRegression) Directory: contains everything (e.g., .Rmd file, .md file, results and plots) generated for the microarray linear regression analysis. 
      + [linear-regression.rmd](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/MicroarrayLinearRegression/MicroarrayLinearRegressionSrc.Rmd)
      + [linear-regression.md](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/MicroarrayLinearRegression/MicroarrayLinearRegressionSrc.md)
      + [Figures](https://github.com/STAT540-UBC-2023/project-zinc/tree/main/MicroarrayLinearRegression/MicroarrayLinearRegressionSrc_files/figure-gfm)
2. [Gene Set Enrichment Analysis](https://github.com/STAT540-UBC-2023/project-zinc/tree/main/GeneSetEnrichmentAnalysis) Directory: contains everything(e.g., .Rmd file, .md file, results and plots) generated for the gene set enrichment analysis. 
      + [Gene-set-enrichment-analysis.rmd](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/GeneSetEnrichmentAnalysis/Gene-Set-Enrichment-Analysis.Rmd)
      + [Gene-set-enrichment-analysis.md](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/GeneSetEnrichmentAnalysis/Gene-Set-Enrichment-Analysis.md)
      + [Figures](https://github.com/STAT540-UBC-2023/project-zinc/tree/main/GeneSetEnrichmentAnalysis/Gene-Set-Enrichment-Analysis_files/figure-gfm)
      + [Results](https://github.com/STAT540-UBC-2023/project-zinc/tree/main/GeneSetEnrichmentAnalysis/Result)
3. [Reports & Presentations](https://github.com/STAT540-UBC-2023/project-zinc/tree/main/Reports%26Presentations) Directory: contains all the reports and presentation slides (e.g., lightning talk, written proposal) created for this project . 

## Summary 

#### Batch Effect Evaluation 
+ Aim: To examine batch-effect in assessing how suitabl it is to combine data from two experiments. Please check the [MicroarrayLinearRegression.md](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/MicroarrayLinearRegression/MicroarrayLinearRegressionSrc.md) for more details. 
+ Results: Moderate batch effect was observed in SVD plots, as shown in these figures ([1](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/MicroarrayLinearRegression/MicroarrayLinearRegressionSrc_files/figure-gfm/unnamed-chunk-7-1.png), [2](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/MicroarrayLinearRegression/MicroarrayLinearRegressionSrc_files/figure-gfm/unnamed-chunk-7-2.png), [3](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/MicroarrayLinearRegression/MicroarrayLinearRegressionSrc_files/figure-gfm/unnamed-chunk-7-3.png)). 

#### Linear Model Fitting 
+ Aim: To examine the differentially expressed genes across a variety of conditions (e.g., T2D status, BMI, and T2D:BMI interaction). Detail analysis can be found in [MicroarrayLinearRegression.md](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/MicroarrayLinearRegression/MicroarrayLinearRegressionSrc.md) for more details. 
+ Results: For the first experiment (sample 1 to 48), 67 genes were found with significant differential expression (FDR<0.05) in T2D vs healthy samples using the simple regression model, as shown in figure 2.  For the second experiment (sample 49 to 77), however, no genes passed the significance cut-off using the interactive model. The linear regression result for experiment 1 showed that T2D is associated with differential gene expression, which answered the first half of the biological question. However, the result from model 2 failed to reject the null hypothesis that BMI is not associated with differential gene expression. No significant differences in gene expression were found between obese and non-obese individuals with T2D.

#### Gene Set Enrichment Analysis 
+ Aim: To determine if a priori defined set of genes (also called gene sets) show statistically significant difference between T2D VS Healthy with the [R markdown file](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/GeneSetEnrichmentAnalysis/Gene-Set-Enrichment-Analysis.Rmd)(and the relevant [markdown file](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/GeneSetEnrichmentAnalysis/Gene-Set-Enrichment-Analysis.md)). 
+ Results: results of different gene sets analyzed can be found in the [Results](https://github.com/STAT540-UBC-2023/project-zinc/tree/main/GeneSetEnrichmentAnalysis/Result) directory. This directory contains a number of tables containing the significant pathways found in different gene sets. For visualizations, please check this [Figure](https://github.com/STAT540-UBC-2023/project-zinc/tree/main/GeneSetEnrichmentAnalysis/Gene-Set-Enrichment-Analysis_files/figure-gfm) directory. 

## Final Report 

Please check out our final written report [here](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/Reports%26Presentations/Final%20Written%20Report.md). 

## Final Presentation 

Please check out our final presentation slides [here](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/Reports%26Presentations/STAT540%20final%20project.pptx). 
