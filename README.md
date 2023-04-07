# STAT540 Group Project

Welcome to the repository for the Zinc group! For this group project, we initially chose this [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5069352/). Data can be found on ArrayExpress with accession numbers: [E-MTAB-5060](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-5060) and [E-MTAB-5061](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-5061?accession=E-MTAB-5061).

We changed our dataset to this one: [GSE41762](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41762). This dataset was used in two papers, one by [Mahdi et al.,(2012)](https://pubmed.ncbi.nlm.nih.gov/23140642/), and the other by [Tang et al., (2014)](https://pubmed.ncbi.nlm.nih.gov/25298321/).

## Team members:

-   Cindy (XiaoYu) Zhang
-   Lucy(Shuxin) Chi
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
3. [Reports & Presentations](https://github.com/STAT540-UBC-2023/project-zinc/tree/main/Reports%26Presentations): contains all the reports and presentation slides (e.g., lightning talk, written proposal) created for this project . 

## Summary 

#### Data Wrangling 
Aim: To re-shape, filter, and wrangle data, and to explore raw data. 
Conclusion: there were 48 samples from one experiment cohort, and the other 29 samples from another, so these two experiments both contributed to the dataset. Please check the [MicroarrayLinearRegression.md](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/MicroarrayLinearRegression/MicroarrayLinearRegressionSrc.md). 

#### Principal Componenet Analysis
Aim: To examine batch-effect. 
Results: Moderate batach effect for U2 vsU1 and U3 vs U1. Visualize results in these figures ([1](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/MicroarrayLinearRegression/MicroarrayLinearRegressionSrc_files/figure-gfm/unnamed-chunk-7-1.png), [2](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/MicroarrayLinearRegression/MicroarrayLinearRegressionSrc_files/figure-gfm/unnamed-chunk-7-2.png), [3](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/MicroarrayLinearRegression/MicroarrayLinearRegressionSrc_files/figure-gfm/unnamed-chunk-7-3.png)). 
Conclusion: moderate batch effect, so two experiments are analyzed separately. 

#### 


## Final Report 

## Final Presentation 
