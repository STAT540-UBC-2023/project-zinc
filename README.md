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
3. [Reports & Presentations](https://github.com/STAT540-UBC-2023/project-zinc/tree/main/Reports%26Presentations) Directory: contains all the reports and presentation slides (e.g., lightning talk, written proposal) created for this project . 

## Summary 

#### Principal Componenet Analysis
+ Aim: To examine batch-effect. Please check the [MicroarrayLinearRegression.md](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/MicroarrayLinearRegression/MicroarrayLinearRegressionSrc.md) for more details. 
+ Results: Moderate batach effect for U2 vsU1 and U3 vs U1. Visualize results in these figures ([1](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/MicroarrayLinearRegression/MicroarrayLinearRegressionSrc_files/figure-gfm/unnamed-chunk-7-1.png), [2](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/MicroarrayLinearRegression/MicroarrayLinearRegressionSrc_files/figure-gfm/unnamed-chunk-7-2.png), [3](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/MicroarrayLinearRegression/MicroarrayLinearRegressionSrc_files/figure-gfm/unnamed-chunk-7-3.png)). 

#### Linear Model Fitting 
+ Aim: To examine the diferentially expressed genes across a variety of conditions (e.g., T2D status, BMI, and T2D;BMI interaction). 
+ Results: 

#### Gene Set Enrichment Analysis 
+ Aim: To determine if a priori defined set of genes (also called gene sets) show statistically significant difference between T2D VS Healthy with the [R markdown file](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/GeneSetEnrichmentAnalysis/Gene-Set-Enrichment-Analysis.Rmd)(and the relevant [markdown file](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/GeneSetEnrichmentAnalysis/Gene-Set-Enrichment-Analysis.md)). 
+ Results: results of different gene sets analyzed can be found in the [Results](https://github.com/STAT540-UBC-2023/project-zinc/tree/main/GeneSetEnrichmentAnalysis/Result) directory. This directory contains a number of tables containing the significant pathways found in different gene sets. For visualizations, please check this [Figure](https://github.com/STAT540-UBC-2023/project-zinc/tree/main/GeneSetEnrichmentAnalysis/Gene-Set-Enrichment-Analysis_files/figure-gfm) directory. 

## Final Report 

Please check out our final written report [here](). 

## Final Presentation 

Please check out our final presentation slides [here](). 
