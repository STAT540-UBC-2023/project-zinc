# Project Proposal 
This is the project proposal for the group zinc. Everyone has read over the proposal and edited everything. Final proposal .md file made, committed, and pushed by Janet. 

## Motivation, Background and Hypothesis (By Lucy)

Type-2 diabetes (T2D) is a common metabolic disorder characterized by an elevation of blood glucose levels due to the lack of functional insulin [1]. Many patients with T2D are commonly portrayed as having accumulated higher fat percentages in their bodies. Obesity is commonly defined as Body Mass Index (BMI) larger than 30 kg/m2 [1,2]. The metabolic complications of obese individuals consist of insulin resistance, culminating in pancreatic cell overload and failure [3], which leads to T2D [4,5]. However, obesity is not the only factor for T2D, as a large population of non-obese T2D patients has also been rising [5]. The driving force for T2D in non-obese patients leans towards a genetic factor, resulting in defective insulin production [6]. Therefore, it is critical to understand the pathogenesis of both conditions for therapies developments [6]. 

The pancreas tightly regulates metabolic homeostasis through its exocrine and endocrine systems, which contains the islets of Langerhans, including 5 cell types (α, β, γ/PP, δ, and ε) [7]. In this paper, Segerstolpe et al. have identified gene expression alterations in obese individuals (vs BMI<30) and T2D patients (vs healthy), respectively [7]. However, given the complication between obesity and T2D, the correlation between BMI and T2D remains to be elucidated. 

Based on the different genetic mechanisms causing two forms of T2D, we hypothesize that BMI, age and or sex are associated with differential gene expression profiles in T2D and healthy subjects. Secondly, the non-obese and obese T2D individuals are expected to have different changes in their genetic compositions.

## Data (By Cindy)

Single-cell RNA sequences (scRNA-seq) [E-MTAB-5061](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-5061?accession=E-MTAB-5061) were harvested from pancreatic islets of post-mortem individuals where 1534 had diabetes and 1980 did not have diabetes, using Illumina HiSeq 2000 sequencing contains 43 bp single-end reads [7]. The scRNA-seq data matrix contains 26271 genes (rows) x 7028 samples (cols), where 3514 cols contain normalized expression values in RPKM and 3514 cols contain raw transcript counts. The corresponding phenotype feature table contains 3514 samples (rows) x 49 characteristics (cols). Bulk RNA-seq [E-MTAB-5060](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-5060?accession=E-MTAB-5060) were harvested from the  whole pancreatic islets of 3 healthy and 4 T2D donors sequenced with Illumina NextSeq 500. The expression matrix consists of  26271 genes (rows) x 14 samples (cols) with the corresponding feature metadata that contains 14 samples (rows) x 40 characteristics (cols). 

The bulk RNA sequencing data enables us to understand the relationship between variables such as BMI and differential gene expression in T2D or healthy subjects via linear regression. The scRNA-seq data allows us to visualize the differential expression at the single cell level. Comparisons can therefore be made for these single cells to evaluate cell-type specific responses and identify cell-markers.

## Analysis Done Previously ~~(By Shunsuke)~~

Here, we briefly summarize the analysis procedures of previous work, to clarify the differences between the original study and our own objectives and corresponding methods. Segerstolpe et al. analyzed the granularity of individual genes of identified cell types include: reporting correlations by Spearman's rho between genetic expressions and obesity in the group of the males without T2D; and finding significantly different gene expressions, comparing the two cohorts by ANOVA [7]. 

## Aim and Analysis Plan 

**Aim 1** *(Cindy)* To understand the effect of BMI, age, or sex on differential gene expressions in T2D and healthy cohorts, the association between BMI and differential gene expression in T2D and healthy individuals will be examined using bulk RNAseq data. The [edegR](https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf) will be employed for normalization of count matrix. The expression data in T2D and healthy individuals will each be subjected to a linear regression model using [limma](https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf), with BMI, sex, and age as explanatory variables. Collinearity will be examined using Spearman’s coefficient.

**Aim 2** *(Lucy&Janet)* Single cell RNA-seq data will be used to perform differential gene expression analysis on obese and non-obese T2D patients. To reject the null hypothesis that differential gene expressions between these patients are insignificant, we will be using [Seurat](https://cran.r-project.org/web/packages/Seurat/index.html) for quality control and data filtering, and [DESeq-2](http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html) for a quantitative analysis of comparative scRNA-seq data. A table of significantly expressed genes (p-adjusted threshold of 0.05) will be provided. Volcano plots are generated to help visualize the significantly up and down-regulated genes in both forms of T2D patients when compared to healthy controls. A comparative two-sample t-test analysis is then conducted on genes with significant expressions from both groups to verify whether the expression profiles are different. 

## Division of Labour 

Group Member|Department|Background/Degree/Affiliation|Division of Labour
------------|----------|---------|----------
Janet(Jia He) Zhang|Pharmacology|BSc.(Pharmacology)|1. Github management 2. Single cell RNA-seq data processing,quality control, and analysis 
Lucy Chi|Biochemistry|BSc. Honour(Biochemistry & Molecular Biology)|1. Github management 2.  literature search on biological background 3. Differential gene analysis (volcano plot for obese vs non-obese T2D)
Cindy (Xiao Yu) Zhang|Bioinformatics|BSc.(Pharm), PharmD, RPh|1. Github management 2. Initial data selection and wrangling 3. Examine relationship between BMI, age, sex and differential gene expression in T2D and healthy subjects using a linear regression model
Shunsuke Ishige|Interdisciplinary Oncology|BSc.(Mathematics, with Philosophy minor), BSc.(Computer Science)|1. Literature survey on methods for quantitative comparison of genetic expression differences 2. Data processing 3. Search for larger data sets(e.g., mice) which we can apply our analysis (if needed) 4. Study data with summary statistics 5. Provide support for coding issues (if needed) 6. Provide typeset drafts in LaTeX (if needed)

## References (By Lucy)

1. Galicia-Garcia, Unai et al. “Pathophysiology of Type 2 Diabetes Mellitus.” *International journal of molecular sciences* vol. 21,17 6275. 30 Aug. 2020, doi:10.3390/ijms21176275
2. Al-Goblan, Abdullah S et al. “Mechanism linking diabetes mellitus and obesity.” *Diabetes, metabolic syndrome and obesity : targets and therapy* vol. 7 587-91. 4 Dec. 2014, doi:10.2147/DMSO.S67400
3. Singla, Parul et al. “Metabolic effects of obesity: A review.” *World journal of diabetes* vol. 1,3 (2010): 76-88. doi:10.4239/wjd.v1.i3.76
4. Prentki, Marc, and Christopher J Nolan. “Islet beta cell failure in type 2 diabetes.” *The Journal of clinical investigation* vol. 116,7 (2006): 1802-12. doi:10.1172/JCI29103
5. Olaogun, Idowu et al. “The Pathophysiology of Type 2 Diabetes Mellitus in Non-obese Individuals: An Overview of the Current Understanding.” *Cureus* vol. 12,4 e7614. 10 Apr. 2020, doi:10.7759/cureus.7614
6. Vaag, Allan, and Søren S Lund. “Non-obese patients with type 2 diabetes and prediabetic subjects: distinct phenotypes requiring special diabetes treatment and (or) prevention?.” *Applied physiology, nutrition, and metabolism = Physiologie appliquee, nutrition et metabolisme* vol. 32,5 (2007): 912-20. doi:10.1139/H07-100
7. Segerstolpe, Åsa et al. “Single-Cell Transcriptome Profiling of Human Pancreatic Islets in Health and Type 2 Diabetes.” *Cell metabolism* vol. 24,4 (2016): 593-607. doi:10.1016/j.cmet.2016.08.020


