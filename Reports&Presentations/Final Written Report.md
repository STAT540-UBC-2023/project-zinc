# A Glimpse into Type II Diabetes(T2D) at the Gene Expression Level 
_Zinc Group_

### Background

Type-2 diabetes (T2D) is a common metabolic disorder characterized by an elevation of blood glucose levels due to the lack of functional insulin [1]. This is due to either the impairment of insulin secretion or the development of insulin resistance by functional cells [1]. Insulin is an important polypeptide hormone responsible for regulating blood glucose level and inducing glucose storage in energy-dependent cells and tissues, such as muscles and liver [2]. The impairment of this system will lead to many symptoms, including fatigue, thirst, and vision loss [3]. Unfortunately, symptoms for early stage T2D are not very obvious, so T2D diagnosis might take several years after onset [3]. When left untreated, it can lead to complications such as kidney impairment, cardiovascular diseases and blindness [4]. 

To date, there is no effective cure. However, many researchers have suggested cost-effective intervention to prevent serious outcome, including physical exercise, and diet control [3]. One of the reasons is because they increasingly noticed a strong association between obesity and diabetes, shown in many published papers [5,6]. Many patients with T2D are commonly portrayed as having accumulated higher fat percentages in their bodies [1]. Obesity is commonly defined as Body Mass Index (BMI) larger than 30 kg/m2 [1,7]. The major metabolic complications of obese individuals consist of insulin resistance, culminating in pancreatic cell overload and failure [8], which leads to T2D [9,10]. However, obesity is not the only factor for T2D, as a large population of non-obese T2D patients has also been rising [10]. The driving force for T2D in non-obese patients leans towards a genetic factor, resulting in defective insulin production [11]. Therefore, it is critical to understand the pathogenesis of both conditions for therapies developments [11].

Since clinical observation shows obesity as a contributing factor for some T2D patience but not the others, it brings up interesting biological questions: does obesity contribute to gene expression change for T2D patients? If so, what is the differential gene expression between obese and non-obese T2D individuals? 

### Hypothesis and Aims

We hypothesize that T2D and obesity are associated with differential gene expression. Moreover, the non-obese and obese T2D individuals are expected to have different gene expression due to the different mechanisms proposed above. 

There are two aims to this analysis. 
+ Aim1: perform linear regression for assessing the effect of T2D and BMI on differential gene expression. 
+ Aim2: perform gene set enrichment analysis with Reactome, Pathway Interaction Database (PID), ImmuneSigDB, and Kyoto Encyclopedia of Genes and Genomes (KEGG) databases for differentially expressed genes. Additionally, we will use DisGeNET to study the types of diseases associated with T2D. 

### Datasets 

To answer the biological question,  the datasets from [GSE41762](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41762) were chosen. This dataset contains microarray data for RNA harvested from islets of cadaver donors and was previously described in Pubmed articles, [PMID23140642](https://pubmed.ncbi.nlm.nih.gov/23140642/) and [PMID25298321](https://pubmed.ncbi.nlm.nih.gov/25298321/)[12,13]. There are 77 samples in total - 53 samples with non-T2D and BMI < 30, 4 samples with non-T2D and BMI > 30, 14 samples with T2D and BMI < 30 and 6 samples with T2D and BMI >30, with expression data for samples 1 to 48 obtained from one experiment and 49 through 77 from another. 

### Methods 

#### [Linear Regression](https://github.com/STAT540-UBC-2023/project-zinc/tree/main/MicroarrayLinearRegression)

###### 1.1 Data Wrangling 

The Robust Multi-Array (RMA) normalized expression matrix and metadata were downloaded from [GSE41762](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41762). The variables of interest were subsetted from metadata. Upon examination, we observed 366 genes and 29 samples with missing data from sample 49 to 77 in the original dataset. However, no detailed description of the discrepancy was found from literature. Next, singular vector decomposition plots were used to visualize any potential batch effects between samples gathered from the first (sample 1 to 48) and second (sample 49 to 77) experiment in this dataset. As shown in figure1 ([1A](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/MicroarrayLinearRegression/MicroarrayLinearRegressionSrc_files/figure-gfm/unnamed-chunk-7-1.png), [1B](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/MicroarrayLinearRegression/MicroarrayLinearRegressionSrc_files/figure-gfm/unnamed-chunk-7-2.png), [1C](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/MicroarrayLinearRegression/MicroarrayLinearRegressionSrc_files/figure-gfm/unnamed-chunk-7-3.png)), a moderate batch effect was observed between the two experiments. Two approaches have been considered in dealing with the batch effect: 

1. Harmonizing data from both experiments using ComBat 
2. Samples from the two experiments were analyzed separately. 

The second approach, albeit more cumbersome, was deemed to reduce risk of inducing bias in data merging, which is limited by unbalanced sample sizes of other confounders under different protocols. Hence, the second approach was chosen. 

###### 1.2 Model Fitting 

The number of samples in each category of experiment 1 and 2 were shown in table 2 and 3, respectively. Since none of the sample had BMI >30, for experiment 1 (sample 49:77) a simple linear model was used with  ~T2D, while experiment 2 (sample 1:48) was fitted with an interaction model, ~T2D* BMI. The [LIMMA](https://bioconductor.org/packages/release/bioc/html/limma.html) package in R was used for linear model fitting to the gene expression data. An empirical Bayes moderation of the standard errors was performed for each model and for the calculation computing of moderated t-statistics. The LIMMA package was well-suited for the purpose of this analysis due to its robustness to noise and statistical power in detecting subtle differences in microarray experiments. Statistical significance was defined as FDR<0.05 by Benjamini-Hochberg Procedure. 

Table 2. Dataset for samples 1:48

+---------------+---------------+--------------------+
|               | BMI<30        | BMI>30             |
+===============+===============+====================+
| non-T2D       | 32            | 0                  |
+---------------+---------------+--------------------+
| T2D           | 8             | 2                  |
+---------------+---------------+--------------------+

Table 3. Dataset for samples 49:77

+---------------+---------------+--------------------+
|               | BMI<30        | BMI>30             |
+===============+===============+====================+
| non-T2D       | 15            | 4                  |
+---------------+---------------+--------------------+
| T2D           | 6             | 4                  |
+---------------+---------------+--------------------+
*Result interpretation please refer to the result section*

#### [Gene Set Enrichment Analysis](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/GeneSetEnrichmentAnalysis/Gene-Set-Enrichment-Analysis.md) 


###### 2.1 Differentially Expressed (DE) Gene ID Conversion and Clean-up

The DE gene were expressed in the transcript cluster IDs form (labeled as “affy_hugene_1_0_st_v1”), which was not compatible with the required IDs for downstream functional analysis. [`biomaRt`](https://bioconductor.org/packages/release/bioc/html/biomaRt.html) package was used to obtained the corresponding ensembl gene IDs, gene symbols, and entrez gene IDs. However, many gene duplicates appeared. R code used the `order` function to sort the rows of gene list by the log_adjP column in descending order. The `head` function is then used to select the top row for each unique gene_symbol value. Similar process was also done using log fold change values for further KEGG analysis. 

###### 2.2 Functional Analysis 

Prior to performing functional analysis, Genome-wide Association Studies (GWAS) was adopted to obtain a gene list that is associated with T2D to assist future analysis. 

Molecular Signature Database (MsigDB) was obtained [14] through [`msigdbr`](http://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp) package in R. Four gene sets (Reactome, Pathway Interaction Database (PID), ImmuneSigDB subcollection of the Immunologic Signature collection, and Kyoto Encyclopedia of Genes and Genomes (KEGG)) were selected to perform functional analysis. 

Specifically, the [`fgsea`](https://bioconductor.org/packages/release/bioc/html/fgsea.html) package in R was used to carry out all ranked-based functional gene set enrichment analysis (fgsea) and produce a list of significant pathways associated to T2D (FDR < 0.05). The [`enrichplot`](https://github.com/YuLab-SMU/enrichplot) and [`ggfortify`](https://cran.r-project.org/web/packages/ggfortify/index.html) packages were used to implement different visualization methods for the results from GSEA analysis. Particularly, all the significant pathways were represented with bar plots using `geom_bar` in ascending order of adjusted p values. While analyzing the KEGG gene sets, additional R packages, [`clusterProfiler`](https://rdrr.io/bioc/clusterProfiler/) and [`ggridge`](https://rdrr.io/cran/ggridges/) were used for graphic visualization of functional profiles (KEGG) for genes and gene clusters. The [`org.Hs.eg.db`](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html) package is an organism specific package used to obtain genome wide annotation for human, primarily based on mapping by Entrez gene IDs. The `gseKEGG` function within the `clusterProfiler` package allowed for a rank-based gene set enrichment analysis. And along with the `enrichplot` package, an enrichment map and ridgeplot for more detailed visualization of the KEGG GSEA analysis were drawn. The genes associated with resulted pathways were compared with the T2D associated gene list obtained from GWAS for result verification. 

###### 2.3 Disease Association

The [‘DisGeNet’](https://www.disgenet.org/home/) database combines information of human genetic disease associations and variant-disease associations together. The [`DOSE`](https://bioconductor.org/packages/release/bioc/html/DOSE.html) package is used to provide an enrichment analysis based on the DisGeNET database. The results were visualized with an enriched bar graph. 

### Results 

#### Linear Regression 

For the first experiment (sample 1 to 48), 67 genes were found with significant differential expression (FDR<0.05) in T2D vs healthy samples using the simple regression model, as shown in [figure 2](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/MicroarrayLinearRegression/MicroarrayLinearRegressionSrc_files/figure-gfm/unnamed-chunk-13-1.png).  For the second experiment (sample 49 to 77), however, no genes passed the significance cut-off using the interactive model. The linear regression result for experiment 1 showed that T2D is associated with differential gene expression, which answered the first half of the biological question. However, the result from model 2 failed to reject the null hypothesis that BMI is not associated with differential gene expression. No significant differences in gene expression were found between obese and non-obese individuals with T2D. One biological explanation is that due to the controversial definition of obesity caused by the current lack of clinical measurements, researchers have grouped individuals with BMI>30 to be "obese" and others as "non-obese". Many research groups have since argued that even with BMI values less than 30, these individuals can still be "metabolically obese" [15]. This suggests there may be better alternatives for categorizing obesity, which could explain why the BMI factor in our analysis showed only one DE gene [15]. Additionally, looking back at our dataset, this result could also be due to the unbalance in samples where 21 samples had BMI<30 and only 8 samples had BMI>30. A larger dataset would be needed for further exploration. Given these results, we decided to only focus on the differentially expressed genes between T2D and healthy individuals in the following gene set enrichment analysis. 

#### Gene Set Enrichment Analysis 

Upon assessing, we decided to perform functional gene set enrichment analysis (fgsea) to further analyze the DE genes with a rank-based method. One benefit is that this ranking method provides a normalization tool for all genes with various expression levels before the enrichment analysis. Although multiple databases were considered, four databases were chosen for further examinations based on the compatibility to our research interest. 

###### The REACTOME Database

First, the REACTOME pathway was used to identify biological pathways associated with T2D. REACTOME provides a biological pathway database and detailed information on molecular events involved in human biological processes [16]. The fgsea analysis with this database showed 82 associated pathways, most of which were consistent with previous literature [17]. For example, PI3K-AKT signaling pathway plays a role in glucose metabolism and insulin signaling with a p-value of 0.000021. Another significant pathway is signaling by receptor tyrosine kinases (RTKs), which involves insulin receptor substrate (IRS) pathway important for T2D development. Interestingly, the bar plot suggested that other than the cellular regulations, many immune pathways are involved, especially the interleukin system, shown in [figure 3](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/GeneSetEnrichmentAnalysis/Gene-Set-Enrichment-Analysis_files/figure-gfm/unnamed-chunk-36-1.png). It is thought that chronic damage of pancreatic beta cell will lead to hyperglycemia caused by the insufficient insulin production. This will subsequently induce immune response impairment [18], explaining the high involvement of many immune pathways in T2D patients. 

###### The PID Database

The [pathway interaction database](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2686461/) presents molecular interactions and events that compose of key cellular events. This database is designed to capture the biological knowledge at different levels of details. Based on fgsea analysis, 25 significant pathways were identified, of which 10 were associated with the interleukin system (example). These results were consistent with those from REACTOME. To further understand this, we decided to implement the immuneSigDB database and identify specific immune pathways associated with Type II diabetes. 

