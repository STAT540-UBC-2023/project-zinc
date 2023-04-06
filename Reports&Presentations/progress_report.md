# Progress Report

## WHAT HAS CHANGED BASED ON THE FINAL PROPOSAL

#### Change in Dataset

The original dataset ([E-MTAB-5060](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-5060)) only had a total of 7 samples with only 1 sample in the T2D+BMI\<30 group (1). Since our hypothesis focuses on whether BMI and T2D status affect differential gene expression, that dataset did not have enough samples to answer our biological question. Additionally, the single cell RNA sequencing data from this paper is no longer feasible for our analysis. The main reason is because our hypothesis focuses on the effects of BMI and T2D on the population level, while single cell data only provide information on the individual level, which means that they are less representative and would not be able to help us answer the biological question.

A new dataset with more samples was therefore chosen: [GSE41762](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41762). This dataset contains microarray data for RNA harvested from islets of cadaver donors and was previously described in Pubmed articles ([PMID23140642](https://pubmed.ncbi.nlm.nih.gov/23140642/) and [PMID25298321](https://pubmed.ncbi.nlm.nih.gov/25298321/)) (2, 3). There are 77 samples in total - 53 samples with non-T2D and BMI \< 30, 4 samples with non-T2D and BMI \> 30, 14 samples with T2D and BMI \< 30 and 6 samples with T2D and BMI \>30.

Specifically, this dataset was generated by the Rosengren Lab from Lund University and used by two of their published papers, Mahdi et al. and Tang et al.(2, 3). Our study is different from previous studies in both the methodology and hypothesis: in the study by [Mahdi et al., (2012)](https://www.cell.com/cell-metabolism/fulltext/S1550-4131(12)00409-3?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS1550413112004093%3Fshowall%3Dtrue), the authors explored the pathophysiology of T2D by analysing the gene coexpression to prove their hypothesis that their pre-defined gene module is associated with T2D, elevated HbA1C and reduced insulin secretion (2). They performed PCA to identify these eignegenes and subsequently a gene set enrichment on these genes in answering this biological question. [Tang et al., (2014)](https://www.science.org/doi/10.1126/scitranslmed.3009934?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed), on the other hand, fitted an additive linear regression model with covariates: age, sex and BMI to investigate the effect of a genetic variant in *ADRA2*, encoding the α2A-adrenergic receptor (α2AAR), on insulin secretion (3).

#### Change in Analysis Plan

Principal Component Analysis (PCA)(see details in [.Rmd file](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/MicroarrayLinearRegressionSrc.Rmd)) was added as exploratory analysis of the dataset prior to the regression analysis.

Also, as mentioned above, the single cell data was disregarded due to the incompatibility with our biological question. Subsequently, the original plan on performing the differential analysis on gene expression using single cell RNA sequencing data (aim2) was removed. Instead, we focus on this new dataset and extend from the analysis performed in aim1.

Specifically, aim1 identified genes that are differentially expressed (DE) under the effect of BMI, T2D, and their interactions (\~BMI\*T2D) (more details about this aim can be found below in the Results section). Aim 2 is to take the DE genes from aim1 and perform Gene Set Enrichment Analysis (GSEA) to further understand the mechanisms and pathways of how these genes are involved. Specifically, we will perform GSEA to discover the potential biological pathways and mechanisms involved that contribute to the phenotypic difference in T2D patients.

We plan to use the MSigDB database, and the fgsea R package. Using the MSigDB database, we can select the gene sets relevant to T2D. Next, we will rank all the genes in our gene expression data based on their differential expressions under the effects of T2D with the fgsea function. The ranking will be performed by sorting the genes in descending order based on their absolute log2FC or effect size, such that the most upregulated genes in the first condition and the most downregulated genes in the second condition are at the top of the list. With the fgsea function, we can compute an enrichment score that reflects the degree to which the gene set is over-represented at the top or bottom of the ranked gene list. The results can be visualized with an 'enrichment map', a network visualization that represents overlaps among enriched pathways, which also allows us to discover potential connections between pathways.

#### Change in Division of Labour

Since our analysis plan has changed, our division of labour has changed slightly as well. Please see the following table for our plan for new division of labour: 
|Name|New Labour Plan|
| ------- |:-----|
| Cindy | (no change) Datasets Selection, Literature Review, Data Wrangling, PCA Visualization, Linear Regression and Result Interpretation|
| Janet |Gene Set Enrichment Analysis|
| Lucy |Gene Set Enrichment Analysis + Principle Component Analysis|

## WHAT IS THE PROGRESS OF THE ANALYSES

#### Current Progress

To achieve aim 1, factor levels have been determined for T2D and BMI variables, and data wrangling was performed for the [metadata](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/MicroarrayLinearRegressionSrc.Rmd). Principal Component Analysis was performed to visualize the  similarities and differences of the overall gene expression among all samples, as well as to conduct qualitative assessment of the overall effects of BMI and T2D. Then, more detailed analysis was conducted through a linear model with the interactive design matrix \~BMI\*T2D status using the limma package, in order to obtain the significantly DE genes in BMI\>30 healthy, BMI\<30 T2D, and BMI\>30 T2D samples ([see detailed analysis](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/MicroarrayLinearRegressionSrc.md)).

#### R Packages Used & Plan to be Used

Differential Gene Expression analysis (aim1): 
1. [knitr](https://www.r-project.org/nosvn/pandoc/knitr.html#:~:text=The%20R%20package%20knitr%20is,my%20everyday%20use%20of%20Sweave) 
2. [tidyverse](https://www.tidyverse.org/packages/)
3. [GEOquery](https://bioconductor.org/packages/release/bioc/html/GEOquery.html)
4. [limma](https://bioconductor.org/packages/release/bioc/html/limma.html)
5. [ggfortify](https://cran.r-project.org/web/packages/ggfortify/index.html)

Gene Set Enrichment Analysis (aim2): 
1. [msigdbr](https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html): curated source of gene sets and functions
2. [fgsea](https://github.com/ctlab/fgsea): tool for gene set enrichment analysis 
3. [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html): annotate gene information 
4. [data.table](https://cran.r-project.org/web/packages/data.table/index.html): for handling large tables
5. [R.utils](https://cran.r-project.org/web/packages/R.utils/index.html): for handling gzipped data

#### Current Reports

Gene Expression Analysis R Markdown file Or [Microarray Linear Regression R Markdown file](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/MicroarrayLinearRegressionSrc.Rmd)

Gene Expression Analysis Markdown file Or [Microarray Linear Regression Markdown file](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/MicroarrayLinearRegressionSrc.md)

PCA Plots [1](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/MicroarrayLinearRegressionSrc_files/figure-gfm/unnamed-chunk-7-1.png), [2](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/MicroarrayLinearRegressionSrc_files/figure-gfm/unnamed-chunk-7-2.png), [3](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/MicroarrayLinearRegressionSrc_files/figure-gfm/unnamed-chunk-7-3.png)

[Table](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/ObvsNonObHealthy.RDS) for Gene Expression Analysis Result

## RESULTS

#### Current Results

Three PCA plots were generated([Plot 1](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/MicroarrayLinearRegressionSrc_files/figure-gfm/unnamed-chunk-7-1.png), [Plot 2](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/MicroarrayLinearRegressionSrc_files/figure-gfm/unnamed-chunk-7-2.png), [Plot 3](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/MicroarrayLinearRegressionSrc_files/figure-gfm/unnamed-chunk-7-3.png)). When plotting PC1 vs PC2 and PC2 vs PC3, the graphs showed that obese individuals are slightly more closely clustered than their non-obese counterparts. On the other hand, both T2D and healthy donors were quite evenly spaced out in the plots.

As shown in [MicroarrayLinearRegressionSrc.md file](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/MicroarrayLinearRegressionSrc.md), one gene probe was associated with DE in obese vs. non-obese healthy individuals, 13 gene probes were associated with DE in T2D vs. healthy individuals with BMI\<30.

#### Do Our Primary Results Support Hypothesis?

Modifying our hypothesis from the written report, sex and age factors are not taken into consideration on the effects of how the gene expressions differ between obese and non-obese T2D patients. Overall, our results showed that BMI, T2D status and their interaction terms have different effects on the DE genes, as shown in our [MicroarrayLinearRegressionSrc.md file](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/MicroarrayLinearRegressionSrc.md). Despite the observed clustering by BMI in the PCA graph, BMI was only found to be significantly associated with one differentially expressed gene and cannot substantially support our hypothesis. One biological explanation is that due to the controversial definition of obesity caused by the current lack of clinical measurements, researchers have grouped individuals with BMI\>30 to be "obese" and others as "non-obese". Many research groups have since argued that even with BMI values less than 30, these individuals can still be "_metabolically_ obese". This suggests there may be better alternatives for categorizing obesity, which could explain why the BMI factor in our analysis showed only one DE gene (4). Additionally, looking back at our dataset, this result could also be due to the unbalance in samples where 67 samples had BMI\<30 and only 10 samples had BMI\>30. A larger dataset would be needed for further exploration.

The results for differential expression analysis can be found [here](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/ObvsNonObHealthy.RDS).

#### Challenges Encountered

Since our original dataset did not have enough samples to test the hypothesis, further literature search was performed to find a suitable dataset that could better answer our [initial biological question](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/Written%20Project%20Proposal.md). During the work for aim 1, we encounted challenges like literature review, evaluating data quality, data wrangling, and missing data handling.

## References

1.  Segerstolpe, Åsa et al. "Single-Cell Transcriptome Profiling of Human Pancreatic Islets in Health and Type 2 Diabetes." Cell metabolism vol. 24, 4 (2016): 593-607. <doi:10.1016/j.cmet.2016.08.020>
2.  Taman, Mahdi et al. "Secreted Frizzled-Related Protein 4 Reduces Insulin Secretion and Is Overexpressed in Type 2 Diabetes." Cell metabolism vol. 16, 5 (2012): 625:633. doi.org/10.1016/j.cmet.2012.10.009
3.  Tang, Yunzhao et al. "Genotype-based treatment of type 2 diabetes with an α2A-adrenergic receptor antagonist." ScienceTranslationalMedicine vol. 6, 257 (2014): 257. DOI: 10.1126/scitranslmed.3009934
4.  Olaogun I, Farag M, Hamid P. The Pathophysiology of Type 2 Diabetes Mellitus in Non-obese Individuals: An Overview of the Current Understanding. Cureus. 2020 Apr 10;12(4):e7614. doi: 10.7759/cureus.7614. PMID: 32399348; PMCID: PMC7213678.