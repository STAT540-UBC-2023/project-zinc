# A Glimpse into Type II Diabetes(T2D) at the Gene Expression Level 
_Team Zinc_

### Background

Type-2 diabetes (T2D) is a common metabolic disorder characterized by an elevation of blood glucose levels due to the lack of functional insulin [1]. This is due to either the impairment of insulin secretion or functional cells developing insulin resistance [1]. Insulin is an important polypeptide hormone responsible for regulating blood glucose level and inducing glucose storage in energy-dependent cells and tissues, such as muscles and the liver [2]. The impairment of this system will lead to many health complications, including fatigue, thirst, and vision impairment [3]. Unfortunately, symptoms for early stage T2D are not very obvious, so T2D diagnosis might come several years after onset [3]. When left untreated, it can lead to more severe complications such as kidney impairment, cardiovascular diseases and blindness [4]. 

To date, there is no effective medical cure. However, many researchers have suggested cost-effective interventions to prevent serious outcomes, including physical exercise, and diet control [3]. One potential reason may be that many published papers have shown a strong association between obesity and diabetes [5,6]. Obesity is commonly defined as having a Body Mass Index (BMI) larger than 30 kg/m2 [1,7]. The major metabolic complication of obese individuals is insulin resistance, which culminates in pancreatic cell overload and failure [8], and can ultimately lead to T2D [9,10]. However, obesity is not the only factor for T2D, as there also exists a large population of non-obese T2D patients [10]. The driving force for T2D in non-obese patients is likely genetic factors, resulting in defective insulin production [11]. Therefore, it is critical to understand the pathogenesis of both conditions for therapy developments [11].

Since clinical observations show obesity as a contributing factor for some T2D patients but not the others, it raises interesting biological questions: does obesity change gene expression for T2D patients? Is there a difference in the gene expression between T2D patients who are obese and non-obese? 

### Hypothesis and Aims

We hypothesize that T2D and obesity are associated with differential gene expression; specifically, non-obese and obese T2D individuals are expected to have different gene expression due to the different mechanisms of their conditions as proposed above. 

There are two aims to this analysis. 
- Aim 1: Perform linear regression for assessing the effect of T2D and BMI on differential gene expression. 
- Aim 2: Perform gene set enrichment analysis with REACTOME, Pathway Interaction Database (PID), ImmuneSigDB, and Kyoto Encyclopedia of Genes and Genomes (KEGG) databases for differentially expressed genes. Additionally, we will use DisGeNET to study the types of diseases associated with T2D. 

### Datasets 

To answer the biological question,  the datasets from [GSE41762](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41762) were chosen. This dataset contains microarray data for RNA harvested from islets of cadaver donors and was previously described in Pubmed articles, [PMID23140642](https://pubmed.ncbi.nlm.nih.gov/23140642/) and [PMID25298321](https://pubmed.ncbi.nlm.nih.gov/25298321/)[12,13]. There are 77 samples in total: 53 samples with non-T2D and BMI < 30, 4 samples with non-T2D and BMI > 30, 14 samples with T2D and BMI < 30 and 6 samples with T2D and BMI >30, with expression data for samples 1 to 48 obtained from one experiment and 49 through 77 from another. 

### Methods 

#### [Linear Regression](https://github.com/STAT540-UBC-2023/project-zinc/tree/main/MicroarrayLinearRegression)

###### 1.1 Data Wrangling 

The Robust Multi-Array (RMA) normalized expression matrix and metadata were downloaded from [GSE41762](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41762). The variables of interest were subsetted from metadata. Upon examination, we observed 366 genes and 29 samples with missing data from sample 49 to 77 in the original dataset. However, no detailed description of the discrepancy was found from literature. Next, singular value decomposition plots were used to visualize any potential batch effects between samples gathered from the first (sample 1 to 48) and second (sample 49 to 77) experiment in this dataset. As shown in figure1 ([1A](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/MicroarrayLinearRegression/MicroarrayLinearRegressionSrc_files/figure-gfm/unnamed-chunk-7-1.png), [1B](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/MicroarrayLinearRegression/MicroarrayLinearRegressionSrc_files/figure-gfm/unnamed-chunk-7-2.png), [1C](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/MicroarrayLinearRegression/MicroarrayLinearRegressionSrc_files/figure-gfm/unnamed-chunk-7-3.png)), a moderate batch effect was observed between the two experiments. Two approaches have been considered in dealing with the batch effect: 

1. Harmonizing data from both experiments using R packages, such as, ComBat 
2. Samples from the two experiments were analyzed separately. 

Due to data limitations, the second approach, albeit more cumbersome, was deemed to reduce risk of bias induced by data merging due to the unbalanced sample sizes and other confounders under different protocols. Hence, the second approach was chosen. 

###### 1.2 Model Fitting 

The number of samples in each category of experiment 1 and 2 were shown in table 1 and 2, respectively. Since none of the sample had BMI >30, for experiment 1 (sample 49:77) a simple linear model was used with  ~T2D, while experiment 2 (sample 1:48) was fitted with an interaction model, ~T2D* BMI. The [LIMMA](https://bioconductor.org/packages/release/bioc/html/limma.html) package in R was used for linear model fitting to the gene expression data. An empirical Bayes moderation of the standard errors was performed for each model and for the calculation computing of moderated t-statistics. LIMMA was well-suited for the purpose of this analysis due to its robustness to noise and statistical power in detecting subtle differences in microarray experiments. Statistical significance was defined as FDR<0.05 by Benjamini-Hochberg Procedure. 

Table 1. Dataset for Experiment 1 - Samples 1:48

| |BMI<30 | BMI>30 |
|------|-----|---------|
|   non-T2D |  32  |  0   | 
| T2D  |  8 |  2   |

Table 2. Dataset for Experiment 2 - Samples 49:77

||BMI<30 | BMI>30 |
|------|-----|---------|
|   non-T2D |  15  |  4   | 
| T2D  |  6 |  4   |

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

First, the [REACTOME](http://www.reactome.org) pathway was used to identify biological pathways associated with T2D. REACTOME provides a biological pathway database and detailed information on molecular events involved in human biological processes [16]. The fgsea analysis with this database showed 82 associated pathways, most of which were consistent with previous literature [17]. For example, PI3K-AKT signaling pathway plays a role in glucose metabolism and insulin signaling with a p-value of 0.000021. Another significant pathway is signaling by receptor tyrosine kinases (RTKs), which involves insulin receptor substrate (IRS) pathway important for T2D development. Interestingly, the bar plot suggested that other than the cellular regulations, many immune pathways were involved, especially the interleukin system, shown in [Figure 3](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/GeneSetEnrichmentAnalysis/Gene-Set-Enrichment-Analysis_files/figure-gfm/unnamed-chunk-30-1.png). To further validate this result and remove the potential database bias, we chose another database, PID, to perform parallel functional analysis.

###### The PID Database

The pathway interaction database presents molecular interactions and events that compose of key cellular events [19]. This database is designed to capture biological knowledge at different levels of detail. Based on this fgsea analysis shown in [Figure 4](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/GeneSetEnrichmentAnalysis/Gene-Set-Enrichment-Analysis_files/figure-gfm/unnamed-chunk-24-1.png), 25 significant pathways were identified, of which 9 were associated with the interleukin system associated with IL-1, IL-2, IL-4, IL-6, IL-12, IL-23, IL-24, and IL-27. These results were consistent with those from REACTOME and indicated high association between T2D and the human immune system. To further understand this, we decided to implement the ImmuneSigDB database and identify specific immune pathways associated with T2D. 
 
###### The ImmuneSigDB Sub-collection of Immunologic Signature Gene Sets

The ImmuneSigDB sub-collection is composed of gene sets that represent immune cell types, immune-related pathways, and gene expression signatures related to the immune system [35]. This sub-collection was generated by integrating gene expression data from multiple sources, including microarray and RNA-seq data, and curating the resulting gene sets [35]. This pathway was chosen for GSEA analysis because the REACTOME and PID GSEA analyses showed a high ratio of immune related significant pathways. 

A total of 1970 (out of 4872) significant pathways were identified based on our results ([Figure 5](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/GeneSetEnrichmentAnalysis/Gene-Set-Enrichment-Analysis_files/figure-gfm/unnamed-chunk-36-1.png)). This large number suggested that many immune associated genes and cellular molecules were related to T2D. These results explained why T2D is known as a chronic metabolic disorder that affects glucose metabolism, and has been found to be associated with chronic low-grade inflammation [20]. It is thought that chronic damage of pancreatic beta cells will lead to hyperglycemia caused by insufficient insulin production. This will subsequently induce immune response impairment [18], explaining the high involvement of many immune pathways in T2D patients. 

Some other studies have found that progression of T2D can cause dysfunction of the immune responses [21]. Other studies have analyzed different immune cells specifically, and found that in patients with T2D, the levels of interleukins (e.g., IL-17) and some T-cells are higher compared to healthy patients [22]. A few studies also looked into the genes and pathways relevant to immune system in patients with T2D, and found that an increased expression of genes are found to be related to immune responses, cell adhesion molecules pathway, tumor necrosis factor superfamily, and infectious disease pathways [23, 24]. These previous publications are consistent with our GSEA analysis results, emphasizing a strong correlation between T2D and the immune system.

###### The KEGG Database 

From our results on the strong association with the immune system, we further proposed that patients with T2D might be more susceptible to other diseases due to the weakening of their immune system. Therefore, we performed analysis using the KEGG database to identify biological and disease pathways that are associated with T2D. The Kyoto Encyclopedia of Genes and Genomes (KEGG) database provides rich information about the molecular pathways and networks within humans, such as the metabolic pathways, signaling pathways, regulatory pathways, and disease associated pathways [25]. fgsea results presented 34 significant pathways associated with T2D (FDR<0.05) shown in [Figure 6](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/GeneSetEnrichmentAnalysis/Gene-Set-Enrichment-Analysis_files/figure-gfm/unnamed-chunk-20-1.png). These results not only showed high involvement of many insulin-related signaling pathways, such as insulin signaling pathways, and pancreatic cancer pathways, they also presented the involvement of many disease pathways. For example, melanoma, prostate cancer, small cell lung cancer, and colorectal cancer. Many other pathways, such as MAPK signaling pathway and cytokine-cytokine receptor interaction, also showed high significance, all of which are involved in insulin pathways. The dysfunction of these pathways might cause the failure of insulin usage, a hallmark for T2D. 

###### Results from REACTOME and KEGG Gene Sets 

Combining the results from REACTOME, PID and KEGG gene sets, we have identified more than 100 significant pathways associated with the risk of T2D. From the first glance, many of the pathways identified are related with cell growth regulation pathways. For example, the processing of capped intron containing per mRNA - a pathway related to cellular processing during transcription; spliceosome pathway - involved in transcription regulation [26]. They provided strong indications that T2D had an effect on cellular growth regulation. These genomic and pathway enrichment evidence suggested the converge on vital pathways including insulin signaling, JAK-STAT signaling, P13-AKT signaling, and TGF-beta signaling, and these vital pathway were involved to play a critical role in pancreatic islets maturity and function, and insulin secretion.
Furthermore, we found that T2D-gene signatures identified in enrichment analysis can elucidate disease conditions when the interlinked genes are taken together with their pathways and functional level interactors. Therefore, we selected KEGG with the most relevant results to disease and plotted them with an enrichment map ([Figure 7](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/GeneSetEnrichmentAnalysis/Gene-Set-Enrichment-Analysis_files/figure-gfm/unnamed-chunk-21-1.png)) to view gene sets as a network and a ridge plot ([Figure 8](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/GeneSetEnrichmentAnalysis/Gene-Set-Enrichment-Analysis_files/figure-gfm/unnamed-chunk-21-2.png)) to identify the enrichment distribution of the significant pathways. As shown in the enrichment map, the KEGG enrichment map showed a network between pancreatic cancer pathway and colorectal cancer pathway, renal cell carcinoma pathway, EGFR tyrosine kinase inhibitor resistance, and melanoma pathway. This result was also confirmed with the ridge plot, where all pathways involving the above diseases were upregulated with an increase in log fold change. In particular, colorectal cancer (CRC) and renal cell carcinoma (RCC) have a close relationship with diabetes [27,28]. Although many studies have shown the association between T2D and colorectal cancer, the exact mechanism is unknown [27]. Colorectal cancer refers to when the cells in the colon grow out of control [29]. One study portrayed the association between a common genetic variance on TCF7L2 and CRC [34]. Previously, we obtained a list of 1693 genes associated with Type II Diabetes using GWAS database. TCF7L2 contained a strong association with T2D (p = 5e-18), which indicated the potential molecular connections between CRC and T2D. Many other significant pathways identified in the previous fgsea were also shared with colorectal cancer, such as the MAPK pathway and the TGF beta signaling pathway.
 
Another interesting disease identified is renal cell carcinoma, which also has been shown to have a strong association with T2D [28]. Previous literature has shown that insulin plays a significant role in RCC development [31]. A commonly known symptom of T2D is hyperinsulinemia, which can lead to an activation of downstream signaling pathways, such as PI3K pathway [30]. This was also verified by the results obtained from REACTOME fgsea analysis. As PI3K pathway is involved in the formation of renal cell carcinoma [32], this shared molecular pathway might be an indication of the associations between these two diseases and it might potentially shine a light on the therapeutics targeting patients with T2D to prevent further progression of its associated diseases, such as renal cell carcinoma. Additionally, others, such as JAK-STAT pathways and MAPK pathways, are also involved in both T2D and renal cell carcinoma regulations [32]. 

The ridgeline plot from KEGG was also generated to further validate the results. As expected, insulin secretion and the beta cell gene expression were downregulated, causing diabetes. Interestingly, many metabolism pathways are also downregulated for T2D patients, such as glutathione metabolism and propanoate metabolism. As T2D is mainly characterized by insulin resistance leading to impairment of glucose uptake and usage and high blood sugar level. This will result in an overproduction of insulin, hyperinsulinemia, which results in many complicated metabolism impairment. The ridge plots also identified branched-chain amino acid catabolism as downregulated. This may result in a build up of branched amino acids in the body, which is associated with insulin resistance - commonly associated with T2D [33]. 

###### Other Disease Associations 

Lastly, [DisGeNET](https://www.disgenet.org/), a large collection of genes and variants associated with human diseases [35], was used to identify diseases associated with T2D based on the DE genes shown in [Figure 9](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/GeneSetEnrichmentAnalysis/Gene-Set-Enrichment-Analysis_files/figure-gfm/unnamed-chunk-38-1.png). As expected, different immune lymphoma diseases and inflammations were identified as having correlations with T2D. Interestingly, about half of the diseases from the top 20 list were exclusively female specific, such as invasive carcinoma of breast and endometrioma. Since there is no specific evidence showing a stronger association between T2D and females from previous studies, this led us to think that T2D might have a stronger effect on female health status if left untreated. 

### Conclusion/Significance

Our linear regression analysis identified 67 differential expressed genes that were associated with T2D. Overall, the GSEA study identified more than 100 unique signaling and metabolic pathways and 1970 immune pathways that were associated with T2D disease. These results signify the importance of genetic variation and their functional roles in Type-2-diabetes. The knowledge of these genetic factors can be used to elucidate the T2D disease mechanisms and better understand its association with other diseases. Particularly, we identified a high connection between T2D and colorectal cancer and melanoma. There is also a high association between T2D and many female focused diseases/symptoms. Collectively, these genetic markers can also be used to provide insights for therapeutic development and bring invaluable medical treatments to those who are in need. 


### References. 

1. Galicia-Garcia, Unai et al. “Pathophysiology of Type 2 Diabetes Mellitus.” _International journal of molecular sciences_ vol. 21,17 6275. 30 Aug. 2020, doi:10.3390/ijms21176275
2. Rahman MS, Hossain KS, Das S, Kundu S, Adegoke EO, Rahman MA, Hannan MA, Uddin MJ, Pang MG. Role of Insulin in Health and Disease: An Update. Int J Mol Sci. 2021 Jun 15;22(12):6403. doi: 10.3390/ijms22126403. PMID: 34203830; PMCID: PMC8232639.
3. World Health Organizations. Diabetes. 2023.
4. Deshpande, Anjali D et al. “Epidemiology of diabetes and diabetes-related complications.” Physical therapy vol. 88,11 (2008): 1254-64. doi:10.2522/ptj.20080020
5. Conway, Baqiyyah et al. The obesity epidemic and rising diabetes incidence in a low-income racially diverse southern US cohort. PLOS ONE. 2018. https://doi.org/10.1371/journal.pone.0190993
6. Magkos, F., Hjorth, M.F. & Astrup, A. Diet and exercise in the prevention and treatment of type 2 diabetes mellitus. Nat Rev Endocrinol 16, 545–555 (2020). https://doi.org/10.1038/s41574-020-0381-5
7. Al-Goblan, Abdullah S et al. “Mechanism linking diabetes mellitus and obesity.” Diabetes, metabolic syndrome and obesity : targets and therapy vol. 7 587-91. 4 Dec. 2014, doi:10.2147/DMSO.S67400
8. Singla, Parul et al. “Metabolic effects of obesity: A review.” _World journal of diabetes_ vol. 1,3 (2010): 76-88. doi:10.4239/wjd.v1.i3.76
9. Prentki, Marc, and Christopher J Nolan. “Islet beta cell failure in type 2 diabetes.” _The Journal of clinical investigation_ vol. 116,7 (2006): 1802-12. doi:10.1172/JCI29103
10. Olaogun, Idowu et al. “The Pathophysiology of Type 2 Diabetes Mellitus in Non-obese Individuals: An Overview of the Current Understanding.” _Cureus_ vol. 12,4 e7614. 10 Apr. 2020, doi:10.7759/cureus.7614
11. Vaag, Allan, and Søren S Lund. “Non-obese patients with type 2 diabetes and prediabetic subjects: distinct phenotypes requiring special diabetes treatment and (or) prevention?.” _Applied physiology, nutrition, and metabolism = Physiologie appliquee, nutrition et metabolisme_ vol. 32,5 (2007): 912-20. doi:10.1139/H07-100
12. Taman, Mahdi et al. "Secreted Frizzled-Related Protein 4 Reduces Insulin Secretion and Is Overexpressed in Type 2 Diabetes." Cell metabolism vol. 16, 5 (2012): 625:633. doi.org/10.1016/j.cmet.2012.10.009
13. Tang, Yunzhao et al. "Genotype-based treatment of type 2 diabetes with an α2A-adrenergic receptor antagonist." ScienceTranslationalMedicine vol. 6, 257 (2014): 257. DOI: 10.1126/scitranslmed.3009934
14. Subramanian, A. et al. (2005) “Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles,” Proceedings of the National Academy of Sciences, 102(43), pp. 15545–15550. Available at: https://doi.org/10.1073/pnas.0506580102.
15. Olaogun I, Farag M, Hamid P. The Pathophysiology of Type 2 Diabetes Mellitus in Non-obese Individuals: An Overview of the Current Understanding. Cureus. 2020 Apr 10;12(4):e7614. doi: 10.7759/cureus.7614. PMID: 32399348; PMCID: PMC7213678.
16. Fabregat A, Jupe S, Matthews L, et al. The Reactome Pathway Knowledgebase. Nucleic Acids Res. 2018 Jan 4;46(D1):D649-D655. doi: 10.1093/nar/gkx1132. PMID: 29145629.
17. Muhammad, Syed, et al. Cellular Signaling Pathways in Insulin Resistance-Systems Biology Analyses of Microarray Dataset Reveals New Drug Target Gene Signatures of Type 2 Diabetes Mellitus. Sec. Systems Biology Archive. 2017; 8. https://doi.org/10.3389/fphys.2017.00013
18. Berbudi, Afiat et al. “Type 2 Diabetes and its Impact on the Immune System.” Current diabetes reviews vol. 16,5 (2020): 442-449. doi:10.2174/1573399815666191024085838
19. Schaefer, Carl F et al. “PID: the Pathway Interaction Database.” Nucleic acids research vol. 37,Database issue (2009): D674-9. doi:10.1093/nar/gkn653
20. Ferlita S, Yegiazaryan A, Noori N, Lal G, Nguyen T, To K, Venketaraman V. Type 2 Diabetes Mellitus and Altered Immune System Leading to Susceptibility to Pathogens, Especially Mycobacterium tuberculosis. Journal of Clinical Medicine. 2019; 8(12):2219. https://doi.org/10.3390/jcm8122219
21. Berbudi, Afiat et al. “Type 2 Diabetes and its Impact on the Immune System.” Current diabetes reviews vol. 16,5 (2020): 442-449. doi:10.2174/1573399815666191024085838
22. Francisco, C O et al. “Cytokine profile and lymphocyte subsets in type 2 diabetes.” Brazilian journal of medical and biological research = Revista brasileira de pesquisas medicas e biologicas vol. 49,4 (2016): e5062. doi:10.1590/1414-431X20155062
23. Wu, Chun et al. “Whole-genome expression analyses of type 2 diabetes in human skin reveal altered immune function and burden of infection.” Oncotarget vol. 8,21 (2017): 34601-34609. doi:10.18632/oncotarget.16118
24. Tonyan, Ziravard N et al. “Overview of Transcriptomic Research on Type 2 Diabetes: Challenges and Perspectives.” Genes vol. 13,7 1176. 30 Jun. 2022, doi:10.3390/genes13071176
25. Kanehisa M, Furumichi M, Sato Y, Ishiguro-Watanabe M, Tanabe M. KEGG: integrating viruses and cellular organisms. Nucleic Acids Res. 2021 Jan 8;49(D1):D545-D551. doi: 10.1093/nar/gkaa970. PMID: 33137290.
26. Matera, A., Wang, Z. A day in the life of the spliceosome. Nat Rev Mol Cell Biol 15, 108–121 (2014). https://doi.org/10.1038/nrm3742
27. Yao, Caroline et al. “Management of colorectal cancer and diabetes.” Journal of the Royal Society of Medicine vol. 107,3 (2014): 103-9. doi:10.1177/0141076813512121
28. Tseng, Chin-Hsiao. “Type 2 Diabetes Mellitus and Kidney Cancer Risk: A Retrospective Cohort Analysis of the National Health Insurance.” PloS one vol. 10,11 e0142480. 11 Nov. 2015, doi:10.1371/journal.pone.0142480
29. Malki, Ahmed et al. “Molecular Mechanisms of Colon Cancer Progression and Metastasis: Recent Insights and Advancements.” International journal of molecular sciences vol. 22,1 130. 24 Dec. 2020, doi:10.3390/ijms22010130
30. Hopkins, Benjamin D et al. “Insulin-PI3K signalling: an evolutionarily insulated metabolic driver of cancer.” Nature reviews. Endocrinology vol. 16,5 (2020): 276-283. doi:10.1038/s41574-020-0329-9
31. Solarek W, Czarnecka AM, Escudier B, Bielecka ZF, Lian F, Szczylik C. Insulin and IGFs in renal cancer risk and progression. Endocr Relat Cancer. 2015 Oct;22(5):R253-64. doi: 10.1530/ERC-15-0135. PMID: 26330483.
32. Guo, Huifang et al. “The PI3K/AKT Pathway and Renal Cell Carcinoma.” Journal of genetics and genomics = Yi chuan xue bao vol. 42,7 (2015): 343-53. doi:10.1016/j.jgg.2015.03.003
33. Cuomo, Paola et al. “Role of Branched-Chain Amino Acid Metabolism in Type 2 Diabetes, Obesity, Cardiovascular Disease and Non-Alcoholic Fatty Liver Disease.” International journal of molecular sciences vol. 23,8 4325. 13 Apr. 2022, doi:10.3390/ijms23084325
34. Wenzel, J., Rose, K., Haghighi, E.B. et al. Loss of the nuclear Wnt pathway effector TCF7L2 promotes migration and invasion of human colorectal cancer cells. Oncogene 39, 3893–3909 (2020). https://doi.org/10.1038/s41388-020-1259-7
35. Godec J, Tan Y, Liberzon A, Tamayo P, Bhattacharya S, Butte A, Mesirov JP, Haining WN, Compendium of Immune Signatures Identifies Conserved and Species-Specific Biology in Response to Inflammation, 2016, Immunity 44(1), 194-206.
