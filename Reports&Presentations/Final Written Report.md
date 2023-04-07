# A Glimpse into Type II Diabetes(T2D) at the Gene Expression Level 
_Zinc Group_

###Background
Type-2 diabetes (T2D) is a common metabolic disorder characterized by an elevation of blood glucose levels due to the lack of functional insulin [1]. This is due to either the impairment of insulin secretion or the development of insulin resistance by functional cells [1]. Insulin is an important polypeptide hormone responsible for regulating blood glucose level and inducing glucose storage in energy-dependent cells and tissues, such as muscles and liver [2]. The impairment of this system will lead to many symptoms, including fatigue, thirst, and vision loss [3]. Unfortunately, symptoms for early stage T2D are not very obvious, so T2D diagnosis might take several years after onset [3]. When left untreated, it can lead to complications such as kidney impairment, cardiovascular diseases and blindness [4]. 

To date, there is no effective cure. However, many researchers have suggested cost-effective intervention to prevent serious outcome, including physical exercise, and diet control [3]. One of the reasons is because they increasingly noticed a strong association between obesity and diabetes, shown in many published papers [5,6]. Many patients with T2D are commonly portrayed as having accumulated higher fat percentages in their bodies [1]. Obesity is commonly defined as Body Mass Index (BMI) larger than 30 kg/m2 [1,7]. The major metabolic complications of obese individuals consist of insulin resistance, culminating in pancreatic cell overload and failure [8], which leads to T2D [9,10]. However, obesity is not the only factor for T2D, as a large population of non-obese T2D patients has also been rising [10]. The driving force for T2D in non-obese patients leans towards a genetic factor, resulting in defective insulin production [11]. Therefore, it is critical to understand the pathogenesis of both conditions for therapies developments [11].

Since clinical observation shows obesity as a contributing factor for some T2D patience but not the others, it brings up interesting biological questions: does obesity contribute to gene expression change for T2D patients? If so, what is the differential gene expression between obese and non-obese T2D individuals? 

###Hypothesis and Aims

We hypothesize that T2D and obesity are associated with differential gene expression. Moreover, the non-obese and obese T2D individuals are expected to have different gene expression due to the different mechanisms proposed above. 

There are two aims to this analysis. 
+ Aim1: perform linear regression for assessing the effect of T2D and BMI on differential gene expression. 
+ Aim2: perform gene set enrichment analysis with Reactome, Pathway Interaction Database (PID), ImmuneSigDB, and Kyoto Encyclopedia of Genes and Genomes (KEGG) databases for differentially expressed genes. Additionally, we will use DisGeNET to study the types of diseases associated with T2D. 

###Datasets 

To answer the biological question,  the datasets from [GSE41762](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41762) were chosen. This dataset contains microarray data for RNA harvested from islets of cadaver donors and was previously described in Pubmed articles, [PMID23140642](https://pubmed.ncbi.nlm.nih.gov/23140642/) and [PMID25298321](https://pubmed.ncbi.nlm.nih.gov/25298321/)[12,13]. There are 77 samples in total - 53 samples with non-T2D and BMI < 30, 4 samples with non-T2D and BMI > 30, 14 samples with T2D and BMI < 30 and 6 samples with T2D and BMI >30, with expression data for samples 1 to 48 obtained from one experiment and 49 through 77 from another. 

###Methods 

####[Linear Regression](https://github.com/STAT540-UBC-2023/project-zinc/tree/main/MicroarrayLinearRegression)




