Gene-Set-Enrichment-Analysis
================
Lucy Chi & Janet Zhang
2023-03-22

## Gene Sets for Downstream Gene Set Enrichment Analysis

### Obtain Molecular Signature Database (KEGG, reactome and GO pathways)

``` r
knitr::kable(msigdbr::msigdbr_species())
```

| species_name                    | species_common_name                                                     |
|:--------------------------------|:------------------------------------------------------------------------|
| Anolis carolinensis             | Carolina anole, green anole                                             |
| Bos taurus                      | bovine, cattle, cow, dairy cow, domestic cattle, domestic cow, ox, oxen |
| Caenorhabditis elegans          | NA                                                                      |
| Canis lupus familiaris          | dog, dogs                                                               |
| Danio rerio                     | leopard danio, zebra danio, zebra fish, zebrafish                       |
| Drosophila melanogaster         | fruit fly                                                               |
| Equus caballus                  | domestic horse, equine, horse                                           |
| Felis catus                     | cat, cats, domestic cat                                                 |
| Gallus gallus                   | bantam, chicken, chickens, Gallus domesticus                            |
| Homo sapiens                    | human                                                                   |
| Macaca mulatta                  | rhesus macaque, rhesus macaques, Rhesus monkey, rhesus monkeys          |
| Monodelphis domestica           | gray short-tailed opossum                                               |
| Mus musculus                    | house mouse, mouse                                                      |
| Ornithorhynchus anatinus        | duck-billed platypus, duckbill platypus, platypus                       |
| Pan troglodytes                 | chimpanzee                                                              |
| Rattus norvegicus               | brown rat, Norway rat, rat, rats                                        |
| Saccharomyces cerevisiae        | baker’s yeast, brewer’s yeast, S. cerevisiae                            |
| Schizosaccharomyces pombe 972h- | NA                                                                      |
| Sus scrofa                      | pig, pigs, swine, wild boar                                             |
| Xenopus tropicalis              | tropical clawed frog, western clawed frog                               |

``` r
#Access the REACTOME gene sets. 
kegg.human.db <- msigdbr::msigdbr(species = "human",
                                  category = "C2",
                                  subcategory = "CP:KEGG")
```

``` r
#Access the REACTOME gene sets. 
#A standardized vocabulary of terms that describe gene function, cellular components, and biological processes.
reactome.human.db <- msigdbr(species = "human", 
                            category = "C2", 
                            subcategory = "CP:REACTOME")
```

``` r
#Access the GO gene sets. 
#Using all GO gene sets: biological process, cellular component, and molecular function. 
GO.human.db <- msigdbr(species = "human", 
                            category = "C5")
```

### Obtain GWAS Catalog Information

``` r
run.if.needed <- function(.file, .code) {
    if(!file.exists(.file)) { .code }
    stopifnot(file.exists(.file))
}
```

``` r
gwas.tidy.file <- "../Files/gwas_catalog_tidy.tsv.gz"
run.if.needed(gwas.tidy.file, {

    gwas.file <- "../Files/gwas_catalog_v1.0-associations_e105_r2022-02-02.tsv.gz"

    run.if.needed(gwas.file, {
        url <- "https://www.ebi.ac.uk/gwas/api/search/downloads/full"
        .file <- str_remove(gwas.file, ".gz$")
        download.file(url, destfile = .file)
        gzip(.file)
        unlink(.file)
    })

    .dt <-
        fread(gwas.file, sep="\t", quote="") %>%
        dplyr::select(`MAPPED_GENE`, `DISEASE/TRAIT`, `PVALUE_MLOG`)

    .dt <- .dt[order(.dt$PVALUE_MLOG, decreasing = TRUE),
               head(.SD, 1),
               by = .(`MAPPED_GENE`, `DISEASE/TRAIT`)]

    .count <- .dt[, .(.N), by = .(`DISEASE/TRAIT`)]
    .dt <- left_join(.count[`N` >= 100, ], .dt)[nchar(`MAPPED_GENE`)> 0,]

    .dt <- .dt[,
               .(gene_symbol = unlist(strsplit(`MAPPED_GENE`, split="[ ,.-]+"))),
               by = .(`DISEASE/TRAIT`, PVALUE_MLOG)]
    .dt[, p.value := 10^(-PVALUE_MLOG)]

    fwrite(.dt, file=gwas.tidy.file)
})
```

### Genes associated with T2D disease

``` r
# set seed.
set.seed(10)

#' Convert a number into a scientifically-formatted string
#' @param x number
num.sci <- function(x) {
    format(x, scientific=TRUE, digits = 2)
}

gwas.db <- fread(gwas.tidy.file)    
gwas.db[, gs_name := `DISEASE/TRAIT`]

t2d_genes<- gwas.db[str_detect(`DISEASE/TRAIT`, "[Dd]iabetes") & !is.na(gene_symbol)] %>%
    mutate(p.value = num.sci(p.value)) %>%
    dplyr::select(`gs_name`, `gene_symbol`, `p.value`, `PVALUE_MLOG`) 

t2d_genes <- t2d_genes %>% filter(t2d_genes$gs_name == "Type 2 diabetes")
## 2364 genes identified
```

### Load DE Gene List

``` r
# read data frame for the DE genes
DEgenes<-readRDS("../Tables/degT2d_DE.RDS")
```

### Data formatting for `fgsea` analysis

``` r
# matched Ensembl IDs to probeIDs
DEgenes$affy_hugene_1_0_st_v1 <- row.names(DEgenes) 
head(DEgenes)
```

    ##             logFC  AveExpr        t      P.Value   adj.P.Val        B
    ## 8003667 1.4034296 7.609690 6.316668 6.816466e-08 0.001958371 7.545931
    ## 8095080 1.0921294 7.404872 5.542533 1.090848e-06 0.007738258 5.146712
    ## 8043995 0.8684121 8.800716 5.529337 1.143097e-06 0.007738258 5.106089
    ## 7938608 1.2241985 7.428053 5.511468 1.217814e-06 0.007738258 5.051107
    ## 7952341 1.0060584 7.362819 5.435595 1.592742e-06 0.007738258 4.817966
    ## 8139087 1.8199810 4.825679 5.431479 1.616065e-06 0.007738258 4.805336
    ##         affy_hugene_1_0_st_v1
    ## 8003667               8003667
    ## 8095080               8095080
    ## 8043995               8043995
    ## 7938608               7938608
    ## 7952341               7952341
    ## 8139087               8139087

``` r
# convert transcript cluster IDs (also desribed as probeIDs) into other gene ids
ensembl_human = useMart("ensembl",dataset="hsapiens_gene_ensembl")
ensemble_genes <- getBM(attributes=c("ensembl_gene_id", "affy_hugene_1_0_st_v1", "hgnc_symbol", "entrezgene_id"),
      filters = "affy_hugene_1_0_st_v1",
      values = rownames(DEgenes),
      useCache=FALSE,
      mart=ensembl_human)
```

``` r
# data clean up and rename
ensemble_genes_clean <-
    ensemble_genes %>%  
    dplyr::rename(gene_symbol = hgnc_symbol) %>% 
    as.data.table()

ensemble_genes_clean$affy_hugene_1_0_st_v1 <- as.character(ensemble_genes_clean$affy_hugene_1_0_st_v1)

# add on the p-values & other information in DEgenes
ensembl_geneIDs <- ensemble_genes_clean %>% left_join(DEgenes, by = "affy_hugene_1_0_st_v1") 

saveRDS(ensembl_geneIDs, file = "../Files/ensembl_geneIDs.RDS")

# Show NA's & other stats
summary(ensemble_genes$entrezgene_id)
```

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
    ##         1     11166    146318  37405545 124900242 125467750     14185

``` r
#filter the list with logFC values for KEGG analysis later on
logFC_genes <- ensembl_geneIDs$logFC
# name the vector
names(logFC_genes) <- ensembl_geneIDs$ensembl_gene_id
# omit any NA values 
gene_list<-na.omit(logFC_genes)
# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
```

``` r
# Create a new dataframe sorted_genes which has only the filtered genes using log(adj.P)
ensembl_geneIDs <- ensembl_geneIDs %>% mutate(log_adjP = -log10(adj.P.Val))

sorted_genes <-
    ensembl_geneIDs[order(ensembl_geneIDs$log_adjP, decreasing = TRUE),
                   head(.SD, 1),
                   by = .(gene_symbol)]

# Create a vector of the gene
deg.ranks <- sorted_genes$log_adjP
# Name vector with ENTREZ ids
names(deg.ranks) <- sorted_genes$gene_symbol
# omit any NA values 
deg.ranks<-na.omit(deg.ranks)
# sort the list in decreasing order (required for clusterProfiler)
deg.ranks = sort(deg.ranks, decreasing = TRUE)
```

## Rank-based Gene Set Enrichment Analysis

### `fgsea` analysis set-up

``` r
# Make a function that can help us to prepare the gene list to feed into the fgsea arguments. 
make.gs.lol <- function(data) {
    data <- as.data.table(data) %>% unique()
    data_list <-
        data[, .(gene = .(gene_symbol)), by = .(gs_name)] %>%
        as.list()
    names <- data_list$gs_name
    return_value <- data_list$gene
    names(return_value) <- names
    return(return_value)
}
```

### KEGG Pathway Analysis

``` r
#Convert the KEGG pathway tibble
kegg.lol <- kegg.human.db %>% dplyr::select(gene_symbol, gs_name) %>% make.gs.lol()
kegg.fgsea <- fgsea::fgsea(pathways = kegg.lol, stats = deg.ranks, scoreType = "pos")
```

    ## Warning in preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam, : There are ties in the preranked stats (84.06% of the list).
    ## The order of those tied genes will be arbitrary, which may produce unexpected results.

### KEGG `fgsea` Analysis

``` r
kegg.fgsea[,
           topGenes := paste0(head(unlist(`leadingEdge`), 5), collapse=", "),
           by = .(pathway)]

kegg.fgsea %>%
    arrange(pval) %>%
    head(10) %>% 
    dplyr::select(-leadingEdge)
```

    ##                                         pathway         pval         padj
    ##  1: KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION 9.407489e-11 1.749793e-08
    ##  2:                     KEGG_PATHWAYS_IN_CANCER 1.804212e-08 1.677917e-06
    ##  3:             KEGG_JAK_STAT_SIGNALING_PATHWAY 2.843616e-08 1.763042e-06
    ##  4:                 KEGG_MAPK_SIGNALING_PATHWAY 3.009503e-06 1.399419e-04
    ##  5:             KEGG_HEMATOPOIETIC_CELL_LINEAGE 8.961631e-06 3.333727e-04
    ##  6:                    KEGG_ALLOGRAFT_REJECTION 1.940800e-05 6.016481e-04
    ##  7:                 KEGG_PPAR_SIGNALING_PATHWAY 7.501648e-05 1.993295e-03
    ##  8:                      KEGG_VIRAL_MYOCARDITIS 8.968111e-05 2.063227e-03
    ##  9:             KEGG_TGF_BETA_SIGNALING_PATHWAY 9.983356e-05 2.063227e-03
    ## 10:                               KEGG_MELANOMA 1.666445e-04 3.047829e-03
    ##       log2err        ES      NES size                                topGenes
    ##  1: 0.8390889 0.6967663 1.350655  255         PDGFRA, IL1R1, IL18R1, IL6, HGF
    ##  2: 0.7337620 0.6591685 1.279391  300            PDGFRA, FGF7, FGF2, IL6, HGF
    ##  3: 0.7337620 0.7153785 1.376867  151      SPRY1, SOCS2, IL6, IL13RA2, SPRED1
    ##  4: 0.6272567 0.6548142 1.269122  252        PDGFRA, IL1R1, FGF7, FGF2, DUSP5
    ##  5: 0.5933255 0.7307650 1.392365   84     IL1R1, IL6, KIT, HLA-DRB3, HLA-DRB1
    ##  6: 0.5756103 0.8305163 1.533336   31 CD86, FAS, HLA-DRB3, HLA-DRB1, HLA-DQA2
    ##  7: 0.5384341 0.7483599 1.411618   64           ACSL4, LPL, ACSL5, ME1, PPARD
    ##  8: 0.5384341 0.7497239 1.407213   58    CD86, HLA-DRB3, HLA-DRB1, ABL1, CD55
    ##  9: 0.5384341 0.7069601 1.347009   84          DCN, BMPR1B, THBS2, MYC, LTBP1
    ## 10: 0.5188481 0.7193180 1.359351   66         PDGFRA, FGF7, FGF2, HGF, CDKN1A

``` r
# filter significant pathways. Set p<0.05. 
kegg_significant_pathways <- kegg.fgsea[padj < 0.05]
kegg_significant_pathways
```

    ##                                               pathway         pval         padj
    ##  1:                             KEGG_ABC_TRANSPORTERS 1.224148e-03 1.264953e-02
    ##  2:                       KEGG_ACUTE_MYELOID_LEUKEMIA 3.089347e-03 2.179191e-02
    ##  3:                            KEGG_ADHERENS_JUNCTION 5.837401e-03 3.619189e-02
    ##  4:                          KEGG_ALLOGRAFT_REJECTION 1.940800e-05 6.016481e-04
    ##  5:          KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION 1.334596e-03 1.275413e-02
    ##  6:                                    KEGG_APOPTOSIS 2.275399e-03 1.840105e-02
    ##  7:                                       KEGG_ASTHMA 7.547702e-03 4.374318e-02
    ##  8:                   KEGG_AUTOIMMUNE_THYROID_DISEASE 1.665939e-03 1.475546e-02
    ##  9:                            KEGG_COLORECTAL_CANCER 7.996066e-03 4.374318e-02
    ## 10:          KEGG_COMPLEMENT_AND_COAGULATION_CASCADES 7.922413e-04 9.209805e-03
    ## 11:       KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION 9.407489e-11 1.749793e-08
    ## 12:                     KEGG_ECM_RECEPTOR_INTERACTION 4.275828e-03 2.840372e-02
    ## 13:                           KEGG_ENDOMETRIAL_CANCER 6.209205e-03 3.725523e-02
    ## 14:                    KEGG_GRAFT_VERSUS_HOST_DISEASE 7.189713e-04 8.915244e-03
    ## 15:                   KEGG_HEMATOPOIETIC_CELL_LINEAGE 8.961631e-06 3.333727e-04
    ## 16: KEGG_INTESTINAL_IMMUNE_NETWORK_FOR_IGA_PRODUCTION 9.480305e-04 1.037257e-02
    ## 17:                   KEGG_JAK_STAT_SIGNALING_PATHWAY 2.843616e-08 1.763042e-06
    ## 18:                       KEGG_MAPK_SIGNALING_PATHWAY 3.009503e-06 1.399419e-04
    ## 19:         KEGG_MATURITY_ONSET_DIABETES_OF_THE_YOUNG 2.275399e-03 1.840105e-02
    ## 20:                                     KEGG_MELANOMA 1.666445e-04 3.047829e-03
    ## 21:    KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY 7.846611e-03 4.374318e-02
    ## 22:                           KEGG_PATHWAYS_IN_CANCER 1.804212e-08 1.677917e-06
    ## 23:                       KEGG_PPAR_SIGNALING_PATHWAY 7.501648e-05 1.993295e-03
    ## 24:                              KEGG_PROSTATE_CANCER 1.802479e-04 3.047829e-03
    ## 25:                            KEGG_PURINE_METABOLISM 2.423389e-03 1.878126e-02
    ## 26:                        KEGG_PYRIMIDINE_METABOLISM 5.168154e-03 3.314747e-02
    ## 27:             KEGG_REGULATION_OF_ACTIN_CYTOSKELETON 3.163342e-03 2.179191e-02
    ## 28:                       KEGG_SMALL_CELL_LUNG_CANCER 1.371412e-03 1.275413e-02
    ## 29:                                  KEGG_SPLICEOSOME 5.449555e-04 7.797055e-03
    ## 30:                 KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS 6.548601e-04 8.700284e-03
    ## 31:                   KEGG_TGF_BETA_SIGNALING_PATHWAY 9.983356e-05 2.063227e-03
    ## 32:                    KEGG_TYPE_II_DIABETES_MELLITUS 3.622760e-04 5.615279e-03
    ## 33:                     KEGG_TYPE_I_DIABETES_MELLITUS 2.978354e-03 2.179191e-02
    ## 34:                            KEGG_VIRAL_MYOCARDITIS 8.968111e-05 2.063227e-03
    ##                                               pathway         pval         padj
    ##       log2err        ES      NES size
    ##  1: 0.4550599 0.7505785 1.395799   39
    ##  2: 0.4317077 0.6991183 1.309008   52
    ##  3: 0.4070179 0.6807228 1.285649   65
    ##  4: 0.5756103 0.8305163 1.533336   31
    ##  5: 0.4550599 0.6824599 1.296661   76
    ##  6: 0.4317077 0.6663575 1.268547   80
    ##  7: 0.4070179 0.7513441 1.379633   26
    ##  8: 0.4550599 0.7218970 1.347642   46
    ##  9: 0.3807304 0.6837331 1.283305   55
    ## 10: 0.4772708 0.6988542 1.321510   67
    ## 11: 0.8390889 0.6967663 1.350655  255
    ## 12: 0.4070179 0.6595044 1.254538   78
    ## 13: 0.4070179 0.6970367 1.298790   47
    ## 14: 0.4772708 0.7809783 1.447523   34
    ## 15: 0.5933255 0.7307650 1.392365   84
    ## 16: 0.4772708 0.7489867 1.395549   42
    ## 17: 0.7337620 0.7153785 1.376867  151
    ## 18: 0.6272567 0.6548142 1.269122  252
    ## 19: 0.4317077 0.7939079 1.443017   23
    ## 20: 0.5188481 0.7193180 1.359351   66
    ## 21: 0.3807304 0.6206217 1.188968  131
    ## 22: 0.7337620 0.6591685 1.279391  300
    ## 23: 0.5384341 0.7483599 1.411618   64
    ## 24: 0.5188481 0.6988932 1.330127   79
    ## 25: 0.4317077 0.6265446 1.206178  152
    ## 26: 0.4070179 0.6477453 1.236896   94
    ## 27: 0.4317077 0.6118705 1.180229  202
    ## 28: 0.4550599 0.6809859 1.295401   78
    ## 29: 0.4772708 0.6533644 1.249572  125
    ## 30: 0.4772708 0.6557024 1.255220  128
    ## 31: 0.5384341 0.7069601 1.347009   84
    ## 32: 0.4984931 0.7634488 1.422495   42
    ## 33: 0.4317077 0.7483863 1.387202   35
    ## 34: 0.5384341 0.7497239 1.407213   58
    ##       log2err        ES      NES size
    ##                                           leadingEdge
    ##  1:           ABCA6,ABCA8,ABCA9,ABCC8,ABCC6,ABCA7,...
    ##  2:                KIT,MYC,STAT5A,RELA,PPARD,TCF7,...
    ##  3:          NLK,CSNK2A1,TGFBR2,YES1,NECTIN3,TCF7,...
    ##  4:     CD86,FAS,HLA-DRB3,HLA-DRB1,HLA-DQA2,HLA-E,...
    ##  5:    HLA-DRB3,HLA-DRB1,CTSL,PSME3,LGMN,HLA-DQA2,...
    ##  6:     IL1R1,IRAK2,IRAK3,TNFRSF10D,FAS,TNFRSF10A,...
    ##  7: HLA-DRB3,HLA-DRB1,HLA-DQA2,CD40,HLA-DQA1,IL10,...
    ##  8:     CD86,FAS,HLA-DRB3,HLA-DRB1,HLA-DQA2,HLA-E,...
    ##  9:              MYC,TGFBR2,TCF7,PIK3R1,RAF1,RAC1,...
    ## 10:                 A2M,C1S,TFPI,C1R,BDKRB1,PLAUR,...
    ## 11:               PDGFRA,IL1R1,IL18R1,IL6,HGF,KIT,...
    ## 12:           THBS2,THBS4,LAMC2,LAMA3,LAMB3,ITGA2,...
    ## 13:             MYC,TCF7,PIK3R1,RAF1,CCND1,CTNNA2,...
    ## 14:       IL6,CD86,FAS,HLA-DRB3,HLA-DRB1,HLA-DQA2,...
    ## 15:         IL1R1,IL6,KIT,HLA-DRB3,HLA-DRB1,IL1R2,...
    ## 16:      IL6,CD86,HLA-DRB3,HLA-DRB1,IL15RA,ICOSLG,...
    ## 17:            SPRY1,SOCS2,IL6,IL13RA2,SPRED1,MYC,...
    ## 18:          PDGFRA,IL1R1,FGF7,FGF2,DUSP5,RASGRP1,...
    ## 19:               IAPP,SLC2A2,PDX1,MAFA,NKX6-1,NKX2-2
    ## 20:            PDGFRA,FGF7,FGF2,HGF,CDKN1A,PIK3R1,...
    ## 21:      TNFRSF10D,FAS,TNFRSF10A,PIK3R1,RAF1,RAC1,...
    ## 22:                  PDGFRA,FGF7,FGF2,IL6,HGF,KIT,...
    ## 23:                ACSL4,LPL,ACSL5,ME1,PPARD,PCK1,...
    ## 24:           PDGFRA,TGFA,RELA,TCF7,CDKN1A,PIK3R1,...
    ## 25:                PDE1A,PDE7B,PDE3A,DCK,NT5E,PNP,...
    ## 26:                 UCK2,DCK,NT5E,PNP,CTPS1,CTPS2,...
    ## 27:           PDGFRA,FGF7,FGF2,BDKRB1,MRAS,PIK3R1,...
    ## 28:             PTGS2,MYC,RELA,LAMC2,PIK3R1,LAMA3,...
    ## 29:         RBM25,TCERG1,EFTUD2,DHX15,DDX42,RBM17,...
    ## 30:            C1S,CD86,C1R,HLA-DRB3,HLA-DRB1,SSB,...
    ## 31:              DCN,BMPR1B,THBS2,MYC,LTBP1,INHBA,...
    ## 32:            SOCS2,SLC2A2,PDX1,SOCS3,MAFA,ABCC8,...
    ## 33:     CD86,FAS,HLA-DRB3,HLA-DRB1,HLA-DQA2,HLA-E,...
    ## 34:         CD86,HLA-DRB3,HLA-DRB1,ABL1,CD55,RAC1,...
    ##                                           leadingEdge
    ##                                         topGenes
    ##  1:            ABCA6, ABCA8, ABCA9, ABCC8, ABCC6
    ##  2:                KIT, MYC, STAT5A, RELA, PPARD
    ##  3:          NLK, CSNK2A1, TGFBR2, YES1, NECTIN3
    ##  4:      CD86, FAS, HLA-DRB3, HLA-DRB1, HLA-DQA2
    ##  5:        HLA-DRB3, HLA-DRB1, CTSL, PSME3, LGMN
    ##  6:          IL1R1, IRAK2, IRAK3, TNFRSF10D, FAS
    ##  7: HLA-DRB3, HLA-DRB1, HLA-DQA2, CD40, HLA-DQA1
    ##  8:      CD86, FAS, HLA-DRB3, HLA-DRB1, HLA-DQA2
    ##  9:              MYC, TGFBR2, TCF7, PIK3R1, RAF1
    ## 10:                  A2M, C1S, TFPI, C1R, BDKRB1
    ## 11:              PDGFRA, IL1R1, IL18R1, IL6, HGF
    ## 12:            THBS2, THBS4, LAMC2, LAMA3, LAMB3
    ## 13:               MYC, TCF7, PIK3R1, RAF1, CCND1
    ## 14:           IL6, CD86, FAS, HLA-DRB3, HLA-DRB1
    ## 15:          IL1R1, IL6, KIT, HLA-DRB3, HLA-DRB1
    ## 16:        IL6, CD86, HLA-DRB3, HLA-DRB1, IL15RA
    ## 17:           SPRY1, SOCS2, IL6, IL13RA2, SPRED1
    ## 18:             PDGFRA, IL1R1, FGF7, FGF2, DUSP5
    ## 19:             IAPP, SLC2A2, PDX1, MAFA, NKX6-1
    ## 20:              PDGFRA, FGF7, FGF2, HGF, CDKN1A
    ## 21:      TNFRSF10D, FAS, TNFRSF10A, PIK3R1, RAF1
    ## 22:                 PDGFRA, FGF7, FGF2, IL6, HGF
    ## 23:                ACSL4, LPL, ACSL5, ME1, PPARD
    ## 24:             PDGFRA, TGFA, RELA, TCF7, CDKN1A
    ## 25:               PDE1A, PDE7B, PDE3A, DCK, NT5E
    ## 26:                  UCK2, DCK, NT5E, PNP, CTPS1
    ## 27:             PDGFRA, FGF7, FGF2, BDKRB1, MRAS
    ## 28:              PTGS2, MYC, RELA, LAMC2, PIK3R1
    ## 29:          RBM25, TCERG1, EFTUD2, DHX15, DDX42
    ## 30:           C1S, CD86, C1R, HLA-DRB3, HLA-DRB1
    ## 31:               DCN, BMPR1B, THBS2, MYC, LTBP1
    ## 32:             SOCS2, SLC2A2, PDX1, SOCS3, MAFA
    ## 33:      CD86, FAS, HLA-DRB3, HLA-DRB1, HLA-DQA2
    ## 34:         CD86, HLA-DRB3, HLA-DRB1, ABL1, CD55
    ##                                         topGenes

There are 34 significant pathways.

### KEGG Pathway Visualization

``` r
# arrange kegg_significant_pathways based on the p.adjust value. 
kegg_ordered <- kegg_significant_pathways %>%
  arrange(padj) 
```

``` r
# Plot the barplot. 
# Plot the 'size' column from the table onto x-axis, the the top 20 pathways on the y-axis
# Using the size because size represents the size of each pathway after removing genes not present in 'names(stats)'
kegg_plot<- ggplot(kegg_ordered, aes(x = size, y = reorder(pathway, padj), fill = padj)) +
  geom_bar(stat = "identity") +
  xlab("number of genes") +
  ylab("pathway") +  
  scale_fill_gradient(low = "gray70", high = "red") +
  ggtitle("KEGG GSEA result") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 4)) 
# scale_y_discrete(labels = lab)
kegg_plot
```

![](Gene-Set-Enrichment-Analysis_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
# #some examples from KEGG pathway
# plot1 <- plotEnrichment(kegg.lol[["KEGG_INSULIN_SIGNALING_PATHWAY"]],
#              deg.ranks) 
# plot2 <- plotEnrichment(kegg.lol[["KEGG_FATTY_ACID_METABOLISM"]],
#               deg.ranks)
```

``` r
geneids <-bitr(names(logFC_genes), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
```

    ## 'select()' returned 1:many mapping between keys and columns

    ## Warning in bitr(names(logFC_genes), fromType = "ENSEMBL", toType = "ENTREZID",
    ## : 14.53% of input gene IDs are fail to map...

``` r
#delete duplicate IDs
dedup_ids = geneids[!duplicated(geneids[c("ENSEMBL")]),]
# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
geneIDs <-
    ensembl_geneIDs[order(ensembl_geneIDs$logFC, decreasing = TRUE),
                   head(.SD, 1),
                   by = .(entrezgene_id)]

# gene_names = ensembl_geneIDs[ensembl_geneIDs$ensembl_gene_id %in% dedup_ids$ENSEMBL,]
# Create a vector of the gene unuiverse
kegg_gene_list <- geneIDs$logFC
# Name vector with ENTREZ ids
names(kegg_gene_list) <- geneIDs$entrezgene_id
# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)
# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

#plotting
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = "hsa",
               minGSSize    = 3,
               maxGSSize    = 80,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")
```

    ## Reading KEGG annotation online: "https://rest.kegg.jp/link/hsa/pathway"...

    ## Reading KEGG annotation online: "https://rest.kegg.jp/list/pathway/hsa"...

    ## Reading KEGG annotation online: "https://rest.kegg.jp/conv/ncbi-geneid/hsa"...

    ## preparing geneSet collections...

    ## GSEA analysis...

    ## Warning in preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam, : There are ties in the preranked stats (6.46% of the list).
    ## The order of those tied genes will be arbitrary, which may produce unexpected results.

    ## Warning in fgseaMultilevel(pathways = pathways, stats = stats, minSize =
    ## minSize, : For some pathways, in reality P-values are less than 1e-10. You can
    ## set the `eps` argument to zero for better estimation.

    ## leading edge analysis...

    ## done...

``` r
#enrichment map
emapplot(pairwise_termsim(kk2), cex.params = list(category_label = 0.5))
```

![](Gene-Set-Enrichment-Analysis_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
# Ridgeplot
ridgeplot(kk2) + labs(x = "enrichment distribution") + theme_ridges(font_size = 5)
```

    ## Picking joint bandwidth of 0.0866

![](Gene-Set-Enrichment-Analysis_files/figure-gfm/unnamed-chunk-20-2.png)<!-- -->

``` r
## isoleucine degradation is down-regulated - a branched amino acid

# https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
```

### Reactome `fgsea` Analysis

``` r
#Convert the REACTOME pathway tibble with the make.gs.lol function made previously. 
REACTOME.lol <- reactome.human.db %>% 
  dplyr :: select(gene_symbol, gs_name) %>% 
  make.gs.lol()
```

``` r
# GSEA with fgsea()
reactome.fgsea <- fgsea::fgsea(pathways = REACTOME.lol, stats = deg.ranks, scoreType = "pos")
```

    ## Warning in preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam, : There are ties in the preranked stats (84.06% of the list).
    ## The order of those tied genes will be arbitrary, which may produce unexpected results.

``` r
#Show top 5 genes for each pathway:
# Here, we are "pasting" the first five genes in the leadingEdge column
reactome.fgsea[,
           topGenes := paste0(head(unlist(`leadingEdge`), 5), collapse=", "),
           by = .(pathway)]
```

``` r
#Show result in table. 
reactome.fgsea %>%
    arrange(pval) %>%
    head(10) %>% 
    dplyr::select(-leadingEdge)
```

    ##                                                                                       pathway
    ##  1:                                                        REACTOME_SIGNALING_BY_INTERLEUKINS
    ##  2:                                              REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM
    ##  3:                                                   REACTOME_MAPK_FAMILY_SIGNALING_CASCADES
    ##  4:                                                                REACTOME_METABOLISM_OF_RNA
    ##  5:                                                                  REACTOME_RRNA_PROCESSING
    ##  6:                                                                       REACTOME_HEMOSTASIS
    ##  7:                                REACTOME_CONSTITUTIVE_SIGNALING_BY_ABERRANT_PI3K_IN_CANCER
    ##  8:                                      REACTOME_NEGATIVE_REGULATION_OF_THE_PI3K_AKT_NETWORK
    ##  9:                                                     REACTOME_PI3K_AKT_SIGNALING_IN_CANCER
    ## 10: REACTOME_DISEASES_OF_SIGNAL_TRANSDUCTION_BY_GROWTH_FACTOR_RECEPTORS_AND_SECOND_MESSENGERS
    ##             pval         padj   log2err        ES      NES size
    ##  1: 6.661038e-14 1.075758e-10 0.9545416 0.6763460 1.319481  434
    ##  2: 1.557360e-13 1.257569e-10 0.9436322 0.6425638 1.259444  668
    ##  3: 9.785811e-12 5.268028e-09 0.8753251 0.6881949 1.336816  308
    ##  4: 5.738441e-10 2.316895e-07 0.8012156 0.6266251 1.228465  656
    ##  5: 5.751561e-09 1.857754e-06 0.7614608 0.6945299 1.340612  205
    ##  6: 1.908050e-08 5.135835e-06 0.7337620 0.6174537 1.210239  624
    ##  7: 6.029913e-08 1.238458e-05 0.7049757 0.7859362 1.489510   73
    ##  8: 6.134778e-08 1.238458e-05 0.7049757 0.7473104 1.422251  104
    ##  9: 1.507012e-07 2.704249e-05 0.6901325 0.7526009 1.431750   94
    ## 10: 2.076085e-07 3.352877e-05 0.6901325 0.6347386 1.236032  397
    ##                               topGenes
    ##  1:    IL1R1, SOCS2, FGF2, IL18R1, IL6
    ##  2:    IL1R1, SOCS2, FGF2, IL18R1, IL6
    ##  3: PDGFRA, FGF7, FGF2, DUSP5, RASGRP1
    ##  4:      NIP7, IMP4, LTV1, BYSL, NOP14
    ##  5:      NIP7, IMP4, LTV1, BYSL, NOP14
    ##  6: A2M, PDE1A, RASGRP1, SERPINE2, HGF
    ##  7:       PDGFRA, FGF7, FGF2, HGF, KIT
    ##  8:       PDGFRA, FGF7, FGF2, HGF, KIT
    ##  9:       PDGFRA, FGF7, FGF2, HGF, KIT
    ## 10:       PDGFRA, FGF7, FGF2, HGF, KIT

``` r
# filter significant pathways. Set p<0.05. 
reactome_significant_pathways <- reactome.fgsea[padj < 0.05]
reactome_significant_pathways
```

    ##                                                                                       pathway
    ##  1:                                                           REACTOME_ADAPTIVE_IMMUNE_SYSTEM
    ##  2: REACTOME_ANTIGEN_ACTIVATES_B_CELL_RECEPTOR_BCR_LEADING_TO_GENERATION_OF_SECOND_MESSENGERS
    ##  3:               REACTOME_ANTI_INFLAMMATORY_RESPONSE_FAVOURING_LEISHMANIA_PARASITE_INFECTION
    ##  4:                             REACTOME_BINDING_AND_UPTAKE_OF_LIGANDS_BY_SCAVENGER_RECEPTORS
    ##  5:                                                     REACTOME_CD22_MEDIATED_BCR_REGULATION
    ##  6:                                   REACTOME_CELL_SURFACE_INTERACTIONS_AT_THE_VASCULAR_WALL
    ##  7:                                                 REACTOME_CHONDROITIN_SULFATE_BIOSYNTHESIS
    ##  8:                                  REACTOME_CHONDROITIN_SULFATE_DERMATAN_SULFATE_METABOLISM
    ##  9:                                                               REACTOME_COMPLEMENT_CASCADE
    ## 10:                                REACTOME_CONSTITUTIVE_SIGNALING_BY_ABERRANT_PI3K_IN_CANCER
    ## 11:                                                 REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS
    ## 12:                                              REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM
    ## 13:                                      REACTOME_DEFECTIVE_B4GALT7_CAUSES_EDS_PROGEROID_TYPE
    ## 14:                             REACTOME_DEFECTIVE_CHST14_CAUSES_EDS_MUSCULOCONTRACTURAL_TYPE
    ## 15:                                                      REACTOME_DEFECTIVE_CHSY1_CAUSES_TPBS
    ## 16:                                                    REACTOME_DERMATAN_SULFATE_BIOSYNTHESIS
    ## 17:                            REACTOME_DISEASES_ASSOCIATED_WITH_GLYCOSAMINOGLYCAN_METABOLISM
    ## 18:                                                        REACTOME_DISEASES_OF_GLYCOSYLATION
    ## 19: REACTOME_DISEASES_OF_SIGNAL_TRANSDUCTION_BY_GROWTH_FACTOR_RECEPTORS_AND_SECOND_MESSENGERS
    ## 20:                                                              REACTOME_EGFR_DOWNREGULATION
    ## 21:                                                          REACTOME_ELASTIC_FIBRE_FORMATION
    ## 22:                                                REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION
    ## 23:                                                 REACTOME_EXTRA_NUCLEAR_ESTROGEN_SIGNALING
    ## 24:                                                            REACTOME_FATTY_ACID_METABOLISM
    ## 25:                                                 REACTOME_FCERI_MEDIATED_CA_2_MOBILIZATION
    ## 26:                                                   REACTOME_FCERI_MEDIATED_MAPK_ACTIVATION
    ## 27:                                                  REACTOME_FCERI_MEDIATED_NF_KB_ACTIVATION
    ## 28:                                     REACTOME_FCGAMMA_RECEPTOR_FCGR_DEPENDENT_PHAGOCYTOSIS
    ## 29:                                                   REACTOME_FCGR3A_MEDIATED_IL10_SYNTHESIS
    ## 30:                                                                  REACTOME_FCGR_ACTIVATION
    ## 31:                                              REACTOME_FC_EPSILON_RECEPTOR_FCERI_SIGNALING
    ## 32:                                                     REACTOME_GLYCOSAMINOGLYCAN_METABOLISM
    ## 33:                                                                       REACTOME_HEMOSTASIS
    ## 34:         REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL
    ## 35:                                                               REACTOME_INFECTIOUS_DISEASE
    ## 36:                                                 REACTOME_INITIAL_TRIGGERING_OF_COMPLEMENT
    ## 37:                                                         REACTOME_INTERLEUKIN_10_SIGNALING
    ## 38:                                                   REACTOME_INTERLEUKIN_1_FAMILY_SIGNALING
    ## 39:                                                          REACTOME_INTERLEUKIN_1_SIGNALING
    ## 40:                                                  REACTOME_INTERLEUKIN_20_FAMILY_SIGNALING
    ## 41:                                       REACTOME_INTERLEUKIN_4_AND_INTERLEUKIN_13_SIGNALING
    ## 42:                                                   REACTOME_INTERLEUKIN_6_FAMILY_SIGNALING
    ## 43:                                                          REACTOME_INTERLEUKIN_7_SIGNALING
    ## 44:                                     REACTOME_INTRACELLULAR_SIGNALING_BY_SECOND_MESSENGERS
    ## 45:                                                     REACTOME_KERATAN_SULFATE_BIOSYNTHESIS
    ## 46:                                                             REACTOME_LEISHMANIA_INFECTION
    ## 47:                                                   REACTOME_MAPK_FAMILY_SIGNALING_CASCADES
    ## 48:                                                      REACTOME_METABOLISM_OF_CARBOHYDRATES
    ## 49:                                                                REACTOME_METABOLISM_OF_RNA
    ## 50:                                                     REACTOME_MET_ACTIVATES_PTK2_SIGNALING
    ## 51:                                                          REACTOME_MET_RECEPTOR_ACTIVATION
    ## 52:                                         REACTOME_MOLECULES_ASSOCIATED_WITH_ELASTIC_FIBRES
    ## 53:                                                                    REACTOME_MRNA_SPLICING
    ## 54:                                              REACTOME_NEGATIVE_REGULATION_OF_MAPK_PATHWAY
    ## 55:                                      REACTOME_NEGATIVE_REGULATION_OF_THE_PI3K_AKT_NETWORK
    ## 56:                                                               REACTOME_PARASITE_INFECTION
    ## 57:                                           REACTOME_PHOSPHOLIPASE_C_MEDIATED_CASCADE_FGFR2
    ## 58:                                                     REACTOME_PI3K_AKT_SIGNALING_IN_CANCER
    ## 59:                                                              REACTOME_PI_3K_CASCADE_FGFR2
    ## 60:                                  REACTOME_PROCESSING_OF_CAPPED_INTRON_CONTAINING_PRE_MRNA
    ## 61:                                                  REACTOME_REGULATION_OF_INSULIN_SECRETION
    ## 62:                                                                 REACTOME_RHO_GTPASE_CYCLE
    ## 63:                                    REACTOME_ROLE_OF_LAT2_NTAL_LAB_ON_CALCIUM_MOBILIZATION
    ## 64:                                            REACTOME_ROLE_OF_PHOSPHOLIPIDS_IN_PHAGOCYTOSIS
    ## 65:                                     REACTOME_RRNA_MODIFICATION_IN_THE_NUCLEUS_AND_CYTOSOL
    ## 66:                                                                  REACTOME_RRNA_PROCESSING
    ## 67:                                                   REACTOME_SCAVENGING_OF_HEME_FROM_PLASMA
    ## 68:                                                                REACTOME_SIGNALING_BY_EGFR
    ## 69:                                                     REACTOME_SIGNALING_BY_ERBB2_IN_CANCER
    ## 70:                                                                REACTOME_SIGNALING_BY_FGFR
    ## 71:                                                               REACTOME_SIGNALING_BY_FGFR2
    ## 72:                                                     REACTOME_SIGNALING_BY_FGFR_IN_DISEASE
    ## 73:                                                        REACTOME_SIGNALING_BY_INTERLEUKINS
    ## 74:                                                                 REACTOME_SIGNALING_BY_MET
    ## 75:                                           REACTOME_SIGNALING_BY_RECEPTOR_TYROSINE_KINASES
    ## 76:                                REACTOME_SIGNALING_BY_RHO_GTPASES_MIRO_GTPASES_AND_RHOBTB3
    ## 77:                                             REACTOME_SIGNALING_BY_THE_B_CELL_RECEPTOR_BCR
    ## 78:                               REACTOME_SYNTHESIS_OF_PROSTAGLANDINS_PG_AND_THROMBOXANES_TX
    ## 79:                      REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_DEATH_RECEPTORS_AND_LIGANDS
    ## 80:                                                     REACTOME_TRANSPORT_OF_SMALL_MOLECULES
    ##                                                                                       pathway
    ##             pval         padj   log2err        ES      NES size
    ##  1: 2.386391e-03 4.817527e-02 0.4317077 0.5582054 1.095261  763
    ##  2: 4.716860e-04 1.692829e-02 0.4984931 0.7185744 1.352728   59
    ##  3: 6.486408e-05 4.554587e-03 0.5384341 0.6486324 1.249932  186
    ##  4: 1.779807e-04 8.710267e-03 0.5188481 0.7088507 1.340675   68
    ##  5: 1.353004e-03 3.524356e-02 0.4550599 0.7659925 1.414541   34
    ##  6: 8.291281e-05 5.356168e-03 0.5384341 0.6594018 1.266777  157
    ##  7: 3.622760e-04 1.500194e-02 0.4984931 0.8400743 1.522160   20
    ##  8: 4.078448e-04 1.568260e-02 0.4984931 0.7387443 1.380181   50
    ##  9: 1.303688e-04 7.260193e-03 0.5188481 0.7055695 1.341787   85
    ## 10: 6.029913e-08 1.238458e-05 0.7049757 0.7859362 1.489510   73
    ## 11: 6.035192e-05 4.554587e-03 0.5573322 0.7796993 1.447679   43
    ## 12: 1.557360e-13 1.257569e-10 0.9436322 0.6425638 1.259444  668
    ## 13: 3.075937e-04 1.307273e-02 0.4984931 0.8445976 1.530356   20
    ## 14: 3.850604e-04 1.554681e-02 0.4984931 0.9402166 1.638798    8
    ## 15: 1.598428e-04 8.327294e-03 0.5188481 0.9488124 1.653781    8
    ## 16: 5.815903e-04 1.947061e-02 0.4772708 0.8930472 1.566974   11
    ## 17: 4.574272e-05 3.693725e-03 0.5573322 0.7912306 1.463694   40
    ## 18: 2.016416e-03 4.539631e-02 0.4317077 0.6367907 1.218728  137
    ## 19: 2.076085e-07 3.352877e-05 0.6901325 0.6347386 1.236032  397
    ## 20: 2.802525e-04 1.257244e-02 0.4984931 0.7982695 1.466730   29
    ## 21: 1.297780e-03 3.435925e-02 0.4550599 0.7327506 1.360508   43
    ## 22: 1.258343e-04 7.257945e-03 0.5188481 0.6217838 1.207523  288
    ## 23: 1.040070e-03 2.999486e-02 0.4550599 0.7034306 1.329399   67
    ## 24: 1.187332e-03 3.300459e-02 0.4550599 0.6287013 1.210025  173
    ## 25: 2.127409e-03 4.539631e-02 0.4317077 0.6979562 1.314161   60
    ## 26: 1.868427e-03 4.437513e-02 0.4550599 0.7006117 1.318913   59
    ## 27: 1.711790e-04 8.639189e-03 0.5188481 0.6759061 1.287152  107
    ## 28: 1.258343e-04 7.257945e-03 0.5188481 0.6890592 1.312823  109
    ## 29: 5.907490e-04 1.947061e-02 0.4772708 0.7124914 1.346709   65
    ## 30: 4.397429e-04 1.630782e-02 0.4984931 0.7662713 1.419808   41
    ## 31: 7.739238e-04 2.386185e-02 0.4772708 0.6369973 1.222979  154
    ## 32: 1.410233e-05 1.339721e-03 0.5933255 0.6948303 1.325905  121
    ## 33: 1.908050e-08 5.135835e-06 0.7337620 0.6174537 1.210239  624
    ## 34: 1.942421e-03 4.539631e-02 0.4550599 0.6350943 1.219276  155
    ## 35: 1.721163e-03 4.256677e-02 0.4550599 0.5594536 1.098601  861
    ## 36: 2.278600e-04 1.051411e-02 0.5188481 0.7447766 1.393006   51
    ## 37: 1.439722e-04 7.750501e-03 0.5188481 0.7615124 1.417612   46
    ## 38: 6.532622e-06 7.033456e-04 0.6105269 0.6962525 1.331016  132
    ## 39: 7.830825e-04 2.386185e-02 0.4772708 0.6742806 1.282903   95
    ## 40: 1.205740e-03 3.300459e-02 0.4550599 0.7911085 1.452155   25
    ## 41: 2.390413e-06 3.509561e-04 0.6272567 0.7266143 1.381933  101
    ## 42: 1.684347e-03 4.250344e-02 0.4550599 0.8048430 1.465007   22
    ## 43: 7.952867e-05 5.351617e-03 0.5384341 0.8117140 1.498974   34
    ## 44: 5.011816e-06 6.745069e-04 0.6105269 0.6351232 1.233298  283
    ## 45: 1.850019e-03 4.437513e-02 0.4550599 0.7817930 1.434700   27
    ## 46: 6.373604e-05 4.554587e-03 0.5384341 0.6267645 1.216358  265
    ## 47: 9.785811e-12 5.268028e-09 0.8753251 0.6881949 1.336816  308
    ## 48: 4.032879e-04 1.568260e-02 0.4984931 0.6098169 1.184153  286
    ## 49: 5.738441e-10 2.316895e-07 0.8012156 0.6266251 1.228465  656
    ## 50: 2.164406e-03 4.539631e-02 0.4317077 0.7794342 1.430726   25
    ## 51: 1.132109e-03 3.207642e-02 0.4550599 0.9578491 1.643615    5
    ## 52: 1.739571e-03 4.256677e-02 0.4550599 0.7509337 1.387768   36
    ## 53: 3.227267e-05 2.743177e-03 0.5573322 0.6556310 1.263333  185
    ## 54: 2.275399e-03 4.651606e-02 0.4317077 0.7225986 1.341659   43
    ## 55: 6.134778e-08 1.238458e-05 0.7049757 0.7473104 1.422251  104
    ## 56: 1.983859e-04 9.423328e-03 0.5188481 0.6968907 1.325413   84
    ## 57: 2.164406e-03 4.539631e-02 0.4317077 0.8269980 1.485829   17
    ## 58: 1.507012e-07 2.704249e-05 0.6901325 0.7526009 1.431750   94
    ## 59: 2.275399e-03 4.651606e-02 0.4317077 0.7987729 1.453958   22
    ## 60: 3.058893e-05 2.743177e-03 0.5573322 0.6399876 1.239341  235
    ## 61: 2.016416e-03 4.539631e-02 0.4317077 0.6756625 1.280519   73
    ## 62: 8.380351e-04 2.506346e-02 0.4772708 0.5823150 1.136123  436
    ## 63: 6.457014e-04 2.044721e-02 0.4772708 0.7417613 1.379172   45
    ## 64: 1.340764e-05 1.339721e-03 0.5933255 0.7697101 1.441887   53
    ## 65: 6.325239e-06 7.033456e-04 0.6105269 0.7711041 1.451889   60
    ## 66: 5.751561e-09 1.857754e-06 0.7614608 0.6945299 1.340612  205
    ## 67: 5.357968e-04 1.841089e-02 0.4772708 0.7626444 1.413087   41
    ## 68: 3.075937e-04 1.307273e-02 0.4984931 0.7478593 1.393024   47
    ## 69: 1.629123e-03 4.176244e-02 0.4550599 0.8013482 1.465955   24
    ## 70: 6.182252e-04 1.996867e-02 0.4772708 0.6915517 1.314099   82
    ## 71: 4.442998e-04 1.630782e-02 0.4984931 0.7067564 1.337100   69
    ## 72: 2.053414e-03 4.539631e-02 0.4317077 0.6990093 1.315896   59
    ## 73: 6.661038e-14 1.075758e-10 0.9545416 0.6763460 1.319481  434
    ## 74: 2.127409e-03 4.539631e-02 0.4317077 0.6746367 1.278575   73
    ## 75: 5.703091e-06 7.033456e-04 0.6105269 0.6109795 1.193453  487
    ## 76: 2.164406e-03 4.539631e-02 0.4317077 0.5627707 1.103326  691
    ## 77: 1.260964e-03 3.394095e-02 0.4550599 0.6483553 1.240827  136
    ## 78: 1.032177e-04 6.411408e-03 0.5384341 0.8956550 1.598061   15
    ## 79: 8.563526e-04 2.514563e-02 0.4772708 0.9099205 1.588603   10
    ## 80: 5.266381e-04 1.841089e-02 0.4772708 0.5705610 1.118964  705
    ##             pval         padj   log2err        ES      NES size
    ##                                                  leadingEdge
    ##  1:              RASGRP1,CD86,SOCS3,RNF4,SMURF1,HLA-DRB3,...
    ##  2:   IGHV3-13,IGHV3-48,PIK3R1,IGHV3-11,IGHV3-53,IGHV3-7,...
    ##  3:                      IAPP,IL6,PTGER2,GLP1R,GGT5,YES1,...
    ##  4: IGHV3-13,IGHV3-48,IGHV3-11,IGHV3-53,IGHV3-7,IGHV3-23,...
    ##  5: IGHV3-13,IGHV3-48,IGHV3-11,IGHV3-53,IGHV3-7,IGHV3-23,...
    ##  6:         SLC16A1,TNFRSF10D,TNFRSF10A,SLC7A6,YES1,GYPC,...
    ##  7:               DCN,CSGALNACT2,CHSY1,VCAN,CHST15,CHST9,...
    ##  8:                  DCN,GPC3,CSGALNACT2,CHSY1,VCAN,DSEL,...
    ##  9:          C1S,C1R,IGHV3-13,IGHV3-48,IGHV3-11,IGHV3-53,...
    ## 10:                        PDGFRA,FGF7,FGF2,HGF,KIT,CD86,...
    ## 11:          C1S,C1R,IGHV3-13,IGHV3-48,IGHV3-11,IGHV3-53,...
    ## 12:                      IL1R1,SOCS2,FGF2,IL18R1,IL6,HGF,...
    ## 13:                                      DCN,GPC3,VCAN,CSPG4
    ## 14:                                           DCN,VCAN,CSPG4
    ## 15:                                     DCN,CHSY1,VCAN,CSPG4
    ## 16:                                      DCN,VCAN,DSEL,CSPG4
    ## 17:                        OGN,PRELP,DCN,GPC3,CHSY1,VCAN,...
    ## 18:                     SPON1,OGN,PRELP,DCN,GPC3,ADAMTS9,...
    ## 19:                        PDGFRA,FGF7,FGF2,HGF,KIT,CD86,...
    ## 20:                    SPRY1,EREG,HBEGF,TGFA,PTPN3,SPRY2,...
    ## 21:                  FBLN1,MFAP4,FBLN5,LTBP1,TGFB3,ITGB5,...
    ## 22:                       FBLN1,MFAP4,A2M,FBLN5,FGF2,DCN,...
    ## 23:                      EREG,HBEGF,TGFA,ESR1,SRF,PIK3R1,...
    ## 24:                    ACSL4,PTGS2,PTGES,HADH,GGT5,ACSL5,...
    ## 25: IGHV3-13,IGHV3-48,IGHV3-11,IGHV3-53,IGHV3-7,IGHV3-23,...
    ## 26: IGHV3-13,IGHV3-48,IGHV3-11,IGHV3-53,IGHV3-7,IGHV3-23,...
    ## 27:     RASGRP1,RELA,IGHV3-13,IGHV3-48,IGHV3-11,IGHV3-53,...
    ## 28:          ABL1,YES1,IGHV3-13,IGHV3-48,PIK3R1,IGHV3-11,...
    ## 29:     YES1,IGHV3-13,IGHV3-48,IGHV3-11,IGHV3-53,IGHV3-7,...
    ## 30:     YES1,IGHV3-13,IGHV3-48,IGHV3-11,IGHV3-53,IGHV3-7,...
    ## 31:       RASGRP1,RELA,IGHV3-13,IGHV3-48,PIK3R1,IGHV3-11,...
    ## 32:                      OGN,PRELP,DCN,GPC3,ST3GAL1,HAS2,...
    ## 33:                  A2M,PDE1A,RASGRP1,SERPINE2,HGF,TFPI,...
    ## 34: IGHV3-13,IGHV3-48,IGHV3-11,IGHV3-53,IGHV3-7,IGHV3-23,...
    ## 35:                  IL1R1,IAPP,IL6,PTGER2,ST3GAL1,GLP1R,...
    ## 36:          C1S,C1R,IGHV3-13,IGHV3-48,IGHV3-11,IGHV3-53,...
    ## 37:                       IL1R1,IL6,CD86,PTGS2,LIF,IL1R2,...
    ## 38:                 IL1R1,IL18R1,IL1RL1,IRAK2,IRAK3,IL33,...
    ## 39:                   IL1R1,IRAK2,IRAK3,IL1R2,RELA,PSMC4,...
    ## 40:                SOCS3,STAT5A,IL24,JAK1,IL22RA1,IFNLR1,...
    ## 41:                       FGF2,IL6,HGF,IL13RA2,PTGS2,MYC,...
    ## 42:                        IL6,SOCS3,LIF,JAK1,CLCF1,OSMR,...
    ## 43:                    SOCS2,HGF,IL7R,STAT5A,PIK3R1,TSLP,...
    ## 44:                       PDGFRA,FGF7,FGF2,PDE1A,HGF,KIT,...
    ## 45:                              OGN,PRELP,ST3GAL1,CHST2,LUM
    ## 46:                      IAPP,IL6,PTGER2,GLP1R,GGT5,NT5E,...
    ## 47:                   PDGFRA,FGF7,FGF2,DUSP5,RASGRP1,IL6,...
    ## 48:                      OGN,PRELP,DCN,GPC3,ST3GAL1,HAS2,...
    ## 49:                      NIP7,IMP4,LTV1,BYSL,NOP14,DHX37,...
    ## 50:                    HGF,LAMC2,LAMA3,LAMB3,ITGA2,LAMA2,...
    ## 51:                                                  HGF,HPN
    ## 52:                  FBLN1,MFAP4,FBLN5,LTBP1,TGFB3,ITGB5,...
    ## 53:               HNRNPH1,EFTUD2,DHX15,CCAR1,PTBP1,DDX42,...
    ## 54:                   DUSP5,PTPN3,RAF1,DUSP6,PEBP1,DUSP1,...
    ## 55:                      PDGFRA,FGF7,FGF2,HGF,KIT,IL1RL1,...
    ## 56:        ABL1,YES1,IGHV3-13,IGHV3-48,IGHV3-11,IGHV3-53,...
    ## 57:                                     FGF7,FGF2,FGF10,FGF9
    ## 58:                        PDGFRA,FGF7,FGF2,HGF,KIT,CD86,...
    ## 59:                   FGF7,FGF2,PIK3R1,FGF10,GRB2,PTPN11,...
    ## 60:               HNRNPH1,EFTUD2,DHX15,CCAR1,PTBP1,DDX42,...
    ## 61:                 ACSL4,SLC2A2,GLP1R,ABCC8,FFAR1,STX1A,...
    ## 62:                 LETM1,SWAP70,GJA1,BASP1,RHOF,ARHGAP6,...
    ## 63:   IGHV3-13,IGHV3-48,PIK3R1,IGHV3-11,IGHV3-53,IGHV3-7,...
    ## 64:   IGHV3-13,IGHV3-48,PIK3R1,IGHV3-11,IGHV3-53,IGHV3-7,...
    ## 65:                      IMP4,NOP14,DHX37,GAR1,RCL1,DKC1,...
    ## 66:                      NIP7,IMP4,LTV1,BYSL,NOP14,DHX37,...
    ## 67: IGHV3-13,IGHV3-48,IGHV3-11,IGHV3-53,IGHV3-7,IGHV3-23,...
    ## 68:                   SPRY1,EREG,HBEGF,TGFA,PTPN3,PIK3R1,...
    ## 69:                 EREG,HBEGF,PIK3R1,PTPN12,CDC37,ERBIN,...
    ## 70:                FGF7,FGF2,SPRED1,HNRNPH1,PTBP1,PIK3R1,...
    ## 71:                 FGF7,FGF2,HNRNPH1,PTBP1,PIK3R1,SPRY2,...
    ## 72:               FGF7,FGF2,STAT5A,PIK3R1,LRRFIP1,POLR2A,...
    ## 73:                      IL1R1,SOCS2,FGF2,IL18R1,IL6,HGF,...
    ## 74:                    HGF,ARF6,LAMC2,PIK3R1,LAMA3,LAMB3,...
    ## 75:                       PDGFRA,SPRY1,FGF7,FGF2,HGF,KIT,...
    ## 76:                 LETM1,SWAP70,GJA1,BASP1,RHOF,ARHGAP6,...
    ## 77:       RASGRP1,RELA,IGHV3-13,IGHV3-48,PIK3R1,IGHV3-11,...
    ## 78:                            PTGS2,PTGES,PTGDS,PTGIS,HPGDS
    ## 79:                TNFRSF10D,FAS,TNFRSF10A,TP53BP2,TNFRSF10B
    ## 80:               A2M,SLC6A14,SLC2A2,ABCA6,ABCA8,SLC15A4,...
    ##                                                  leadingEdge
    ##                                            topGenes
    ##  1:              RASGRP1, CD86, SOCS3, RNF4, SMURF1
    ##  2:  IGHV3-13, IGHV3-48, PIK3R1, IGHV3-11, IGHV3-53
    ##  3:                  IAPP, IL6, PTGER2, GLP1R, GGT5
    ##  4: IGHV3-13, IGHV3-48, IGHV3-11, IGHV3-53, IGHV3-7
    ##  5: IGHV3-13, IGHV3-48, IGHV3-11, IGHV3-53, IGHV3-7
    ##  6:     SLC16A1, TNFRSF10D, TNFRSF10A, SLC7A6, YES1
    ##  7:            DCN, CSGALNACT2, CHSY1, VCAN, CHST15
    ##  8:              DCN, GPC3, CSGALNACT2, CHSY1, VCAN
    ##  9:          C1S, C1R, IGHV3-13, IGHV3-48, IGHV3-11
    ## 10:                    PDGFRA, FGF7, FGF2, HGF, KIT
    ## 11:          C1S, C1R, IGHV3-13, IGHV3-48, IGHV3-11
    ## 12:                 IL1R1, SOCS2, FGF2, IL18R1, IL6
    ## 13:                          DCN, GPC3, VCAN, CSPG4
    ## 14:                                DCN, VCAN, CSPG4
    ## 15:                         DCN, CHSY1, VCAN, CSPG4
    ## 16:                          DCN, VCAN, DSEL, CSPG4
    ## 17:                    OGN, PRELP, DCN, GPC3, CHSY1
    ## 18:                    SPON1, OGN, PRELP, DCN, GPC3
    ## 19:                    PDGFRA, FGF7, FGF2, HGF, KIT
    ## 20:                 SPRY1, EREG, HBEGF, TGFA, PTPN3
    ## 21:               FBLN1, MFAP4, FBLN5, LTBP1, TGFB3
    ## 22:                  FBLN1, MFAP4, A2M, FBLN5, FGF2
    ## 23:                    EREG, HBEGF, TGFA, ESR1, SRF
    ## 24:                 ACSL4, PTGS2, PTGES, HADH, GGT5
    ## 25: IGHV3-13, IGHV3-48, IGHV3-11, IGHV3-53, IGHV3-7
    ## 26: IGHV3-13, IGHV3-48, IGHV3-11, IGHV3-53, IGHV3-7
    ## 27:     RASGRP1, RELA, IGHV3-13, IGHV3-48, IGHV3-11
    ## 28:          ABL1, YES1, IGHV3-13, IGHV3-48, PIK3R1
    ## 29:    YES1, IGHV3-13, IGHV3-48, IGHV3-11, IGHV3-53
    ## 30:    YES1, IGHV3-13, IGHV3-48, IGHV3-11, IGHV3-53
    ## 31:       RASGRP1, RELA, IGHV3-13, IGHV3-48, PIK3R1
    ## 32:                  OGN, PRELP, DCN, GPC3, ST3GAL1
    ## 33:              A2M, PDE1A, RASGRP1, SERPINE2, HGF
    ## 34: IGHV3-13, IGHV3-48, IGHV3-11, IGHV3-53, IGHV3-7
    ## 35:               IL1R1, IAPP, IL6, PTGER2, ST3GAL1
    ## 36:          C1S, C1R, IGHV3-13, IGHV3-48, IGHV3-11
    ## 37:                    IL1R1, IL6, CD86, PTGS2, LIF
    ## 38:             IL1R1, IL18R1, IL1RL1, IRAK2, IRAK3
    ## 39:                IL1R1, IRAK2, IRAK3, IL1R2, RELA
    ## 40:              SOCS3, STAT5A, IL24, JAK1, IL22RA1
    ## 41:                  FGF2, IL6, HGF, IL13RA2, PTGS2
    ## 42:                    IL6, SOCS3, LIF, JAK1, CLCF1
    ## 43:                SOCS2, HGF, IL7R, STAT5A, PIK3R1
    ## 44:                  PDGFRA, FGF7, FGF2, PDE1A, HGF
    ## 45:                 OGN, PRELP, ST3GAL1, CHST2, LUM
    ## 46:                  IAPP, IL6, PTGER2, GLP1R, GGT5
    ## 47:              PDGFRA, FGF7, FGF2, DUSP5, RASGRP1
    ## 48:                  OGN, PRELP, DCN, GPC3, ST3GAL1
    ## 49:                   NIP7, IMP4, LTV1, BYSL, NOP14
    ## 50:                 HGF, LAMC2, LAMA3, LAMB3, ITGA2
    ## 51:                                        HGF, HPN
    ## 52:               FBLN1, MFAP4, FBLN5, LTBP1, TGFB3
    ## 53:            HNRNPH1, EFTUD2, DHX15, CCAR1, PTBP1
    ## 54:                DUSP5, PTPN3, RAF1, DUSP6, PEBP1
    ## 55:                    PDGFRA, FGF7, FGF2, HGF, KIT
    ## 56:        ABL1, YES1, IGHV3-13, IGHV3-48, IGHV3-11
    ## 57:                         FGF7, FGF2, FGF10, FGF9
    ## 58:                    PDGFRA, FGF7, FGF2, HGF, KIT
    ## 59:                 FGF7, FGF2, PIK3R1, FGF10, GRB2
    ## 60:            HNRNPH1, EFTUD2, DHX15, CCAR1, PTBP1
    ## 61:              ACSL4, SLC2A2, GLP1R, ABCC8, FFAR1
    ## 62:                LETM1, SWAP70, GJA1, BASP1, RHOF
    ## 63:  IGHV3-13, IGHV3-48, PIK3R1, IGHV3-11, IGHV3-53
    ## 64:  IGHV3-13, IGHV3-48, PIK3R1, IGHV3-11, IGHV3-53
    ## 65:                  IMP4, NOP14, DHX37, GAR1, RCL1
    ## 66:                   NIP7, IMP4, LTV1, BYSL, NOP14
    ## 67: IGHV3-13, IGHV3-48, IGHV3-11, IGHV3-53, IGHV3-7
    ## 68:                 SPRY1, EREG, HBEGF, TGFA, PTPN3
    ## 69:              EREG, HBEGF, PIK3R1, PTPN12, CDC37
    ## 70:              FGF7, FGF2, SPRED1, HNRNPH1, PTBP1
    ## 71:              FGF7, FGF2, HNRNPH1, PTBP1, PIK3R1
    ## 72:             FGF7, FGF2, STAT5A, PIK3R1, LRRFIP1
    ## 73:                 IL1R1, SOCS2, FGF2, IL18R1, IL6
    ## 74:                 HGF, ARF6, LAMC2, PIK3R1, LAMA3
    ## 75:                  PDGFRA, SPRY1, FGF7, FGF2, HGF
    ## 76:                LETM1, SWAP70, GJA1, BASP1, RHOF
    ## 77:       RASGRP1, RELA, IGHV3-13, IGHV3-48, PIK3R1
    ## 78:               PTGS2, PTGES, PTGDS, PTGIS, HPGDS
    ## 79:   TNFRSF10D, FAS, TNFRSF10A, TP53BP2, TNFRSF10B
    ## 80:              A2M, SLC6A14, SLC2A2, ABCA6, ABCA8
    ##                                            topGenes

There are 80 significant pathways from REACTOME ‘fgsea’ analysis.

### Reactome Visualization

``` r
# get the top 20 pathways from reactome_significant_pathways based on the p.adjusted value. 
top20_reactome <- reactome_significant_pathways %>%
  arrange(padj) %>%
  head(20)
```

``` r
lab1 = c('REACTOME_DISEASES_ASSOCIATED_WITH_ \n GLYCOSAMINOGLYCAN_METABOLISM', 'REACTOME_PROCESSING_OF_CAPPED_INTRON_ \n CONTAINING_PRE_MRNA', 'REACTOME_MRNA_SPLICING', 'REACTOME_ROLE_OF_PHOSPHOLIPIDS_IN_PHAGOCYTOSIS', 'REACTOME_GLYCOSAMINOGLYCAN_METABOLISM', 'REACTOME_SIGNALING_BY_RECEPTOR_TYROSINE_KINASES', 'REACTOME_RRNA_MODIFICATION_IN_THE_NUCLEUS_AND_CYTOSOL', 'REACTOME_INTERLEUKIN_1_FAMILY_SIGNALING', 'REACTOME_INTRACELLULAR_SIGNALING_BY_SECOND_ \n MESSENGERS', 'REACTOME_INTERLEUKIN_4_AND_INTERLEUKIN_13_SIGNALING', 'REACTOME_DISEASES_OF_SIGNAL_TRANSDUCTION_ \n BY_GROWTH_FACTOR_RECEPTORS_AND_SECOND_MESSENGERS', 'REACTOME_PI3K_AKT_SIGNALING_IN_CANCER', 'REACTOME_NEGATIVE_REGULATION_OF_THE_PI3K_ \n AKT_NETWORK', 'REACTOME_CONSTITUTIVE_SIGNALING_BY_ABERRANT_ \n P13K_IN_CANCER', 'REACTOME_HEMOSTASIS', 'REACTOME_RRNA_PROCESSING', 'REACTOME_METABOLISM_OF_RNA', 'REACTOME_MAPK_FAMILY_SIGNALING_CASCADES','REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM', 'REACTOME_SIGNALING_BY_INTERLEUKINS')

# Plot the barplot. 
# Plot the 'size' column from the table onto x-axis, the the top 20 pathways on the y-axis
# Using the size because size represents the size of each pathway after removing genes not present in 'names(stats)'
reactome_plot<- ggplot(top20_reactome, aes(x = size, y = reorder(pathway, padj), fill = padj)) +
  geom_bar(stat = "identity") +
  xlab("number of genes") +
  ylab("pathway") +  
  scale_fill_gradient(low = "gray70", high = "red") +
  ggtitle("REACTOME GSEA result") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 6)) + 
  scale_y_discrete(labels = lab1)
reactome_plot
```

![](Gene-Set-Enrichment-Analysis_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

``` r
reactome_gene_list <- geneIDs$logFC
# Name vector with ENTREZ ids
names(reactome_gene_list) <- geneIDs$entrezgene_id
# omit any NA values 
reactome_gene_list<-na.omit(reactome_gene_list)
# sort the list in decreasing order (required for clusterProfiler)
reactome_gene_list = sort(reactome_gene_list, decreasing = TRUE)
```

``` r
#plotting
reactome <- gsePathway(geneList = reactome_gene_list, 
                       organism = "human", 
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH", 
              verbose=FALSE)
```

    ## Warning in preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam, : There are ties in the preranked stats (6.46% of the list).
    ## The order of those tied genes will be arbitrary, which may produce unexpected results.

    ## Warning in fgseaMultilevel(pathways = pathways, stats = stats, minSize =
    ## minSize, : For some pathways, in reality P-values are less than 1e-10. You can
    ## set the `eps` argument to zero for better estimation.

### GO `fgsea` Analysis

``` r
#Convert the GO tibble with the make.gs.lol function made previously. 
GO.lol <- GO.human.db %>% 
  dplyr :: select(gene_symbol, gs_name) %>% 
  make.gs.lol()
```

``` r
# GSEA with fgsea()
GO.fgsea <- fgsea::fgsea(pathways = GO.lol, stats = deg.ranks, scoreType = "pos")
```

    ## Warning in preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam, : There are ties in the preranked stats (84.06% of the list).
    ## The order of those tied genes will be arbitrary, which may produce unexpected results.

``` r
#Show top 5 genes for each pathway:
# Here, we are "pasting" the first five genes in the leadingEdge column
GO.fgsea[,
           topGenes := paste0(head(unlist(`leadingEdge`), 5), collapse=", "),
           by = .(pathway)]
```

``` r
#Show result in table. 
GO.fgsea %>%
    arrange(pval) %>%
    head(10) %>% 
    dplyr::select(-leadingEdge) %>% 
    knitr::kable()
```

| pathway                                           | pval | padj |   log2err |        ES |      NES | size | topGenes                              |
|:--------------------------------------------------|-----:|-----:|----------:|----------:|---------:|-----:|:--------------------------------------|
| GOBP_REGULATION_OF_IMMUNE_SYSTEM_PROCESS          |    0 |    0 | 1.1866510 | 0.6331397 | 1.247325 | 1313 | IL1R1, SLIT2, MFAP4, A2M, NR4A3       |
| GOBP_NEGATIVE_REGULATION_OF_RESPONSE_TO_STIMULUS  |    0 |    0 | 1.0864405 | 0.6119568 | 1.207090 | 1567 | SERPINF1, SFRP4, PDGFRA, SPRY1, SLIT2 |
| GOBP_IMMUNE_EFFECTOR_PROCESS                      |    0 |    0 | 1.0672100 | 0.6692079 | 1.312346 |  594 | IL1R1, MFAP4, A2M, NR4A3, IL18R1      |
| GOBP_POSITIVE_REGULATION_OF_IMMUNE_SYSTEM_PROCESS |    0 |    0 | 1.0672100 | 0.6456084 | 1.269593 |  825 | IL1R1, MFAP4, A2M, NR4A3, IL18R1      |
| GOBP_REGULATION_OF_CELLULAR_COMPONENT_MOVEMENT    |    0 |    0 | 1.0276699 | 0.6274816 | 1.235432 | 1036 | SERPINF1, PDGFRA, IL1R1, SLIT2, FBLN1 |
| GOBP_CELL_MIGRATION                               |    0 |    0 | 1.0276699 | 0.6122145 | 1.206981 | 1457 | SERPINF1, PDGFRA, IL1R1, SLIT2, CH25H |
| GOBP_REGULATION_OF_CELL_ACTIVATION                |    0 |    0 | 1.0073180 | 0.6630151 | 1.300317 |  598 | PDGFRA, NR4A3, RASGRP1, IL6, SERPINE2 |
| GOBP_REGULATION_OF_CELL_POPULATION_PROLIFERATION  |    0 |    0 | 1.0073180 | 0.6020832 | 1.187810 | 1630 | SERPINF1, SFRP4, PDGFRA, SPRY1, OGN   |
| GOBP_LYMPHOCYTE_ACTIVATION                        |    0 |    0 | 1.0073180 | 0.6491489 | 1.274578 |  705 | IL18R1, RASGRP1, IL6, KIT, GPR183     |
| GOBP_LOCOMOTION                                   |    0 |    0 | 0.9969862 | 0.5980782 | 1.180241 | 1808 | SERPINF1, PDGFRA, IL1R1, SLIT2, CH25H |

``` r
# filter significant pathways. Set p<0.05. 
GO_significant_pathways <- GO.fgsea[padj < 0.05]
GO_significant_pathways
```

    ##                                            pathway         pval         padj
    ##    1:       GOBP_ACTIN_CYTOSKELETON_REORGANIZATION 3.570318e-03 4.397801e-02
    ##    2:            GOBP_ACTIN_FILAMENT_BASED_PROCESS 5.922388e-05 2.234762e-03
    ##    3:             GOBP_ACTIN_FILAMENT_ORGANIZATION 2.460386e-03 3.426160e-02
    ##    4:           GOBP_ACTIVATION_OF_IMMUNE_RESPONSE 4.725465e-07 4.404076e-05
    ##    5: GOBP_ACTIVATION_OF_PROTEIN_KINASE_C_ACTIVITY 1.518675e-03 2.414740e-02
    ##   ---                                                                       
    ## 1326:                                   HP_UVEITIS 1.316188e-03 2.175507e-02
    ## 1327:                 HP_VASCULAR_SKIN_ABNORMALITY 1.150517e-03 2.002210e-02
    ## 1328:                         HP_VISUAL_IMPAIRMENT 2.793365e-03 3.741658e-02
    ## 1329:                               HP_WEIGHT_LOSS 7.098125e-04 1.448748e-02
    ## 1330:                         HP_WIDE_NASAL_BRIDGE 3.394917e-04 8.403642e-03
    ##         log2err        ES      NES size                             leadingEdge
    ##    1: 0.4317077 0.6483794 1.236872  101    PDGFRA,FGF7,KIT,HMCN1,ANXA1,ABL1,...
    ##    2: 0.5573322 0.5780178 1.136151  760    PDGFRA,SLIT2,FGF7,KIT,SFRP1,ARF6,...
    ##    3: 0.4317077 0.5827269 1.137876  425  SLIT2,SFRP1,ARF6,SWAP70,HMCN1,RHOF,...
    ##    4: 0.6749629 0.6441243 1.252081  316     MFAP4,A2M,NR4A3,C1S,C1R,RABGEF1,...
    ##    5: 0.4550599 0.9634775 1.657519    4                               CD86,ABL1
    ##   ---                                                                          
    ## 1326: 0.4550599 0.7355416 1.375650   44 RASGRP1,IL6,FAS,HLA-DRB1,TCF4,MALT1,...
    ## 1327: 0.4550599 0.5789632 1.132771  503     RASGRP1,C1S,KIT,C1R,GJA1,ANTXR2,...
    ## 1328: 0.4317077 0.5499370 1.082959 1079  DCN,SEMA3A,DLL1,GJA1,CHRDL1,TWIST2,...
    ## 1329: 0.4772708 0.6067918 1.177331  275          IL6,KIT,GPC3,PDX1,GJA1,FAS,...
    ## 1330: 0.4984931 0.5819746 1.140220  544     KIT,GPC3,SOX6,LETM1,GJA1,TWIST2,...
    ##                                topGenes
    ##    1:   PDGFRA, FGF7, KIT, HMCN1, ANXA1
    ##    2:   PDGFRA, SLIT2, FGF7, KIT, SFRP1
    ##    3: SLIT2, SFRP1, ARF6, SWAP70, HMCN1
    ##    4:       MFAP4, A2M, NR4A3, C1S, C1R
    ##    5:                        CD86, ABL1
    ##   ---                                  
    ## 1326: RASGRP1, IL6, FAS, HLA-DRB1, TCF4
    ## 1327:      RASGRP1, C1S, KIT, C1R, GJA1
    ## 1328:   DCN, SEMA3A, DLL1, GJA1, CHRDL1
    ## 1329:        IL6, KIT, GPC3, PDX1, GJA1
    ## 1330:      KIT, GPC3, SOX6, LETM1, GJA1

There are 1262 significant pathways present in the GO gene sets.

### GO Analysis Visualization

``` r
# get the top 20 pathways from GO_significant_pathways based on the p.adjust value. 
top20_GO <- GO_significant_pathways %>%
  arrange(padj) %>%
  head(20)
```

``` r
# Names of the pathways are too long, so breaking them into two lines over here for grpah later. 
labs = c('GOBP_REGULATION_OF_CELL_DIFFERENTIATION', 'GOBP_NEGATIVE_REGULATION_OF_MULTICELLULAR_\n
ORGANISMAL_PROCESS',
         'GOBP_REGULATION_OF_RESPONSE_TO_STRESS',
         'GOBP_POSITIVE_REGULATION_OF_DEVELOPMENTAL_PROCESS',
         'GOBP_B_CELL_ACTIVATION',
         'GOBP_REGULATION_OF_RESPONSE_TO_EXTERNAL_STIMULUS',
         'GOBP_RIBOSOME_BIOGENESIS',
         'GOBP_IMMUNE_RESPONSE',
         'GOBP_CELL_ACTIVATION',
         'GOBP_IMMUNE_SYSTEM_DEVELOPMENT',
         'GOBP_REGULATION_OF_CELL_POPULATION_PROLIFERATION',
         'GOBP_REGULATION_OF_CELL_ACTIVATION',
         'GOBP_LYMPHOCYTE_ACTIVATION',
         'GOBP_LOCOMOTION',
         'GOBP_REGULATION_OF_CELLULAR_COMPONENT_MOVEMENT',
         'GOBP_POSITIVE_REGULATION_OF_IMMUNE_SYSTEM_PROCESS',
         'GOBP_IMMUNE_EFFECTOR_PROCESS',
         'GOBP_CELL_MIGRATION',
         'GOBP_NEGATIVE_REGULATION_OF_RESPONSE_TO_STIMULUS',
         'GOBP_REGULATION_OF_IMMUNE_SYSTEM_PROCESS')
      
# Plot the barplot. 
# Plot the 'size' column from the table onto x-axis, the the top 20 pathways on the y-axis
# Using the size because size represents the size of each pathway after removing genes not present in 'names(stats)'
go_barplot<-ggplot(top20_GO, aes(x = size, y = reorder(pathway, padj), fill = padj)) +
  geom_bar(stat = "identity") +
  xlab("number of genes") +
  ylab("pathway") +  
  scale_fill_gradient(low = "gray70", high = "red") +
  ggtitle("GO GSEA result") +
  theme(axis.text.y = element_text(size = 6)) +
  scale_y_discrete(labels = labs)
go_barplot
```

![](Gene-Set-Enrichment-Analysis_files/figure-gfm/unnamed-chunk-36-1.png)<!-- -->

### Disease Association with T2D based on the DE genes Using DisGeNET

``` r
FDR <- 0.05
filtered_genes <- ensembl_geneIDs %>% filter(P.Value < FDR)

filtered_genes <- filtered_genes[order(filtered_genes$P.Value),]
filtered_genes$entrezgene_id[1:30]
```

    ##  [1]  5176  5156  3554 10418 79827  6424 10252  4253  3375  8013  4969 54498
    ## [13]  5136  9353  8835 56605  2252  4239  2182 10516  2192     2  9023  2247
    ## [25] 11254 23551  8809 84935  1847 10514

``` r
# de <- names(geneList)[abs(geneList) > 2]
unique(filtered_genes$entrezgene_id)
```

    ##    [1]      5176      5156      3554     10418     79827      6424     10252
    ##    [8]      4253      3375      8013      4969     54498      5136      9353
    ##   [15]      8835     56605      2252      4239      2182     10516      2192
    ##   [22]         2      9023      2247     11254     23551      8809     84935
    ##   [29]      1847     10514      6926     79987     10125       716      3569
    ##   [36]     25890      5270      3082      9501      3815     24147     23460
    ##   [43]      4929      6514     10351      5549     83716      9173      1880
    ##   [50]      6422     27115     64094    401145      1634      5732      3656
    ##   [57]      8406      7035      2719     11213      5996     26002     85027
    ##   [64]     11228    121260    387763     56944    204851       942     54663
    ##   [71]      3598    284716      1734      3651     10371       715      4883
    ##   [78]      5743      2069     59271        NA    375757     56999      6566
    ##   [85]       623    161742     90865      4673     22822      3202       658
    ##   [92]     23212      9536      6482      7373      2740      7371     28514
    ##   [99]     24145      3489      6003      7058     23414     55553      3182
    ##  [106]      3033      9358     10350      2028     89795     84912      4609
    ##  [113]       382      5329      1588    345275     55022     51301      3954
    ##  [120]     23413     23076      9021     51388     23075      5272      4023
    ##  [127]     80273    254428     23258      3976      8073      8808     51491
    ##  [134]      2687     92856      2014      2697      5139     51701      1633
    ##  [141]     55661    389558     91851      4052     54491      6781      5372
    ##  [148]     27342      3037     23560    117581 100033814    116039     79698
    ##  [155]     84549     51703       558    389692     84946      2117      9770
    ##  [162]     79974       687     96626      3987 100288695     58517     57214
    ##  [169]      5225       705    338557      1024      3164       660      8602
    ##  [176]     55829      4907      6047      9047     83872     10409      8793
    ##  [183]     56652     55454     57647     79642      1601      2274      5730
    ##  [190]     26207     11320     23223     54733    118429     55107      9188
    ##  [197]    283927     79582      8682     57497      3624    144195     56919
    ##  [204]     54509      1839     10752      7039       355     23743     22856
    ##  [211]    203859 102723370      8797      4283    347731      9595      4747
    ##  [218]     90161      5552      4354      1910      9507      1051     84519
    ##  [225]     25925     57154      6833      3187      8676     27248       347
    ##  [232]     58475     84689      9435      1009    140738    122786      5638
    ##  [239]     57419     55612     10915      8324       395      4319      9057
    ##  [246]      3125      3123      6103      1105      6707      7850    165918
    ##  [253]    284695      2042 124900260     64092     55379     57393    285203
    ##  [260]      9343    353355      5774      3575     51114    143458      6776
    ##  [267]      1973       604      8404       654 100129842      1665      6541
    ##  [274]      8706    283106      1457     55749       489      9592      1462
    ##  [281]      7060       301     51650    203111     55196    133308     29960
    ##  [288]     79691        25     23082     23600      5858      5793      5502
    ##  [295]     57157    130074     96764    387893     57608 100288413      8942
    ##  [302]     23481    143879     54943     11199    150094 102724428     10156
    ##  [309] 100271927      7048     24137      8886     83440     11009    134147
    ##  [316]     79134      2099     92126      3428      1794     23268      2328
    ##  [323]      7525    406947    114614      9448       976     84870      9267
    ##  [330]      5725      5970     11112      2995      9095      2741      9671
    ##  [337]     55788      4199     65983 100130302    170067     55177      5467
    ##  [344]      3918       196      5420     29114     84883     29015    221409
    ##  [351]     65985     51365     57511      5570     54433      9945      5271
    ##  [358]    342184     10171     22808    144165     30817      6418    646817
    ##  [365]     55646     10228     54865     51260      1736     28958     23328
    ##  [372]     64005      4860      6715     25945      9120     25921     65062
    ##  [379]      6722      1805      8568      2286 114841035     10512     56964
    ##  [386]     27000     10691      6059      6741     51018    283349      2004
    ##  [393]     11177       275     89894     23017      6421    132864     25937
    ##  [400]      5122     54517     57602      7052      6932       639     27166
    ##  [407]     10630      1236     79649     84217      5261      1955      4864
    ##  [414]     51205    643834      5222    643847 124906959 124906960     91316
    ##  [421]    728411    387036      7030      1026     23514      6515      5295
    ##  [428]     23246      9697      1503     51569     85480      2553      3159
    ##  [435]    128178      9816      7538     57478    124491    154810      4825
    ##  [442]      9662     26574      1438     10253      9732      1514     11325
    ##  [449]     81575    221692     26253     55720      7159    654502     29970
    ##  [456] 100505385     54888    347733      4869     83931      3936      9208
    ##  [463]     83860     23764      6695      3604     80218     23559      4170
    ##  [470]    574406     57685     55127    151647     26156      4232     27097
    ##  [477]     23160      3707     56902      4338      1604     27440      5376
    ##  [484]      6367    402055      6092     27065     79634      7832     83607
    ##  [491]     10885    326340     54852      5740      3098     57669    391059
    ##  [498]      5142      3035     91039    164284     83857     10184        38
    ##  [505]     10438    403341 100033808     27067      9325      5396    286749
    ##  [512]     49854    223082      4478    148345    116150     10969      3716
    ##  [519]      2034      7805     10123     56474     10103      3762     57586
    ##  [526]      8877      1662     23760      3909      2521     10040      2012
    ##  [533]      5704     58527      2864       663      5172      5052     79178
    ##  [540]     10921     81887    222234      1965     29964      8495      5791
    ##  [547]       730    134285      5833      1382 100033437      1969      1902
    ##  [554]     23433     84879       956      4254      5806     27111     55216
    ##  [561]      5822 102724159      9889      2620      5208       699    162394
    ##  [568]      6667 100033445    266722      5105      5937      5757     22934
    ##  [575]     81127     79892      3185    153443    148398      3561     87178
    ##  [582]      9308      5305      6943      6804    148789      8613     84365
    ##  [589]     11214    143684     85450      8565     55137      9976     55003
    ##  [596]     28969     57453     55824     23034     25824     83795      4256
    ##  [603]    200728      4691      4215    219654      3512    401665      3059
    ##  [610]     23604     27306      8508     78987     55033      2869     51118
    ##  [617]     23645    387921      5894     23151     51167 123956252 100033801
    ##  [624]     23529     79034      3096      9459     84991     58985     51602
    ##  [631]    654429      7050     57050     22821    140901 100033815    143903
    ##  [638]    196527      7798    317772 100287896      6374      3773      8828
    ##  [645]    129049     54732      9180     65975     51593    448831    441581
    ##  [652]     10978     80176    150465     22943     23519     51050      4234
    ##  [659]     55233    124801      9221    140699     54543 107986346      6925
    ##  [666]    220108     93349     65083      3667      2005    221322     26577
    ##  [673]     79901     23219      3914     51429      1848      9001      9388
    ##  [680]     54973     10901    317749     25827      5037      8496     84276
    ##  [687]     84858       390     11021      9649      5902    284099      3673
    ##  [694]    113235     92002     55260     84617    126433     26155      7436
    ##  [701]     81929      5575 100033818    448834      4489     55650      9262
    ##  [708] 100506581     81790      5782      1936      6583    203238    145581
    ##  [715]      4899      5321    445815      3840     10658      3075     84970
    ##  [722]     57705      5690       968     10425      5879    169200     25841
    ##  [729]     51363      6698     59084    122060     84916     84866     25837
    ##  [736]      2048      2669    152519      2889     26094     79905      3908
    ##  [743]    117246      1984    143244    119391      4725     84154     78989
    ##  [750]     89857     22874     55243     11013    387700      9219     10451
    ##  [757]    161357     65084    345557      5217     10938     22843      7430
    ##  [764]    148811      1843     10801     10046     10186     56906     54716
    ##  [771]    115004     29914      7756     85477      8683     10615      6318
    ##  [778]    200539     51239     79745        97     23098       157      1802
    ##  [785]     26610      9669      8445     51208     79000       845     10529
    ##  [792]      9820      8844     53919      3589    339665      2788     10916
    ##  [799]      6875     29767      1909      5214     23132     28232    253827
    ##  [806]     79018     23515     54796      2118    642987      2124     51116
    ##  [813]      8531     84790 100033809     54906      3755      9508     83539
    ##  [820]     64921      3572     54617      3467      9277      7980      4973
    ##  [827]      9146      4707      3557      1996    115209      8661     51201
    ##  [834]     10197     83706      6924     51421      2153     26097     23099
    ##  [841]     79159     23590      4314    253012     51726      9187     63920
    ##  [848]     64425      7043      4060     54806     83732     10468      1464
    ##  [855]      3099     11140     10892       684      1429      7345      9149
    ##  [862]      8514      4588      6641      2531     26064     64859      2533
    ##  [869]    140576      7716     84293      8795      1556     79022      4593
    ##  [876]      9019      6863     51227     80117     84058      8572     10019
    ##  [883]      9844      6014      4345     22914 100528032      9903       308
    ##  [890]    125113      6846     65980     57179       846     51339     55573
    ##  [897]     64081     51232     79635      5884     54434    124808     29889
    ##  [904]     79711      1723     57818      5972      2888      3601      5210
    ##  [911]     64149     90843    130271      9275      5717      3480     51274
    ##  [918]     11237     57552      7458    163259      6723      4247 100033444
    ##  [925]       384    256586    117154      3693      5649     54815      2107
    ##  [932]     80183      9674      4839      4616      8710 124901197     57134
    ##  [939]     11004     51330     83481     54883     11343     10590     10492
    ##  [946]    653857     81035     22824      2214      2215 124905743     23479
    ##  [953]      1622      7975      6342      1666     80036     57544      6281
    ##  [960]     23101     56259     10189      4821      4512      4513      4509
    ##  [967]     10322      2819    284098     79183    282679      5496     25936
    ##  [974]    246243      1819      2554     54206      5641     54861     84319
    ##  [981]     57731       486 100533181      5137     26297      3074     79703
    ##  [988]     81627     55183     27077      1654     55131     23704      6332
    ##  [995]      2300      6670     51700    152926     79873      3487     10946
    ## [1002]      1469     10383      4313      9819     28987    118425      3939
    ## [1009]      2938      3687      9659 124904395      9122     57687     23049
    ## [1016]    114789        36    154075     64651      2202     65260      6627
    ## [1023]     80155      2713       374      8368    554313     25932      7272
    ## [1030]    285016     83690    256691      9689      2983     11051 102723996
    ## [1037]     23308    114897     10007     10849     93134      2710 102724908
    ## [1044]     55010     51435     23389      1268     10586     79660     79834
    ## [1051]     84668      7913    257240    339803        18    134701      7456
    ## [1058]     84954     78991     64395     51621    158160     51478     51267
    ## [1065]      7223    254263      5430      2053      3609    346171     23266
    ## [1072]     55344      3669    152189      8139     29098 100033805      4258
    ## [1079]     11333    386680     55283    388882     10993      1660      8747
    ## [1086]    145241     55341     51776    391253     79805       421    401474
    ## [1093]      9314      7307 102724594     11010      2865      2866      4678
    ## [1100]     29085    389813    440173      4093     22981      1962    693149
    ## [1107]    131616     54957     57125    113220     26958     90627     51278
    ## [1114]    126248      7306    648987     10166     57405     91647    140885
    ## [1121]      5991      3479    259266     57584      6480      8932      7803
    ## [1128]      7750      3419    222658     51084     83743     30836     11279
    ## [1135]     55844     78992      9551     92912     26130      6258    152330
    ## [1142]      7574     64109     54039 100033822      2057      5087      2653
    ## [1149]        27     10940     23194     51154    150590 100033441 100033443
    ## [1156] 100033453 100033603 100036565 100033804 100033816     55592      4643
    ## [1163]    163702      3117      3118     84134    146850    127018     83604
    ## [1170]       966     22807 100616475     11218    127435      1675     27087
    ## [1177]     90639     27340     56675       483      5128     51760     81848
    ## [1184]     23210     10482      5871      3566      1800      5501     55966
    ## [1191]     22889    677771    284029     80169    348980    114795      2355
    ## [1198]      2913     55198     58504    254251      1316    195814      4320
    ## [1205]     23276    112616      9842 100132354      3620     55500      8799
    ## [1212] 100131211      3133     27161 100996928     51631     92291     84081
    ## [1219] 100033440 100033442 100033446 100033448 100033449 100033455 100033456
    ## [1226] 100033458 100033460 100033799 100033802 100033803 100033810 100033813
    ## [1233] 100033817      5360     57480      2255     23528    116285     84450
    ## [1240]     10568       908      9524      2006      9705    338433 100033447
    ## [1247] 100033450 100033451 100033454 100033812      8541      5775      4245
    ## [1254]     58472      6700     51062     29068       587 100534599     57461
    ## [1261]      3313     84133      6626      4149     56262      6615     59272
    ## [1268]     23122    160760    150221    440804     85376      1193      5950
    ## [1275]     85478     54762    201176     81624       958     56129      3692
    ## [1282]      9880     10813     10203    123879    440093     85315       595
    ## [1289]    139189       328     10846     90632     27098      2590      7851
    ## [1296]     84316     10150      8388       550    138240     27241      9498
    ## [1303]     10311     57523      1075      4071     79772     10979      7305
    ## [1310]      9948    117177    285753     54660     56288     56967     93587
    ## [1317]     64841      9693      5266    129285     79641      9415     55008
    ## [1324]      8884      4790      3078      7291       995      8940     64682
    ## [1331]     55625      4750     10170      3576     64850     11337      6235
    ## [1338]     51379    401398 101928283     55658      4928     79677     10721
    ## [1345]     64783    121512      5784     54928      3977      5530      6455
    ## [1352]      8662    170949 100287083     57488    221830       187     65095
    ## [1359]      9899     51390      2642     57228     10640      1027      1958
    ## [1366]    113612     10893     57003       972     54941     10236    728239
    ## [1373]     81557      2791    221749      4885    767605    767606      5713
    ## [1380]      9684    113263    677849 100507246     26503     56265     56145
    ## [1387]      9752     56139     56138     57467      6840     84919      3184
    ## [1394]    493753 100036564    345778      9861      5411      3417     80017
    ## [1401]     56683 110091775     23603      9929     55856     65243      9404
    ## [1408]      3586     51056       701      9510       773     27292     54805
    ## [1415]     83858     10420    283489        34     55005    387742      4175
    ## [1422]     55228    113146      6317     50618     59348     60412      6356
    ## [1429]     57101      6319      3758     11156    144453      3559    406912
    ## [1436]    388815     54882    404734    285343     55666      6675     54873
    ## [1443]     80331       501     23214     55239    245915      1496      7001
    ## [1450]      7675      4921      9128      4670      5228      5522      6427
    ## [1457]      6257      8787     84890     23271     83464      4070     60686
    ## [1464]     55607    731157    254042     11100    406998    115701     84284
    ## [1471]    692312     56342 100861548    259291    259289    150084    127534
    ## [1478]     57536     55142      6876      3655    440888     57126       814
    ## [1485]     10785    199675     51321      4329     10324     51304      5606
    ## [1492]     81027     79707     54084    134637      6375      1543      8882
    ## [1499]     27091     10178      3038      5027    653604    440686    333932
    ## [1506]    126961    360023     83540     11102 109703458    144501     64499
    ## [1513]      7177    140465     57696    161823    145942     79954     60468
    ## [1520]      1997     57510     55526     84886       152     25804      8767
    ## [1527]     55819      8869     64147      3276     79039     27089     79155
    ## [1534]     79627     50651     64756      7567 124903710    643236      8761
    ## [1541]     26219    283129     29881      5036     56832      9935      7071
    ## [1548]      3208    205428      3659    124989      6789    414236    151295
    ## [1555]      9258      8833     57532      3937      7091      9462    127495
    ## [1562]     29081      7378     55054    134728     10435     65065      6560
    ## [1569]     81797      3726      9901     10320      8939     85025     10245
    ## [1576]      9987     64854     84153      3122     23597     85458      2954
    ## [1583]      4773      8968    142679    127077      4088      2939      9720
    ## [1590]     91369    221393     23473     84803      5934      2023    440712
    ## [1597]      8669    400073     29937      8555     63904      9322    148581
    ## [1604]      9792     54982     56931     25900      1454 102800317     80209
    ## [1611]      5655      8667      2166      1174    406897    768215    406898
    ## [1618]      1010    541472     26354      1797     55759    166647      3776
    ## [1625]    347688     54881     51015    374308      2114    729614    116984
    ## [1632]      6456     23350      6366      9572    728661      1719      6414
    ## [1639]      3139    114791    118813    157574     85463     29959     29888
    ## [1646]     11103     11130    285521      2151     81875      9829    127062
    ## [1653]     27042     64130      4793     23141      5728    283383      8382
    ## [1660]     57539 124900220    128977     54407    284021      3443     55816
    ## [1667]     10975     84285      6311      7450     93627      3823    348938
    ## [1674]     23040     51585      2940    340542     55920     54954     10396
    ## [1681]      6792      8986     25847      6314     27316    494115       650
    ## [1688]     51754     55296     54463     25980      4650     29922      9562
    ## [1695]      4904    286827     55723      6403     51575 100528030      8693
    ## [1702]    154661     64077     22850    169611     79966    404672      3704
    ## [1709]       920      5444    284434     92689     51614     83482      2742
    ## [1716]     54504    284805    440957    220323     51605       813      5930
    ## [1723]     56672      5158     54921 113455421    390082      5898    342908
    ## [1730]    200403    283093     55322     23504      9849      9295     57472
    ## [1737]    390637     57465      2185      5914    137970 100131980 124900627
    ## [1744]    728957 100132396    440077     10550    286016     57205     63895
    ## [1751]    284611      9675     25780        41    219333      9489      7414
    ## [1758]       381    284023     51759 105379252 100133036     22803     23313
    ## [1765]     54510     84332     54532     51337     23137     56144      9194
    ## [1772]     51004     55830     80177     51019      1510     55654     23229
    ## [1779]      3242     25764     10169      7128     51119      2811 101928906
    ## [1786] 102723678     57636      2555      2189      7994    117283      9414
    ## [1793]     57117    653784     79777    113451      8288    340533     84455
    ## [1800]     57380      6510      1020     11091      9590    645832     90268
    ## [1807]      3190      5798    266655    339512    166824     57192      4925
    ## [1814]      8632    339327     51733      4594      2110     80014      5328
    ## [1821]      2055     10440     93974      3843    202299    326625     26273
    ## [1828]     83734     84197     10528    692201     26793     84456     26499
    ## [1835]     23013      6638      8926      2593      4154    285971      5597
    ## [1842]      2962      5623       831      9667     92140       190     55288
    ## [1849]       640    257415      1948     23325     56261    134218      2893
    ## [1856]    149175     54681      6157     84078      8780     53335     29951
    ## [1863]      9748      8027     91544      6935    135295    123096      5890
    ## [1870]      1861    162073 100124535    132299    285855      1159    548596
    ## [1877]      3635     10949    221223      4495       123     55432     51528
    ## [1884]     64388      4140      3786     55347     55904     10234     10563
    ## [1891]      3821    255877      4116     92609     57732     93624     64397
    ## [1898]     10209     23348     80762     55713     84838      9790    221294
    ## [1905]     55282      9457    399668      4681 100532736      8128      7171
    ## [1912]       368    730013    653190     51667     55599      9061    152185
    ## [1919]      6387     84698     27183     26128     84259    285550     10055
    ## [1926]     11092      7410      3161    283932      5709      1404     54429
    ## [1933]      5797      2872       202     26524    151258    151987      3751
    ## [1940]      3767      9469     10347     54939      8573    158067     55735
    ## [1947]     10154    151246     55914      7019       182      6187     26784
    ## [1954]     51309     23243     84823    147687      2572     89927     81533
    ## [1961]      5288     10105     29851      5971    130399     29087    645332
    ## [1968]     23431     27079     57470     23338      8707     10882    157922
    ## [1975]     83594    259282      5332    340547     10307      2773     10645
    ## [1982]    284207    222546      7078       133     25948     55696     11116
    ## [1989]    170627     83982      4194     80833    116931     51176     79177
    ## [1996]      6988    138805     64801     11266    338949     84674    135228
    ## [2003]     10493     10894     54986     23197    148753      4437     57191
    ## [2010]     55213      5922    168448     85014     51385    169026     49855
    ## [2017]    284040 100533496      5923      6354     78988     49856      6461
    ## [2024] 100526842    497661      9982      6634     79717       407    196074
    ## [2031]    439934    201895     84135     84681    168455    221264     10730
    ## [2038]      6462      1611     23256     81614     55691     23164     57477
    ## [2045]      9203     51308     51109      9148      9702     10606      6830
    ## [2052]     84818     56681    285605      9467    347736     57188     57161
    ## [2059]        21    338321    196996    286097      5608     54471       468
    ## [2066]       571 100379661      8318 100033807    114803      1795     26189
    ## [2073]     26278     54439      3658      5264      1777    126393      7769
    ## [2080]    140458     55793     51752     24149    401262     57653 100499483
    ## [2087]    692084      8363      8362     55001      2113     84820    548644
    ## [2094]    246721     55758     81831    283377     23592     10623       523
    ## [2101]      3291     26580      9791    387522      7335     29099    137964
    ## [2108]    728621      1874      3444      3448     57684     84959      8850
    ## [2115]      5341     57522     51361     23507      6579      6415     92745
    ## [2122]      6738       329      2585    285456      6799      7327    643155
    ## [2129]     79783      6309     80829    619208      3903     51147    123016
    ## [2136]     11015     22941    285368     92949 100652871     25774      2941
    ## [2143]      5209     84066      6588     51382     57514      4321     23541
    ## [2150]     10549     23321     22801      9094      5546      7225     55013
    ## [2157]     51412      8804     51079    221687      6006     11191      4661
    ## [2164]      9092    118491     60493     79731     23613       205      1006
    ## [2171]       583     11280     23511    349152      6505       969     64284
    ## [2178]    150142      9360      2977     57701     23314      7046      6094
    ## [2185]     23530     94160     57150     63914     81849      7037    129401
    ## [2192]       379      9052     90835     55741 111089941    113419      9780
    ## [2199]     79036     27244    113444    340719     56834     84232     50486
    ## [2206]     79414      2026    401647     27071    158055    285335     55998
    ## [2213]     51705     51808      5948     83988     10628    401027      9986
    ## [2220]     57540     57698      4212    259197     84752      9666      4312
    ## [2227]     23195     11189    168433       868    285555      6628      7321
    ## [2234]       528    130752      7766     79750     27072      8204      7862
    ## [2241]    374659      1410      5825    134265    147841

``` r
edo <- enrichDGN(unique(filtered_genes$entrezgene_id))

barplot(edo, showCategory=25, font.size = 7) 
```

![](Gene-Set-Enrichment-Analysis_files/figure-gfm/unnamed-chunk-37-1.png)<!-- -->

``` r
# there seems to be a higher association with female diseases. T2D has a stronger influence on females?
```
