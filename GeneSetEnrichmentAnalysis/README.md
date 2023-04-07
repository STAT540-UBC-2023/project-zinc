# Introduction 

This directory contains the R markdown file and the markdown file generated from the Gene Set Enrichment Analysis(GSEA) of the differentially expressed genes. 

We analyzed a total of 4 pathways from [MSigDB](http://www.gsea-msigdb.org/gsea/msigdb/human/collection_details.jsp#IMMUNESIGDB) database. And the following pathways are being selected and analyzed: 
1. [KEGG](http://www.pathway.jp) Pathway
2. [REACTOME](http://www.reactome.org) Pathway 
3. [Pathway Interaction Database](http://www.ndexbio.org)(PID) Pathway 
4. [ImmuneSigDB](https://www.cell.com/immunity/fulltext/S1074-7613(15)00532-4) Pathway 

We also used the [GWAS Catalog](https://www.ebi.ac.uk/gwas/) to identify genes associated with T2D. And [DisGeNET](https://www.disgenet.org/) was used for other disease found to be associated with T2D. 

The ['degT2d_DEcopy.RDS'](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/GeneSetEnrichmentAnalysis/degT2d_DEcopy.RDS) contains the list of differentially expressed genes in samples with Type-2 diabetes, and it was used for GSEA analysis. 

The [Result](https://github.com/STAT540-UBC-2023/project-zinc/tree/main/GeneSetEnrichmentAnalysis/Result) sub-directory contains all the significant pathways (FDR<0.05) identified and genes found to be associated with T2D on GWAS catalog. 
  + [KEGG](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/GeneSetEnrichmentAnalysis/Result/kegg_significant_pathways.csv) 
  + [REACTOME](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/GeneSetEnrichmentAnalysis/Result/reactome_significant_pathways.csv)
  + [PID](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/GeneSetEnrichmentAnalysis/Result/pid_significant_pathways.csv) 
  + [ImmuneSigDB](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/GeneSetEnrichmentAnalysis/Result/IMMUNE_significant_pathways.csv)
  + [GWAS](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/GeneSetEnrichmentAnalysis/Result/tsd.genes.csv)


The [Figure](https://github.com/STAT540-UBC-2023/project-zinc/tree/main/GeneSetEnrichmentAnalysis/Gene-Set-Enrichment-Analysis_files/figure-gfm) sub-directory contains all the plots generated. 
  1. Barplots: 
  - [KEGG](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/GeneSetEnrichmentAnalysis/Gene-Set-Enrichment-Analysis_files/figure-gfm/unnamed-chunk-20-1.png)
  - [REACTOME](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/GeneSetEnrichmentAnalysis/Gene-Set-Enrichment-Analysis_files/figure-gfm/unnamed-chunk-30-1.png)
  - [PID](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/GeneSetEnrichmentAnalysis/Gene-Set-Enrichment-Analysis_files/figure-gfm/unnamed-chunk-24-1.png)
  - [ImmuneSigDB](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/GeneSetEnrichmentAnalysis/Gene-Set-Enrichment-Analysis_files/figure-gfm/unnamed-chunk-36-1.png) 
  - [DisGeNet](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/GeneSetEnrichmentAnalysis/Gene-Set-Enrichment-Analysis_files/figure-gfm/unnamed-chunk-38-1.png) 
  2. Additional visualizations for KEGG analysis results: 
  - [Enrichment plot](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/GeneSetEnrichmentAnalysis/Gene-Set-Enrichment-Analysis_files/figure-gfm/unnamed-chunk-21-1.png) 
  - [Ridge plot](https://github.com/STAT540-UBC-2023/project-zinc/blob/main/GeneSetEnrichmentAnalysis/Gene-Set-Enrichment-Analysis_files/figure-gfm/unnamed-chunk-21-2.png)
    
