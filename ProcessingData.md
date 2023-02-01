Processing/Sorting data
================
Janet
31/01/2023

## Set up

Installing all the packages we will need for data analysis.

``` r
library(BiocManager)
BiocManager::install("limma")
```

    ## Bioconductor version 3.16 (BiocManager 1.30.19), R 4.2.1 (2022-06-23)

    ## Warning: package(s) not installed when version(s) same as or greater than current; use
    ##   `force = TRUE` to re-install: 'limma'

    ## Old packages: 'BiasedUrn', 'codetools', 'conflicted', 'future', 'nlme',
    ##   'sourcetools', 'tinytex', 'utf8', 'xfun'

Installing additional packages.

``` r
BiocManager :: install(c("edgeR", "Glimma", "org.Mm.eg.db", "RColorBrewer", "NMF", "BiasedUrn"))
```

    ## Bioconductor version 3.16 (BiocManager 1.30.19), R 4.2.1 (2022-06-23)

    ## Warning: package(s) not installed when version(s) same as or greater than current; use
    ##   `force = TRUE` to re-install: 'edgeR' 'Glimma' 'org.Mm.eg.db' 'RColorBrewer'
    ##   'NMF'

    ## Installing package(s) 'BiasedUrn'

    ## 
    ##   There is a binary version available but the source version is later:
    ##           binary source needs_compilation
    ## BiasedUrn  2.0.8  2.0.9              TRUE

    ## installing the source package 'BiasedUrn'

    ## Old packages: 'codetools', 'conflicted', 'future', 'nlme', 'sourcetools',
    ##   'tinytex', 'utf8', 'xfun'

Loading the packages.

``` r
library(limma)
```

    ## Warning: package 'limma' was built under R version 4.2.2

``` r
library(edgeR)
```

    ## Warning: package 'edgeR' was built under R version 4.2.2

``` r
library(Glimma)
library(org.Mm.eg.db)
```

    ## Loading required package: AnnotationDbi

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following object is masked from 'package:limma':
    ## 
    ##     plotMA

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
    ##     table, tapply, union, unique, unsplit, which.max, which.min

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: IRanges

    ## Loading required package: S4Vectors

    ## Warning: package 'S4Vectors' was built under R version 4.2.2

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## 

``` r
library(RColorBrewer)
library(NMF)
```

    ## Loading required package: registry

    ## Loading required package: rngtools

    ## Loading required package: cluster

    ## NMF - BioConductor layer [OK] | Shared memory capabilities [NO: bigmemory] | Cores 3/4

    ##   To enable shared memory capabilities, try: install.extras('
    ## NMF
    ## ')

    ## 
    ## Attaching package: 'NMF'

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     nrun

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2
    ## ──

    ## ✔ ggplot2 3.4.0     ✔ purrr   1.0.1
    ## ✔ tibble  3.1.8     ✔ dplyr   1.1.0
    ## ✔ tidyr   1.3.0     ✔ stringr 1.5.0
    ## ✔ readr   2.1.3     ✔ forcats 1.0.0
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::collapse()   masks IRanges::collapse()
    ## ✖ dplyr::combine()    masks Biobase::combine(), BiocGenerics::combine()
    ## ✖ dplyr::desc()       masks IRanges::desc()
    ## ✖ tidyr::expand()     masks S4Vectors::expand()
    ## ✖ dplyr::filter()     masks stats::filter()
    ## ✖ dplyr::first()      masks S4Vectors::first()
    ## ✖ dplyr::lag()        masks stats::lag()
    ## ✖ ggplot2::Position() masks BiocGenerics::Position(), base::Position()
    ## ✖ purrr::reduce()     masks IRanges::reduce()
    ## ✖ dplyr::rename()     masks S4Vectors::rename()
    ## ✖ dplyr::select()     masks AnnotationDbi::select()
    ## ✖ dplyr::slice()      masks IRanges::slice()

## Reading the data

Reading the raw data into R first.

``` r
wholeisletsdata <- read.table("pancreas_refseq_rpkms_wholeislet.txt", sep = "\t")

countsdata <- read.table("pancreas_refseq_rpkms_counts_3514sc.txt", sep = "\t")
```

These tables contain the RNA-sequencing data for single-cell and
whole-islet collected by the research group. Just want to explore the
raw data first.

``` r
head(wholeisletsdata)
```

    ##        V1                                  V2         V3        V4         V5
    ## 1   SGIP1                           NM_032291  1.8346344  2.903260  0.8923227
    ## 2   AZIN2              NM_052998+NM_001293562  4.3659185  3.596074  4.2254892
    ## 3   CLIC4                           NM_013943 19.8483416 13.735953 19.2217038
    ## 4   AGBL4                           NM_032785  0.9446949  1.147623  1.4361525
    ## 5  NECAP2 NM_001145277+NM_001145278+NM_018090 17.8139553 18.669569 15.5400754
    ## 6 SLC45A1                        NM_001080397  3.1712442  3.291781  2.3400335
    ##           V6         V7         V8        V9  V10  V11  V12  V13  V14  V15  V16
    ## 1  2.0496800  0.5858715  1.4801373  1.451944  360  616  204  567  171  270  268
    ## 2  3.2662948  5.5111072  5.0144054  3.118402  388  344  434  405  729  410  259
    ## 3 19.6130108 11.7552715 20.7124839 16.357313 3679 2753 4151 5125 3241 3569 2852
    ## 4  0.7780117  1.5554878  0.7320145  0.859605  118  155  209  137  289   85  101
    ## 5 22.5866481 14.5487809 23.0083211 19.210497 1522 1737 1562 2754 1870 1826 1572
    ## 6  2.3032176  2.3674838  2.7086801  1.871828  335  376  288  343  372  266  186
