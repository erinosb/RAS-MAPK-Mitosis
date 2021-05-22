# Analysis: Heatmaps of key gene signature sets

The purpose of this analysis is to assess how and whether published gene
signature sets are changing the Pollak et al dataset. To do this, we
will calculate the fold-enrichment of Here is the Pollak et al dataset:

-   Pollak J, Rai KG, Funk CC, Arora S, Lee E, Zhu J, Price ND, Paddison
    PJ, Ramirez JM, Rostomily RC. [Ion channel expression patterns in
    glioblastoma stem cells with functional and therapeutic implications
    for malignancy](https://doi.org/10.1371/journal.pone.0172884) *PLoS
    One* **2017**; Mar 6;12(3):e0172884.

For this analysis, we used the following publications:

-   Creighton CJ. [A gene transcription signature of the Akt/mTOR
    pathway in clinical breast
    tumors](https://doi.org/10.1038/sj.onc.1210245). *Oncogene*
    **2007**;26:4648–55.

-   Dry JR, Pavey S, Pratilas CA, Harbron C, Runswick S, Hodgson D, et
    al. [Transcriptional Pathway Signatures Predict MEK Addiction and
    Response to Selumetinib
    (AZD6244)](https://doi.org/10.1158/0008-5472.can-09-1577). *Cancer
    Res* **2010**;70:2264–73.

-   Loboda A, Nebozhyn M, Klinghoffer R, Frazier J, Chastain M, Arthur
    W, et al. [A gene expression signature of RAS pathway dependence
    predicts response to PI3K and RAS pathway inhibitors and expands the
    population of RAS pathway activated
    tumors.](https://doi.org/10.1186/1755-8794-3-26) *Bmc Med Genomics*
    **2010**;3:26.

-   Pratilas CA, Taylor BS, Ye Q, Viale A, Sander C, Solit DB, et
    al. [V600EBRAF is associated with disabled feedback inhibition of
    RAF–MEK signaling and elevated transcriptional output of the
    pathway.](https://doi.org/10.1073/pnas.0900780106) *Proc National
    Acad Sci* **2009**;106:4519–24.

-   Wagle M-C, Kirouac D, Klijn C, Liu B, Mahajan S, Junttila M, et
    al. [A transcriptional MAPK Pathway Activity Score (MPAS) is a
    clinically relevant biomarker in multiple cancer
    types](https://doi.org/10.1038/s41698-018-0051-4). *Npj Precis
    Oncol* **2018**;2:7.

## Step 1: Load packages

``` r
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(reshape2)
library(pryr)
library(tinytex)
```

## Step 2: Set Today’s date:

    ## [1] "210522"

## Step 3: Import

Import the lists of changing genes from the Pollak paper generated in
the DESeq2 analysis.

## Step 4: Import the datasets from other published papers -

## Step 5: Writing functions 1 -

Write a function called *make_contingency_table* that creates
contingency tables for statistical analyses:

Write a function called *vector_of_p\_values* that creates a vector of
p-values generated from comparisons with the DESeq2 differentially
expressed genes expressing the likeliness that a given gene signature
set of genes are **enriched** amongst the list of genes.

Write a function called *vector_of_depletion_p\_values* that creates a
vector of p=values generated from comparisons with the DESeq2
differentially expressed genes expressing the likeliness that a given
gene signature set of genes are **depleted** amongst the list of genes.

Make contingency table creates contingency tables like so:

``` r
exampletable <- make_contingency_table(lobodaUp_ensbl_geneID , upIn1502_2$EnsembleGeneName, all_genes$ensembl_gene_id)
exampletable
```

    ##                 published Notpublished
    ## RNAseqChanging         14          523
    ## RNAseqUnchanged        91        14873

P values for enrichment are calculated from these contingency tables
using Fisher exact test:

``` r
examplePvalue <- fisher.test(exampletable, alternative = "greater")$p.value
examplePvalue
```

    ## [1] 1.533953e-05

## Step 6: Writing functions 2 -

Write a function called *calculate_percent_enrichment* that calculates
the over-representation of a published gene signature set of genes
within the gene sets from GME RNA-seq analysis. The gene sets lists here
are up-regulated genes published in Pollak et al. that we determined
were over-expressed in each GME cell line.

Write a function called *calculate_fold_enrichment* that calcualets the
fold enrichment between a published genes signature set and the genesets
from our RNA-seq analysis.

Percent enrichment is calculated as:
$$\\text{Percent enrichment}= \\displaystyle \\left(\\frac{\\text{set and list intersection}}{\\text{GME list length}}\\right)100\\%$$

Fold enrichemnt is the percent calculation above divided by a similar
calculation performed for the genome as a whole (instead of \# in GME
list)

Fold enrichment calculated as:
$$\\text{Fold enrichment}= \\displaystyle \\frac{\\text{percent enrichment in GME}}{\\text{percent enrichment genomewide}}$$

``` r
exampletable
```

    ##                 published Notpublished
    ## RNAseqChanging         14          523
    ## RNAseqUnchanged        91        14873

``` r
examplePercent <- sum(exampletable[1,1])/sum(exampletable[1,])*100
examplePercent
```

    ## [1] 2.607076

``` r
percentGenomewide <- sum(exampletable[,1])/sum(exampletable)*100
percentGenomewide
```

    ## [1] 0.6773757

``` r
foldEnrichmentExample <- examplePercent/ percentGenomewide

foldEnrichmentExample
```

    ## [1] 3.84879


    ## Step 7: Gather statistics


    ```r
    ######################
    # P-values
    ######################

    # Generate an empty dataframe for the p-values calculating statistical significance of enrichment:
    Pvalue_df <- data.frame(UpIn827 = numeric(),
                            UpIn1502 = numeric(),
                            UpIn131 = numeric(),
                            UpInG166 = numeric(),
                            UpInG179 = numeric(),
                            DownIn827 = numeric(),
                            DownIn1502 = numeric(),
                            DownIn131 = numeric(),
                            DownInG166 = numeric(),
                            DownInG179 = numeric())

    # Save P-values of statistical enrichment to a dataframe called Pvalue_df
    Pvalue_df[1,] <- vector_of_p_values(mek18_ensbl_geneID)
    Pvalue_df[2,] <- vector_of_p_values(Pratilas_ensbl_geneID)
    Pvalue_df[3,] <- vector_of_p_values(wagle_ensbl_geneID)
    Pvalue_df[4,] <- vector_of_p_values(creighton_ensbl_geneID)
    Pvalue_df[5,] <- vector_of_p_values(lobodaUp_ensbl_geneID)

    rownames(Pvalue_df)[1] <- "Dry_et_al_mek_18"
    rownames(Pvalue_df)[2] <- "Pratilas_et_al_mek_50"
    rownames(Pvalue_df)[3] <- "Wagle_et_al_mek_10"
    rownames(Pvalue_df)[4] <- "Creighton_et_al_akt_57"
    rownames(Pvalue_df)[5] <- "Loboda_ras_101"

    #Pvalue_df

    ######################
    # Percent
    ######################

    # Generate an empty dataframe for the percentages:
    percent_df <- data.frame(UpIn827 = numeric(),
                            UpIn1502 = numeric(),
                            UpIn131 = numeric(),
                            UpInG166 = numeric(),
                            UpInG179 = numeric(),
                            DownIn827 = numeric(),
                            DownIn1502 = numeric(),
                            DownIn131 = numeric(),
                            DownInG166 = numeric(),
                            DownInG179 = numeric(),
                            Genomewide = numeric())

    # Save percents to a dataframe called percent_df
    percent_df[1,] <- calculate_percent_enrichment(mek18_ensbl_geneID)
    percent_df[2,] <- calculate_percent_enrichment(Pratilas_ensbl_geneID)
    percent_df[3,] <- calculate_percent_enrichment(wagle_ensbl_geneID)
    percent_df[4,] <- calculate_percent_enrichment(creighton_ensbl_geneID)
    percent_df[5,] <- calculate_percent_enrichment(lobodaUp_ensbl_geneID)

    rownames(percent_df)[1] <- "Dry_et_al_mek_18"
    rownames(percent_df)[2] <- "Pratilas_et_al_mek_50"
    rownames(percent_df)[3] <- "Wagle_et_al_mek_10"
    rownames(percent_df)[4] <- "Creighton_et_al_akt_57"
    rownames(percent_df)[5] <- "Loboda_ras_101"

    #percent_df

    ######################
    # Fold Enrichment
    ######################

    # Generate an empty dataframe for the fold enrichment:
    df_enrich <- data.frame(UpIn827 = numeric(),
                            UpIn1502 = numeric(),
                            UpIn131 = numeric(),
                            UpInG166 = numeric(),
                            UpInG179 = numeric(),
                            DownIn827 = numeric(),
                            DownIn1502 = numeric(),
                            DownIn131 = numeric(),
                            DownInG166 = numeric(),
                            DownInG179 = numeric())

    # Gather statistics: Fold enrichment - enrichment each list
    df_enrich[1,] <- calculate_fold_enrichment(mek18_ensbl_geneID)
    df_enrich[2,] <- calculate_fold_enrichment(Pratilas_ensbl_geneID)
    df_enrich[3,] <- calculate_fold_enrichment(wagle_ensbl_geneID)
    df_enrich[4,] <- calculate_fold_enrichment(creighton_ensbl_geneID)
    df_enrich[5,] <- calculate_fold_enrichment(lobodaUp_ensbl_geneID)

    rownames(df_enrich)[1] <- "Dry_et_al_mek_18"
    rownames(df_enrich)[2] <- "Pratilas_et_al_mek_50"
    rownames(df_enrich)[3] <- "Wagle_et_al_mek_10"
    rownames(df_enrich)[4] <- "Creighton_et_al_akt_57"
    rownames(df_enrich)[5] <- "Loboda_ras_101"

    #df_enrich

    ######################
    # Fold Depletion
    ######################

    # Generate an empty dataframe for the fold depletion:
    df_deplete <- data.frame(UpIn827 = numeric(),
                            UpIn1502 = numeric(),
                            UpIn131 = numeric(),
                            UpInG166 = numeric(),
                            UpInG179 = numeric(),
                            DownIn827 = numeric(),
                            DownIn1502 = numeric(),
                            DownIn131 = numeric(),
                            DownInG166 = numeric(),
                            DownInG179 = numeric())

    # Gather statistics: Fold depletion - depletion in each list
    df_deplete[1,] <- calculate_fold_depletion(mek18_ensbl_geneID)
    df_deplete[2,] <- calculate_fold_depletion(Pratilas_ensbl_geneID)
    df_deplete[3,] <- calculate_fold_depletion(wagle_ensbl_geneID)
    df_deplete[4,] <- calculate_fold_depletion(creighton_ensbl_geneID)
    df_deplete[5,] <- calculate_fold_depletion(lobodaUp_ensbl_geneID)

    rownames(df_deplete)[1] <- "Dry_et_al_mek_18"
    rownames(df_deplete)[2] <- "Pratilas_et_al_mek_50"
    rownames(df_deplete)[3] <- "Wagle_et_al_mek_10"
    rownames(df_deplete)[4] <- "Creighton_et_al_akt_57"
    rownames(df_deplete)[5] <- "Loboda_ras_101"

    #df_deplete

## Step 8: Barplots of enrichment/Depletion

![](210516_HeatmapsOfRelevantGenes_Pollak_files/figure-markdown_github/unnamed-chunk-10-1.png)![](210516_HeatmapsOfRelevantGenes_Pollak_files/figure-markdown_github/unnamed-chunk-10-2.png)![](210516_HeatmapsOfRelevantGenes_Pollak_files/figure-markdown_github/unnamed-chunk-10-3.png)![](210516_HeatmapsOfRelevantGenes_Pollak_files/figure-markdown_github/unnamed-chunk-10-4.png)![](210516_HeatmapsOfRelevantGenes_Pollak_files/figure-markdown_github/unnamed-chunk-10-5.png)

## Step 9: Dot plot:

    ## No id variables; using all as measure variables

![](210516_HeatmapsOfRelevantGenes_Pollak_files/figure-markdown_github/unnamed-chunk-11-1.png)

    ## quartz_off_screen 
    ##                 2

    ## [1] "Here is a matrix of P-values"

    ##                          UpIn827     UpIn1502      UpIn131    UpInG166
    ## Dry_et_al_mek_18       0.1383638 4.700563e-01 4.256552e-01 0.105334819
    ## Pratilas_et_al_mek_50  0.1206343 5.421819e-01 7.988593e-01 0.480668917
    ## Wagle_et_al_mek_10     1.0000000 1.000000e+00 1.000000e+00 1.000000000
    ## Creighton_et_al_akt_57 1.0000000 8.664615e-01 8.276531e-01 0.529602246
    ## Loboda_ras_101         0.7418419 1.533953e-05 1.627833e-08 0.005206878
    ##                            UpInG179 DownIn827   DownIn1502    DownIn131
    ## Dry_et_al_mek_18       0.5004543314 0.1383638 4.700563e-01 4.256552e-01
    ## Pratilas_et_al_mek_50  0.8656489018 0.1206343 5.421819e-01 7.988593e-01
    ## Wagle_et_al_mek_10     1.0000000000 1.0000000 1.000000e+00 1.000000e+00
    ## Creighton_et_al_akt_57 1.0000000000 1.0000000 8.664615e-01 8.276531e-01
    ## Loboda_ras_101         0.0001621189 0.7418419 1.533953e-05 1.627833e-08
    ##                         DownInG166   DownInG179
    ## Dry_et_al_mek_18       0.105334819 0.5004543314
    ## Pratilas_et_al_mek_50  0.480668917 0.8656489018
    ## Wagle_et_al_mek_10     1.000000000 1.0000000000
    ## Creighton_et_al_akt_57 0.529602246 1.0000000000
    ## Loboda_ras_101         0.005206878 0.0001621189

    ## [1] "Here is a matrix of negative log-transformed P-values"

    ##                          UpIn827   UpIn1502    UpIn131  UpInG166  UpInG179
    ## Creighton_et_al_akt_57 0.0000000  0.1433376  0.1891611 0.6356290 0.0000000
    ## Pratilas_et_al_mek_50  2.1149913  0.6121537  0.2245704 0.7325766 0.1442759
    ## Dry_et_al_mek_18       1.9778691  0.7549028  0.8541257 2.2506113 0.6922389
    ## Wagle_et_al_mek_10     0.0000000  0.0000000  0.0000000 0.0000000 0.0000000
    ## Loboda_ras_101         0.2986192 11.0850774 17.9334313 5.2577749 8.7271804

    ## [1] "Here is a matrix of fold-enrichment values"

    ##                          UpIn827  UpIn1502   UpIn131 UpInG166  UpInG179
    ## Creighton_et_al_akt_57 0.0000000 0.5064197 0.5786114 1.135480 0.0000000
    ## Pratilas_et_al_mek_50  2.1104152 1.1102278 0.6342471 1.244660 0.5086965
    ## Dry_et_al_mek_18       3.0483776 1.6036623 1.8322695 3.595685 1.4695677
    ## Wagle_et_al_mek_10     0.0000000 0.0000000 0.0000000 0.000000 0.0000000
    ## Loboda_ras_101         0.7838685 3.8487896 5.3397568 2.773814 3.2750366

    ## [1] "Here is a matrix of log-transformed fold-enrichment values"

    ##                           UpIn827   UpIn1502    UpIn131  UpInG166   UpInG179
    ## Creighton_et_al_akt_57       -Inf -0.6803895 -0.5471241 0.1270551       -Inf
    ## Pratilas_et_al_mek_50   0.7468847  0.1045652 -0.4553166 0.2188627 -0.6759037
    ## Dry_et_al_mek_18        1.1146095  0.4722900  0.6055554 1.2797346  0.3849683
    ## Wagle_et_al_mek_10           -Inf       -Inf       -Inf      -Inf       -Inf
    ## Loboda_ras_101         -0.2435140  1.3477587  1.6751801 1.0202234  1.1863290

    ## [1] "Here is a matrix of fold-depletion values"

    ##                          UpIn827 UpIn1502   UpIn131  UpInG166  UpInG179
    ## Dry_et_al_mek_18       1.3038102 0.000000 0.0000000 1.1700634 4.0054264
    ## Pratilas_et_al_mek_50  0.2256595 0.000000 0.8178221 0.8100439 0.6932469
    ## Wagle_et_al_mek_10     1.1734292 3.263368 0.0000000 0.0000000 0.0000000
    ## Creighton_et_al_akt_57 0.2058648 0.000000 0.0000000 0.0000000 0.0000000
    ## Loboda_ras_101         1.8998378 0.621594 0.8100333 0.8023292 1.0299668
    ##                        DownIn827 DownIn1502 DownIn131 DownInG166 DownInG179
    ## Dry_et_al_mek_18       1.3038102   0.000000 0.0000000  1.1700634          1
    ## Pratilas_et_al_mek_50  0.2256595   0.000000 0.8178221  0.8100439          1
    ## Wagle_et_al_mek_10     1.1734292   3.263368 0.0000000  0.0000000        NaN
    ## Creighton_et_al_akt_57 0.2058648   0.000000 0.0000000  0.0000000        NaN
    ## Loboda_ras_101         1.8998378   0.621594 0.8100333  0.8023292          1

Session info:

``` r
sessionInfo()
```

    ## R version 4.0.3 (2020-10-10)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Catalina 10.15.7
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] tinytex_0.31       pryr_0.1.4         reshape2_1.4.4     pheatmap_1.0.12   
    ## [5] dplyr_1.0.3        ggplot2_3.3.3      RColorBrewer_1.1-2
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.6        highr_0.8         pillar_1.4.7      compiler_4.0.3   
    ##  [5] plyr_1.8.6        tools_4.0.3       digest_0.6.27     evaluate_0.14    
    ##  [9] lifecycle_0.2.0   tibble_3.0.6      gtable_0.3.0      pkgconfig_2.0.3  
    ## [13] rlang_0.4.10      DBI_1.1.1         yaml_2.2.1        xfun_0.20        
    ## [17] withr_2.4.1       stringr_1.4.0     knitr_1.31        generics_0.1.0   
    ## [21] vctrs_0.3.6       grid_4.0.3        tidyselect_1.1.0  glue_1.4.2       
    ## [25] R6_2.5.0          rmarkdown_2.6     farver_2.0.3      purrr_0.3.4      
    ## [29] magrittr_2.0.1    codetools_0.2-18  scales_1.1.1      ellipsis_0.3.1   
    ## [33] htmltools_0.5.1.1 assertthat_0.2.1  colorspace_2.0-0  labeling_0.4.2   
    ## [37] stringi_1.5.3     munsell_0.5.0     crayon_1.4.0

``` r
citation()
```

    ## 
    ## To cite R in publications use:
    ## 
    ##   R Core Team (2020). R: A language and environment for statistical
    ##   computing. R Foundation for Statistical Computing, Vienna, Austria.
    ##   URL https://www.R-project.org/.
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Manual{,
    ##     title = {R: A Language and Environment for Statistical Computing},
    ##     author = {{R Core Team}},
    ##     organization = {R Foundation for Statistical Computing},
    ##     address = {Vienna, Austria},
    ##     year = {2020},
    ##     url = {https://www.R-project.org/},
    ##   }
    ## 
    ## We have invested a lot of time and effort in creating R, please cite it
    ## when using it for data analysis. See also 'citation("pkgname")' for
    ## citing R packages.
