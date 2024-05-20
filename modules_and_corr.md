WGCNA & Age Correlation
================
Aaron Mohammed

``` r
library(hdWGCNA)
library(Seurat)
library(GeneOverlap)
library(WGCNA) 
library(dplyr)
library(writexl)
```

``` r
tissuez <- c("Germline", "Hypodermis", "Intestine", "Muscle", 
             "Neuron", "Pharynx", "Vulva_uterus")
```

``` r
modules.correlations.func <- function (tissue_name) {
  
  ############################################################################
                                # Subset Tissue
  ############################################################################
  
  tissue_dir <- file.path(main_path, tissue_name)
  dir.create(tissue_dir)
  
  # Strings that specify genotype and tissue
  N2_t_name <- paste('N2',sep= '_', tissue_name)
  LIPL4_t_name <- paste('LIPL4',sep= '_', tissue_name)
  DAF2_t_name <- paste('DAF2',sep= '_', tissue_name)
  RSKS1_t_name <- paste('RSKS1',sep= '_', tissue_name)
  
  # Subset the tissue
  Idents(seurat_data)<-'tissue'
  T_sub <- subset(seurat_data,idents = tissue_name)
  
  # Subset the samples
  Idents(T_sub) <- 'orig.ident'
  T_sub <- subset(T_sub,idents= c('N2D1', 'N2D1R', 'N2OPD1', 
                                  'N2D6', 'N2D12', 'N2D14', 
                                  'LIPL4D1', 'LIPL4D6', 'DAF2D1', 
                                  'DAF2D6', 'RSKS1D1', 'RSKS1D6'))
   
  # Add a column in meta.data that specifies the genotype
  strings <- T_sub@meta.data[['orig.ident']]
  strings <- gsub('D1.*','',strings)
  strings <- gsub('D6.*','',strings)
  strings <- gsub('OP','',strings)
  
  # unique(strings) # Output should be: 'LIPL4' 'N2' 'RSKS1' 'DAF2' 
  
  T_sub <- AddMetaData(T_sub, strings, col.name= 'geno')
  
  ############################################################################
                                    # WGCNA
  ############################################################################
  
  # Begin setup for WGCNA
  T_sub <- SetupForWGCNA(T_sub, 
                         gene_select = 'fraction', 
                         fraction = 0.05, 
                         wgcna_name = tissue_name)
  
  # Create metacells for each genotype
  T_sub <- MetacellsByGroups(seurat_obj = T_sub,
                             group.by = 'geno',
                             k = 25,
                             max_shared = 10, 
                             ident.group = 'geno')
  
  # Normalize the metacells
  T_sub <- NormalizeMetacells(T_sub)
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  wgcna.func <- function(seurat, group_name, group.by= 'geno', assay = 'RNA', 
                         slot = 'data', networkType = 'signed', geno_tiss, 
                         tom_outdir, overwrite_tom = TRUE) {
  
    # Setup expression matrix 
    S <- SetDatExpr(seurat,
                    group_name = group_name,
                    group.by= group.by, 
                    assay = assay, 
                    slot = slot)
    
    # Test different soft powers
    S <- TestSoftPowers(S,
                         networkType = networkType)
    
    # Construct the co-expression network
    S <- ConstructNetwork(S, 
                           tom_name = geno_tiss,
                           tom_outdir = tom_outdir,
                           overwrite_tom = overwrite_tom) 
    
    # Scale and center the data
    S <- ScaleData(S, features = VariableFeatures(S))
    
    # Compute the module eigengenes 
    S <- ModuleEigengenes(S)
    
    # Compute eigengene-based connectivity (kME)
    S <- ModuleConnectivity(S)
    
    return(S)
    
    }
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  
  # Run WGCNA for each genotype
  N2 <- wgcna.func(T_sub, 
                   group_name = 'N2', 
                   geno_tiss = N2_t_name, 
                   tom_outdir = tissue_dir)
  
  LIPL4 <- wgcna.func(T_sub, 
                      group_name = 'LIPL4', 
                      geno_tiss = LIPL4_t_name,
                      tom_outdir = tissue_dir)
  
  DAF2 <- wgcna.func(T_sub, 
                     group_name = 'DAF2', 
                     geno_tiss = DAF2_t_name,
                     tom_outdir = tissue_dir)
  
  RSKS1 <- wgcna.func(T_sub, 
                      group_name = 'RSKS1', 
                      geno_tiss = RSKS1_t_name,
                      tom_outdir = tissue_dir)

  ############################################################################
                                  # Age Correlation
  ############################################################################
  
  age <- T_sub@meta.data[['genotype']]
  age <- gsub('N2', '', age)
  age <- gsub('LIPL4', '', age)
  age <- gsub('DAF2', '', age)
  age <- gsub('RSKS1', '', age)
  age <- gsub('D', '', age)
  
  N2 <- AddMetaData(N2, as.numeric(age), col.name= 'age')
  LIPL4 <- AddMetaData(LIPL4, as.numeric(age), col.name= 'age')
  DAF2 <- AddMetaData(DAF2, as.numeric(age), col.name= 'age')
  RSKS1 <- AddMetaData(RSKS1, as.numeric(age), col.name= 'age')
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  age.corr <- function (seurat, ident = "geno", geno, traits) {
    
    Idents(seurat)<- ident
    
    # Calculate the correlation between the gene modules and age
    seurat <- ModuleTraitCorrelation(seurat,
                                     subset_by = ident,
                                     subset_groups = geno,
                                     traits = traits)
    
    # Extract age correlation and p-value results
    MT_corr <- GetModuleTraitCorrelation(seurat)
    corr_pval <- as.data.frame(MT_corr[['cor']][[geno]]['age',])
    corr_pval[, 2] <- MT_corr[['pval']][[geno]]['age',]
    colnames(corr_pval) <- c('corr','pval')
    
    # Concatenate the name of the genotype with the names of the module colors
    # and place these strings in a column of the data frame 
    corr_pval$geno_color <- paste0(geno, sep='_', rownames(corr_pval))
    
    # Filter the data frame keeping the rows with a p-value of 0
    corr_pval <- corr_pval[corr_pval$pval == 0, ]
    
    return(corr_pval)
    
  }
    
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  # Each row of the resulting data frames below include the modules that were
  # significantly correlated with age
  N2_CP <- age.corr(N2, geno = "N2", traits = c('age','nCount_RNA'))
  LIPL4_CP <- age.corr(LIPL4, geno = "LIPL4", traits = c('age','nCount_RNA'))
  DAF2_CP <- age.corr(DAF2, geno = "DAF2", traits = c('age','nCount_RNA'))
  RSKS1_CP <- age.corr(RSKS1, geno = "RSKS1", traits = c('age','nCount_RNA'))
  
  # Concatenate the geno_color columns of each data frame into a single vector
  sig.geno_modz <- c(N2_CP$geno_color, LIPL4_CP$geno_color, 
                     DAF2_CP$geno_color, RSKS1_CP$geno_color)
  
  ############################################################################
                    # Get Overlapping Genes Between Modules
  ############################################################################
  
  # Extract the gene module data from each genotype
  N2_modules <- GetModules(N2)
  LIPL4_modules <- GetModules(LIPL4)
  DAF2_modules <- GetModules(DAF2)
  RSKS1_modules <- GetModules(RSKS1)
  
  # Place each data set into a list
  all_modz_L <- list(N2_modules, LIPL4_modules, DAF2_modules, RSKS1_modules)
  names(all_modz_L) <- c('N2', 'LIPL4', 'DAF2', 'RSKS1')
  
  # Create an empty list that will be used to capture the significant module 
  # data of each genotype. Each significant module from each genotype will be
  # an index of this list
  sig_modz_L <- {}
  
  # Extract the significant gene module data
  for (i in 1:length(sig.geno_modz)) {
    geno <- gsub('_.*', '', sig.geno_modz[i])
    color <- gsub('.*_', '', sig.geno_modz[i])
    
    mod <- all_modz_L[[geno]][all_modz_L[[geno]]$module == color,] 
  
    sig_modz_L[[i]] <- mod
  }
    
  names(sig_modz_L) <- sig.geno_modz
  
  genome.size <- nrow(N2_modules) # This will be the same for each genotype
  overlapz <- {} # List for GeneOverlap results
  namez <- {} # List for the names of the comparisons ex. N2_yellow & LIPL4_blue
  index <- 0
  
  for (i in 1:length(sig_modz_L)) {
    A <- sig.geno_modz[i] # ex. N2_yellow
    
    for (j in 1:length(sig_modz_L)) {
      geno_A <- gsub('_.*','',sig.geno_modz[i]) # ex. N2
      geno_B <- gsub('_.*','',sig.geno_modz[j]) # ex. LIPL4
      
      if (geno_A != geno_B) { # ex. Don't run if comparing N2 & N2
      B <- sig.geno_modz[j] # ex. LIPL4_blue
      
      mod_1 <- sig_modz_L[[A]] # ex. Gene module data for N2_yellow
      mod_2 <- sig_modz_L[[B]] # ex. Gene module data for LIPL4_blue
      
      index = index + 1
      
      # Get the gene overlap between the two modules and store in the list
      overlapz[[index]] <- testGeneOverlap(newGeneOverlap(mod_1$gene_name,
                                                          mod_2$gene_name, 
                                                          genome.size = genome.size))
      
      # Keep track of the genotype modules being compared for this index
      namez[[index]] <- paste(A, sep='-', B)
      }
    }
  }
  
  # Add the names of the comparisons
  names(overlapz) <- namez
  
  # 
  overlap <- as.data.frame(matrix(nrow= length(overlapz), ncol= 6))
  names(overlap) <- c('M1','M2','M1_size','M2_size','Intersection_size','pval')
  overlap_names <- names(overlapz)
  
  for (i in 1:length(overlapz)) {
    name_mod_1 <- gsub('-.*','',overlap_names[i])
    name_mod_2 <- gsub('.*-','',overlap_names[i])
    
    overlap[i,'M1'] <- name_mod_1
    overlap[i,'M2'] <- name_mod_2
    overlap[i,'M1_size'] <- length(overlapz[[i]]@listA)
    overlap[i,'M2_size'] <- length(overlapz[[i]]@listB)
    overlap[i,'Intersection_size'] <- length(overlapz[[i]]@intersection)
    overlap[i,'pval'] <- overlapz[[i]]@pval
  }
  
  write.csv(overlap, 
            file.path(tissue_dir,
                      paste(tissue_name, 
                            sep= '_', 'overlap.csv')))
  
  # Get the genes that are a part of each module
  genez <- {}
  
  for (i in 1:length(sig_modz_L)) {
    genez[[i]] <- sig_modz_L[[i]][,1:2]
  }
  
  names(genez) <- names(sig_modz_L)
  
  write_xlsx(genez, file.path(tissue_dir,
                              paste(tissue_name, 
                                    sep= '_', 'gene_lists.xlsx')))
  
  # This is just for this markdown file to show an example of the result
  if (tissue_name == "Neuron") {return(overlap)}
  
}
```

``` r
for (i in 1:length(tissuez)) {modules.correlations.func(tissuez[i])}
```

## Below is an example of the gene overlap results between the modules of each genotype for neurons:

    ##                 M1              M2 M1_size M2_size Inter_size          pval
    ## 1        N2_yellow      LIPL4_blue     339     677        275 9.803770e-179
    ## 2        N2_yellow     LIPL4_brown     339     149         14  2.849248e-01
    ## 3        N2_yellow      DAF2_green     339     201         77  5.094222e-36
    ## 4        N2_yellow       DAF2_blue     339     310        220 1.780775e-203
    ## 5        N2_yellow      DAF2_brown     339     243          1  1.000000e+00
    ## 6        N2_yellow RSKS1_turquoise     339     601        267 1.227035e-183
    ## 7        N2_yellow     RSKS1_brown     339     330          0  1.000000e+00
    ## 8         N2_brown      LIPL4_blue     389     677         32  9.999985e-01
    ## 9         N2_brown     LIPL4_brown     389     149         57  2.960469e-23
    ## 10        N2_brown      DAF2_green     389     201         11  9.797214e-01
    ## 11        N2_brown       DAF2_blue     389     310          3  1.000000e+00
    ## 12        N2_brown      DAF2_brown     389     243          5  9.999990e-01
    ## 13        N2_brown RSKS1_turquoise     389     601         45  9.383061e-01
    ## 14        N2_brown     RSKS1_brown     389     330         13  9.999354e-01
    ## 15    N2_turquoise      LIPL4_blue     872     677         87  1.000000e+00
    ## 16    N2_turquoise     LIPL4_brown     872     149          6  1.000000e+00
    ## 17    N2_turquoise      DAF2_green     872     201         23  9.997597e-01
    ## 18    N2_turquoise       DAF2_blue     872     310         11  1.000000e+00
    ## 19    N2_turquoise      DAF2_brown     872     243         52  3.550733e-01
    ## 20    N2_turquoise RSKS1_turquoise     872     601         53  1.000000e+00
    ## 21    N2_turquoise     RSKS1_brown     872     330         98  1.490869e-05
    ## 22        N2_green      LIPL4_blue      94     677          0  1.000000e+00
    ## 23        N2_green     LIPL4_brown      94     149          0  1.000000e+00
    ## 24        N2_green      DAF2_green      94     201          0  1.000000e+00
    ## 25        N2_green       DAF2_blue      94     310          0  1.000000e+00
    ## 26        N2_green      DAF2_brown      94     243         13  2.167173e-03
    ## 27        N2_green RSKS1_turquoise      94     601          0  1.000000e+00
    ## 28        N2_green     RSKS1_brown      94     330          6  7.398978e-01
    ## 29      LIPL4_blue       N2_yellow     677     339        275 9.803770e-179
    ## 30      LIPL4_blue        N2_brown     677     389         32  9.999985e-01
    ## 31      LIPL4_blue    N2_turquoise     677     872         87  1.000000e+00
    ## 32      LIPL4_blue        N2_green     677      94          0  1.000000e+00
    ## 33      LIPL4_blue      DAF2_green     677     201        119  1.213525e-47
    ## 34      LIPL4_blue       DAF2_blue     677     310        239 7.626527e-143
    ## 35      LIPL4_blue      DAF2_brown     677     243         10  1.000000e+00
    ## 36      LIPL4_blue RSKS1_turquoise     677     601        364 2.959245e-174
    ## 37      LIPL4_blue     RSKS1_brown     677     330         13  1.000000e+00
    ## 38     LIPL4_brown       N2_yellow     149     339         14  2.849248e-01
    ## 39     LIPL4_brown        N2_brown     149     389         57  2.960469e-23
    ## 40     LIPL4_brown    N2_turquoise     149     872          6  1.000000e+00
    ## 41     LIPL4_brown        N2_green     149      94          0  1.000000e+00
    ## 42     LIPL4_brown      DAF2_green     149     201          7  5.517590e-01
    ## 43     LIPL4_brown       DAF2_blue     149     310         12  3.889256e-01
    ## 44     LIPL4_brown      DAF2_brown     149     243          2  9.984997e-01
    ## 45     LIPL4_brown RSKS1_turquoise     149     601         23  3.367095e-01
    ## 46     LIPL4_brown     RSKS1_brown     149     330          4  9.976080e-01
    ## 47      DAF2_green       N2_yellow     201     339         77  5.094222e-36
    ## 48      DAF2_green        N2_brown     201     389         11  9.797214e-01
    ## 49      DAF2_green    N2_turquoise     201     872         23  9.997597e-01
    ## 50      DAF2_green        N2_green     201      94          0  1.000000e+00
    ## 51      DAF2_green      LIPL4_blue     201     677        119  1.213525e-47
    ## 52      DAF2_green     LIPL4_brown     201     149          7  5.517590e-01
    ## 53      DAF2_green RSKS1_turquoise     201     601        134  2.154763e-70
    ## 54      DAF2_green     RSKS1_brown     201     330          0  1.000000e+00
    ## 55       DAF2_blue       N2_yellow     310     339        220 1.780775e-203
    ## 56       DAF2_blue        N2_brown     310     389          3  1.000000e+00
    ## 57       DAF2_blue    N2_turquoise     310     872         11  1.000000e+00
    ## 58       DAF2_blue        N2_green     310      94          0  1.000000e+00
    ## 59       DAF2_blue      LIPL4_blue     310     677        239 7.626527e-143
    ## 60       DAF2_blue     LIPL4_brown     310     149         12  3.889256e-01
    ## 61       DAF2_blue RSKS1_turquoise     310     601        233 1.093832e-148
    ## 62       DAF2_blue     RSKS1_brown     310     330          0  1.000000e+00
    ## 63      DAF2_brown       N2_yellow     243     339          1  1.000000e+00
    ## 64      DAF2_brown        N2_brown     243     389          5  9.999990e-01
    ## 65      DAF2_brown    N2_turquoise     243     872         52  3.550733e-01
    ## 66      DAF2_brown        N2_green     243      94         13  2.167173e-03
    ## 67      DAF2_brown      LIPL4_blue     243     677         10  1.000000e+00
    ## 68      DAF2_brown     LIPL4_brown     243     149          2  9.984997e-01
    ## 69      DAF2_brown RSKS1_turquoise     243     601          4  1.000000e+00
    ## 70      DAF2_brown     RSKS1_brown     243     330          5  9.999787e-01
    ## 71 RSKS1_turquoise       N2_yellow     601     339        267 1.227035e-183
    ## 72 RSKS1_turquoise        N2_brown     601     389         45  9.383061e-01
    ## 73 RSKS1_turquoise    N2_turquoise     601     872         53  1.000000e+00
    ## 74 RSKS1_turquoise        N2_green     601      94          0  1.000000e+00
    ## 75 RSKS1_turquoise      LIPL4_blue     601     677        364 2.959245e-174
    ## 76 RSKS1_turquoise     LIPL4_brown     601     149         23  3.367095e-01
    ## 77 RSKS1_turquoise      DAF2_green     601     201        134  2.154763e-70
    ## 78 RSKS1_turquoise       DAF2_blue     601     310        233 1.093832e-148
    ## 79 RSKS1_turquoise      DAF2_brown     601     243          4  1.000000e+00
    ## 80     RSKS1_brown       N2_yellow     330     339          0  1.000000e+00
    ## 81     RSKS1_brown        N2_brown     330     389         13  9.999354e-01
    ## 82     RSKS1_brown    N2_turquoise     330     872         98  1.490869e-05
    ## 83     RSKS1_brown        N2_green     330      94          6  7.398978e-01
    ## 84     RSKS1_brown      LIPL4_blue     330     677         13  1.000000e+00
    ## 85     RSKS1_brown     LIPL4_brown     330     149          4  9.976080e-01
    ## 86     RSKS1_brown      DAF2_green     330     201          0  1.000000e+00
    ## 87     RSKS1_brown       DAF2_blue     330     310          0  1.000000e+00
    ## 88     RSKS1_brown      DAF2_brown     330     243          5  9.999787e-01

``` r
sessionInfo()
```

    ## R version 4.2.1 (2022-06-23)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur ... 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] writexl_1.4.2         dplyr_1.1.4           WGCNA_1.72-1         
    ## [4] fastcluster_1.2.3     dynamicTreeCut_1.63-1 GeneOverlap_1.32.0   
    ## [7] SeuratObject_4.1.3    Seurat_4.3.0.1        hdWGCNA_0.2.03       
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] backports_1.4.1        Hmisc_5.0-1            plyr_1.8.8            
    ##   [4] igraph_1.5.1           lazyeval_0.2.2         sp_2.0-0              
    ##   [7] splines_4.2.1          listenv_0.9.0          scattermore_1.2       
    ##  [10] GenomeInfoDb_1.32.4    ggplot2_3.4.3          digest_0.6.33         
    ##  [13] foreach_1.5.2          htmltools_0.5.6        GO.db_3.15.0          
    ##  [16] fansi_1.0.6            checkmate_2.2.0        magrittr_2.0.3        
    ##  [19] memoise_2.0.1          tensor_1.5             cluster_2.1.4         
    ##  [22] doParallel_1.0.17      ROCR_1.0-11            globals_0.16.2        
    ##  [25] Biostrings_2.64.1      matrixStats_1.0.0      spatstat.sparse_3.0-2 
    ##  [28] colorspace_2.1-0       blob_1.2.4             ggrepel_0.9.3         
    ##  [31] xfun_0.40              RCurl_1.98-1.12        crayon_1.5.2          
    ##  [34] jsonlite_1.8.8         progressr_0.14.0       spatstat.data_3.0-1   
    ##  [37] impute_1.70.0          survival_3.5-7         zoo_1.8-12            
    ##  [40] iterators_1.0.14       glue_1.7.0             polyclip_1.10-4       
    ##  [43] gtable_0.3.4           zlibbioc_1.42.0        XVector_0.36.0        
    ##  [46] leiden_0.4.3           future.apply_1.11.0    BiocGenerics_0.42.0   
    ##  [49] abind_1.4-5            scales_1.2.1           DBI_1.1.3             
    ##  [52] spatstat.random_3.1-5  miniUI_0.1.1.1         Rcpp_1.0.11           
    ##  [55] htmlTable_2.4.1        viridisLite_0.4.2      xtable_1.8-4          
    ##  [58] reticulate_1.35.0      foreign_0.8-85         bit_4.0.5             
    ##  [61] preprocessCore_1.58.0  Formula_1.2-5          stats4_4.2.1          
    ##  [64] htmlwidgets_1.6.2      httr_1.4.7             gplots_3.1.3          
    ##  [67] RColorBrewer_1.1-3     ellipsis_0.3.2         ica_1.0-3             
    ##  [70] pkgconfig_2.0.3        nnet_7.3-19            uwot_0.1.16           
    ##  [73] deldir_1.0-9           utf8_1.2.4             tidyselect_1.2.0      
    ##  [76] rlang_1.1.3            reshape2_1.4.4         later_1.3.1           
    ##  [79] AnnotationDbi_1.58.0   munsell_0.5.0          tools_4.2.1           
    ##  [82] cachem_1.0.8           cli_3.6.2              generics_0.1.3        
    ##  [85] RSQLite_2.3.1          ggridges_0.5.4         evaluate_0.21         
    ##  [88] stringr_1.5.1          fastmap_1.1.1          yaml_2.3.7            
    ##  [91] goftest_1.2-3          knitr_1.43             bit64_4.0.5           
    ##  [94] fitdistrplus_1.1-11    caTools_1.18.2         purrr_1.0.2           
    ##  [97] RANN_2.6.1             KEGGREST_1.36.3        pbapply_1.7-2         
    ## [100] future_1.33.0          nlme_3.1-162           mime_0.12             
    ## [103] compiler_4.2.1         rstudioapi_0.14        plotly_4.10.2         
    ## [106] png_0.1-8              spatstat.utils_3.0-3   tibble_3.2.1          
    ## [109] stringi_1.8.3          lattice_0.21-8         Matrix_1.5-4          
    ## [112] vctrs_0.6.5            pillar_1.9.0           lifecycle_1.0.4       
    ## [115] spatstat.geom_3.2-4    lmtest_0.9-40          RcppAnnoy_0.0.21      
    ## [118] data.table_1.14.8      cowplot_1.1.1          bitops_1.0-7          
    ## [121] irlba_2.3.5.1          httpuv_1.6.11          patchwork_1.1.3       
    ## [124] R6_2.5.1               promises_1.2.1         KernSmooth_2.23-20    
    ## [127] gridExtra_2.3          IRanges_2.30.1         parallelly_1.36.0     
    ## [130] codetools_0.2-19       MASS_7.3-60            gtools_3.9.4          
    ## [133] sctransform_0.3.5      GenomeInfoDbData_1.2.8 S4Vectors_0.34.0      
    ## [136] parallel_4.2.1         rpart_4.1.19           grid_4.2.1            
    ## [139] tidyr_1.3.0            rmarkdown_2.24         Rtsne_0.16            
    ## [142] spatstat.explore_3.2-1 base64enc_0.1-3        Biobase_2.56.0        
    ## [145] shiny_1.7.5
