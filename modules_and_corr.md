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

    ## pickSoftThreshold: will use block size 4297.
    ##  pickSoftThreshold: calculating connectivity for given powers...
    ##    ..working on genes 1 through 4297 of 4297
    ##    Power SFT.R.sq slope truncated.R.sq  mean.k. median.k.  max.k.
    ## 1      1   0.0605  7.60          0.921 2.17e+03  2.18e+03 2330.00
    ## 2      2   0.1910 -6.14          0.937 1.13e+03  1.12e+03 1330.00
    ## 3      3   0.5450 -7.45          0.947 5.95e+02  5.87e+02  807.00
    ## 4      4   0.7930 -6.91          0.977 3.22e+02  3.13e+02  511.00
    ## 5      5   0.8480 -5.77          0.981 1.77e+02  1.70e+02  336.00
    ## 6      6   0.8890 -4.79          0.986 9.98e+01  9.33e+01  228.00
    ## 7      7   0.9040 -4.14          0.984 5.74e+01  5.21e+01  159.00
    ## 8      8   0.9240 -3.59          0.978 3.37e+01  2.93e+01  113.00
    ## 9      9   0.9340 -3.19          0.971 2.02e+01  1.68e+01   82.80
    ## 10    10   0.9310 -2.89          0.958 1.24e+01  9.79e+00   62.20
    ## 11    12   0.9630 -2.34          0.966 4.98e+00  3.40e+00   37.00
    ## 12    14   0.9800 -1.97          0.975 2.20e+00  1.22e+00   23.40
    ## 13    16   0.9840 -1.74          0.979 1.06e+00  4.55e-01   15.80
    ## 14    18   0.9760 -1.61          0.970 5.56e-01  1.74e-01   11.50
    ## 15    20   0.9690 -1.52          0.962 3.15e-01  6.80e-02    8.50
    ## 16    22   0.9730 -1.42          0.967 1.90e-01  2.72e-02    6.39
    ## 17    24   0.9760 -1.35          0.973 1.21e-01  1.11e-02    4.86
    ## 18    26   0.9360 -1.35          0.936 7.98e-02  4.55e-03    3.81
    ## 19    28   0.9370 -1.32          0.950 5.46e-02  1.93e-03    3.10
    ## 20    30   0.8810 -1.37          0.920 3.84e-02  8.17e-04    2.55
    ## Soft power not provided. Automatically using the lowest power that meets 0.8 scale-free topology fit. Using soft_power = 5
    ##  Calculating consensus modules and module eigengenes block-wise from all genes
    ##  Calculating topological overlaps block-wise from all genes
    ##    Flagging genes and samples with too many missing values...
    ##     ..step 1
    ##     TOM calculation: adjacency..
    ##     ..will not use multithreading.
    ##      Fraction of slow calculations: 0.000000
    ##     ..connectivity..
    ##     ..matrix multiplication (system BLAS)..
    ##     ..normalization..
    ##     ..done.
    ##  ..Working on block 1 .
    ##  ..Working on block 1 .
    ##  ..merging consensus modules that are too close..
    ## [1] "yellow"
    ## [1] "grey"
    ## [1] "blue"
    ## [1] "brown"
    ## [1] "turquoise"
    ## [1] "red"
    ## [1] "green"
    ## pickSoftThreshold: will use block size 4297.
    ##  pickSoftThreshold: calculating connectivity for given powers...
    ##    ..working on genes 1 through 4297 of 4297
    ##    Power SFT.R.sq  slope truncated.R.sq  mean.k. median.k.  max.k.
    ## 1      1   0.0508  14.50          0.934 2.16e+03  2.16e+03 2220.00
    ## 2      2   0.3260 -21.70          0.981 1.10e+03  1.10e+03 1190.00
    ## 3      3   0.6240 -21.90          0.955 5.66e+02  5.64e+02  657.00
    ## 4      4   0.7280 -16.40          0.964 2.97e+02  2.94e+02  372.00
    ## 5      5   0.8200 -12.30          0.972 1.57e+02  1.55e+02  219.00
    ## 6      6   0.8880 -10.70          0.987 8.48e+01  8.28e+01  140.00
    ## 7      7   0.9310  -8.28          0.988 4.63e+01  4.46e+01   94.00
    ## 8      8   0.9550  -6.54          0.989 2.57e+01  2.43e+01   65.70
    ## 9      9   0.9690  -5.33          0.989 1.44e+01  1.33e+01   47.60
    ## 10    10   0.9770  -4.48          0.991 8.27e+00  7.39e+00   35.60
    ## 11    12   0.9810  -3.39          0.988 2.85e+00  2.33e+00   21.30
    ## 12    14   0.9770  -2.75          0.984 1.06e+00  7.57e-01   13.60
    ## 13    16   0.9830  -2.29          0.991 4.33e-01  2.53e-01    9.14
    ## 14    18   0.9730  -2.00          0.978 1.94e-01  8.66e-02    6.34
    ## 15    20   0.9740  -1.79          0.979 9.51e-02  3.05e-02    4.52
    ## 16    22   0.9690  -1.64          0.969 5.10e-02  1.09e-02    3.29
    ## 17    24   0.9720  -1.52          0.967 2.95e-02  4.03e-03    2.44
    ## 18    26   0.9860  -1.44          0.986 1.82e-02  1.50e-03    1.84
    ## 19    28   0.9650  -1.39          0.959 1.18e-02  5.72e-04    1.41
    ## 20    30   0.9570  -1.36          0.948 7.96e-03  2.23e-04    1.09
    ## Soft power not provided. Automatically using the lowest power that meets 0.8 scale-free topology fit. Using soft_power = 5
    ##  Calculating consensus modules and module eigengenes block-wise from all genes
    ##  Calculating topological overlaps block-wise from all genes
    ##    Flagging genes and samples with too many missing values...
    ##     ..step 1
    ##     TOM calculation: adjacency..
    ##     ..will not use multithreading.
    ##      Fraction of slow calculations: 0.000000
    ##     ..connectivity..
    ##     ..matrix multiplication (system BLAS)..
    ##     ..normalization..
    ##     ..done.
    ##  ..Working on block 1 .
    ##  ..Working on block 1 .
    ##  ..merging consensus modules that are too close..
    ## [1] "blue"
    ## [1] "turquoise"
    ## [1] "grey"
    ## [1] "red"
    ## [1] "brown"
    ## [1] "green"
    ## [1] "black"
    ## [1] "yellow"
    ## pickSoftThreshold: will use block size 4297.
    ##  pickSoftThreshold: calculating connectivity for given powers...
    ##    ..working on genes 1 through 4297 of 4297
    ##    Power SFT.R.sq  slope truncated.R.sq  mean.k. median.k.  max.k.
    ## 1      1   0.0062   8.31          0.998 2.15e+03  2.15e+03 2200.00
    ## 2      2   0.0130  -6.55          0.948 1.09e+03  1.09e+03 1140.00
    ## 3      3   0.0523  -8.50          0.914 5.56e+02  5.55e+02  602.00
    ## 4      4   0.2610 -13.60          0.943 2.87e+02  2.85e+02  325.00
    ## 5      5   0.5100 -12.10          0.924 1.49e+02  1.48e+02  178.00
    ## 6      6   0.7460 -11.70          0.923 7.82e+01  7.72e+01  107.00
    ## 7      7   0.8890  -9.49          0.945 4.15e+01  4.05e+01   66.70
    ## 8      8   0.9530  -7.58          0.957 2.22e+01  2.14e+01   43.90
    ## 9      9   0.9810  -6.16          0.978 1.20e+01  1.14e+01   30.40
    ## 10    10   0.9920  -5.08          0.990 6.61e+00  6.11e+00   22.60
    ## 11    12   0.9610  -3.78          0.950 2.08e+00  1.79e+00   14.10
    ## 12    14   0.9750  -2.88          0.967 7.06e-01  5.38e-01    9.89
    ## 13    16   0.4200  -3.20          0.273 2.63e-01  1.65e-01    7.48
    ## 14    18   0.4510  -2.77          0.369 1.10e-01  5.20e-02    5.91
    ## 15    20   0.4480  -2.40          0.379 5.26e-02  1.68e-02    4.81
    ## 16    22   0.4080  -2.59          0.367 2.85e-02  5.51e-03    4.01
    ## 17    24   0.4090  -2.39          0.381 1.72e-02  1.85e-03    3.39
    ## 18    26   0.4080  -2.54          0.369 1.13e-02  6.32e-04    2.91
    ## 19    28   0.3700  -2.00          0.210 7.86e-03  2.20e-04    2.52
    ## 20    30   0.3570  -2.16          0.181 5.74e-03  7.77e-05    2.21
    ## Soft power not provided. Automatically using the lowest power that meets 0.8 scale-free topology fit. Using soft_power = 7
    ##  Calculating consensus modules and module eigengenes block-wise from all genes
    ##  Calculating topological overlaps block-wise from all genes
    ##    Flagging genes and samples with too many missing values...
    ##     ..step 1
    ##     TOM calculation: adjacency..
    ##     ..will not use multithreading.
    ##      Fraction of slow calculations: 0.000000
    ##     ..connectivity..
    ##     ..matrix multiplication (system BLAS)..
    ##     ..normalization..
    ##     ..done.
    ##  ..Working on block 1 .
    ##  ..Working on block 1 .
    ##  ..merging consensus modules that are too close..
    ## [1] "green"
    ## [1] "blue"
    ## [1] "grey"
    ## [1] "turquoise"
    ## [1] "brown"
    ## [1] "black"
    ## [1] "red"
    ## [1] "yellow"
    ##   ..Excluding 1 genes from the calculation due to too many missing samples or zero variance.
    ## pickSoftThreshold: will use block size 4296.
    ##  pickSoftThreshold: calculating connectivity for given powers...
    ##    ..working on genes 1 through 4296 of 4296
    ##    Power SFT.R.sq  slope truncated.R.sq  mean.k. median.k.  max.k.
    ## 1      1   0.1200  33.90          0.974 2.16e+03  2.16e+03 2220.00
    ## 2      2   0.0792  15.20          0.995 1.10e+03  1.10e+03 1160.00
    ## 3      3   0.0368  -6.38          0.975 5.64e+02  5.63e+02  624.00
    ## 4      4   0.3640 -12.00          0.811 2.93e+02  2.92e+02  355.00
    ## 5      5   0.6130 -10.10          0.710 1.54e+02  1.52e+02  213.00
    ## 6      6   0.7780  -7.92          0.736 8.16e+01  8.03e+01  134.00
    ## 7      7   0.8770  -6.31          0.851 4.39e+01  4.26e+01   89.80
    ## 8      8   0.9010  -5.11          0.926 2.39e+01  2.28e+01   63.80
    ## 9      9   0.9080  -4.18          0.936 1.32e+01  1.23e+01   47.60
    ## 10    10   0.8990  -3.37          0.926 7.44e+00  6.69e+00   36.80
    ## 11    12   0.9110  -2.41          0.927 2.51e+00  2.01e+00   24.00
    ## 12    14   0.9320  -1.97          0.930 9.39e-01  6.25e-01   16.80
    ## 13    16   0.9430  -1.69          0.928 4.02e-01  1.99e-01   12.20
    ## 14    18   0.9620  -1.49          0.953 1.99e-01  6.50e-02    9.09
    ## 15    20   0.9720  -1.34          0.971 1.12e-01  2.17e-02    6.87
    ## 16    22   0.9640  -1.24          0.967 6.98e-02  7.42e-03    5.25
    ## 17    24   0.9560  -1.18          0.969 4.64e-02  2.59e-03    4.05
    ## 18    26   0.9620  -1.13          0.988 3.23e-02  9.22e-04    3.14
    ## 19    28   0.9580  -1.10          0.987 2.31e-02  3.33e-04    2.44
    ## 20    30   0.9380  -1.08          0.970 1.68e-02  1.23e-04    1.91
    ## Soft power not provided. Automatically using the lowest power that meets 0.8 scale-free topology fit. Using soft_power = 7
    ##  Calculating consensus modules and module eigengenes block-wise from all genes
    ##  Calculating topological overlaps block-wise from all genes
    ##    Flagging genes and samples with too many missing values...
    ##     ..step 1
    ##     TOM calculation: adjacency..
    ##     ..will not use multithreading.
    ##      Fraction of slow calculations: 0.000000
    ##     ..connectivity..
    ##     ..matrix multiplication (system BLAS)..
    ##     ..normalization..
    ##     ..done.
    ##  ..Working on block 1 .
    ##  ..Working on block 1 .
    ##  ..merging consensus modules that are too close..
    ## [1] "turquoise"
    ## [1] "grey"
    ## [1] "brown"
    ## [1] "red"
    ## [1] "green"
    ## [1] "black"
    ## [1] "blue"
    ## [1] "yellow"
    ## [1] "subsetting"
    ## [1] "subsetting"
    ## [1] "subsetting"
    ## [1] "subsetting"

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
    ##  [64] htmlwidgets_1.6.2      httr_1.4.7             FNN_1.1.3.2           
    ##  [67] gplots_3.1.3           RColorBrewer_1.1-3     ellipsis_0.3.2        
    ##  [70] ica_1.0-3              pkgconfig_2.0.3        nnet_7.3-19           
    ##  [73] uwot_0.1.16            deldir_1.0-9           utf8_1.2.4            
    ##  [76] tidyselect_1.2.0       rlang_1.1.3            reshape2_1.4.4        
    ##  [79] later_1.3.1            AnnotationDbi_1.58.0   munsell_0.5.0         
    ##  [82] tools_4.2.1            cachem_1.0.8           cli_3.6.2             
    ##  [85] generics_0.1.3         RSQLite_2.3.1          ggridges_0.5.4        
    ##  [88] evaluate_0.21          stringr_1.5.1          fastmap_1.1.1         
    ##  [91] yaml_2.3.7             goftest_1.2-3          knitr_1.43            
    ##  [94] bit64_4.0.5            fitdistrplus_1.1-11    caTools_1.18.2        
    ##  [97] purrr_1.0.2            RANN_2.6.1             KEGGREST_1.36.3       
    ## [100] pbapply_1.7-2          future_1.33.0          nlme_3.1-162          
    ## [103] mime_0.12              compiler_4.2.1         rstudioapi_0.14       
    ## [106] plotly_4.10.2          png_0.1-8              spatstat.utils_3.0-3  
    ## [109] tibble_3.2.1           stringi_1.8.3          lattice_0.21-8        
    ## [112] Matrix_1.5-4           vctrs_0.6.5            pillar_1.9.0          
    ## [115] lifecycle_1.0.4        spatstat.geom_3.2-4    lmtest_0.9-40         
    ## [118] RcppAnnoy_0.0.21       data.table_1.14.8      cowplot_1.1.1         
    ## [121] bitops_1.0-7           irlba_2.3.5.1          httpuv_1.6.11         
    ## [124] patchwork_1.1.3        R6_2.5.1               promises_1.2.1        
    ## [127] KernSmooth_2.23-20     gridExtra_2.3          IRanges_2.30.1        
    ## [130] parallelly_1.36.0      codetools_0.2-19       MASS_7.3-60           
    ## [133] gtools_3.9.4           withr_3.0.0            sctransform_0.3.5     
    ## [136] GenomeInfoDbData_1.2.8 S4Vectors_0.34.0       parallel_4.2.1        
    ## [139] rpart_4.1.19           grid_4.2.1             tidyr_1.3.0           
    ## [142] rmarkdown_2.24         Rtsne_0.16             spatstat.explore_3.2-1
    ## [145] base64enc_0.1-3        Biobase_2.56.0         shiny_1.7.5
