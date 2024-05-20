Module Age Correlation & Circos Plots
================
Aaron Mohammed

``` r
library(hdWGCNA)
library(Seurat)
library(GeneOverlap)
library(WGCNA) 
library(dplyr)
library(writexl)
library(circlize)
library(tidyr)

modules_path <- file.path(getwd(), "Modules")
dir.create(modules_path)

plots_path <- file.path(modules_path, 'Circos_plots')
dir.create(plots_path)
```

``` r
tissuez <- c("Germline", "Hypodermis", "Intestine", "Muscle", 
             "Neuron", "Pharynx", "Vulva_uterus")
```

<p align="center">
<h1>
WGCNA & Module Trait Correlation
</h1>
</p>
<p style="font-size: 50px">
In this section, I wrote a function that does the following:
<ul>
<li>
Subsets the full seurat object by tissue
</li>
<li>
Performs WGCNA on each genotype
</li>
<li>
Calculates each gene module’s correlation with age
</li>
<li>
Calculates the gene intersection sizes between the modules of each
genotype
</li>
</ul>
</p>

``` r
modules.correlations.func <- function (tissue_name) {
  
  ############################################################################
                                # Subset Tissue
  ############################################################################
  
  tissue_dir <- file.path(modules_path, tissue_name)
  dir.create(tissue_dir)
  
  # Strings that specify genotype and tissue
  N2_t_name <- paste('N2',sep= '_', tissue_name)
  LIPL4_t_name <- paste('LIPL4',sep= '_', tissue_name)
  DAF2_t_name <- paste('DAF2',sep= '_', tissue_name)
  RSKS1_t_name <- paste('RSKS1',sep= '_', tissue_name)
  
  # Subset the tissue
  Idents(full_seurat)<-'tissue'
  T_sub <- subset(full_seurat,idents = tissue_name)
  
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
  age.corr <- function (seurat, ident = "geno", geno, traits,
                        tissue_dir = tissue_dir, tissue_name = tissue_name) {
    
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
    
    saveRDS(corr_pval, file.path(tissue_dir, 
                                 paste0('CP', sep= '_', 
                                        geno, sep= '_',
                                        tissue_name, sep= '',
                                        '.RData')))    
    
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
                    # Test Gene Overlap Between Modules
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

## Below is an example of the overlap results between the modules of each genotype for neurons:

    ##                 M1              M2 M1_size M2_size Intersection_size
    ## 1        N2_yellow      LIPL4_blue     339     677               275
    ## 2        N2_yellow     LIPL4_brown     339     149                14
    ## 3        N2_yellow      DAF2_green     339     201                77
    ## 4        N2_yellow       DAF2_blue     339     310               220
    ## 5        N2_yellow      DAF2_brown     339     243                 1
    ## 6        N2_yellow RSKS1_turquoise     339     601               267
    ## 7        N2_yellow     RSKS1_brown     339     330                 0
    ## 8         N2_brown      LIPL4_blue     389     677                32
    ## 9         N2_brown     LIPL4_brown     389     149                57
    ## 10        N2_brown      DAF2_green     389     201                11
    ## 11        N2_brown       DAF2_blue     389     310                 3
    ## 12        N2_brown      DAF2_brown     389     243                 5
    ## 13        N2_brown RSKS1_turquoise     389     601                45
    ## 14        N2_brown     RSKS1_brown     389     330                13
    ## 15    N2_turquoise      LIPL4_blue     872     677                87
    ## 16    N2_turquoise     LIPL4_brown     872     149                 6
    ## 17    N2_turquoise      DAF2_green     872     201                23
    ## 18    N2_turquoise       DAF2_blue     872     310                11
    ## 19    N2_turquoise      DAF2_brown     872     243                52
    ## 20    N2_turquoise RSKS1_turquoise     872     601                53
    ## 21    N2_turquoise     RSKS1_brown     872     330                98
    ## 22        N2_green      LIPL4_blue      94     677                 0
    ## 23        N2_green     LIPL4_brown      94     149                 0
    ## 24        N2_green      DAF2_green      94     201                 0
    ## 25        N2_green       DAF2_blue      94     310                 0
    ## 26        N2_green      DAF2_brown      94     243                13
    ## 27        N2_green RSKS1_turquoise      94     601                 0
    ## 28        N2_green     RSKS1_brown      94     330                 6
    ## 29      LIPL4_blue       N2_yellow     677     339               275
    ## 30      LIPL4_blue        N2_brown     677     389                32
    ## 31      LIPL4_blue    N2_turquoise     677     872                87
    ## 32      LIPL4_blue        N2_green     677      94                 0
    ## 33      LIPL4_blue      DAF2_green     677     201               119
    ## 34      LIPL4_blue       DAF2_blue     677     310               239
    ## 35      LIPL4_blue      DAF2_brown     677     243                10
    ## 36      LIPL4_blue RSKS1_turquoise     677     601               364
    ## 37      LIPL4_blue     RSKS1_brown     677     330                13
    ## 38     LIPL4_brown       N2_yellow     149     339                14
    ## 39     LIPL4_brown        N2_brown     149     389                57
    ## 40     LIPL4_brown    N2_turquoise     149     872                 6
    ## 41     LIPL4_brown        N2_green     149      94                 0
    ## 42     LIPL4_brown      DAF2_green     149     201                 7
    ## 43     LIPL4_brown       DAF2_blue     149     310                12
    ## 44     LIPL4_brown      DAF2_brown     149     243                 2
    ## 45     LIPL4_brown RSKS1_turquoise     149     601                23
    ## 46     LIPL4_brown     RSKS1_brown     149     330                 4
    ## 47      DAF2_green       N2_yellow     201     339                77
    ## 48      DAF2_green        N2_brown     201     389                11
    ## 49      DAF2_green    N2_turquoise     201     872                23
    ## 50      DAF2_green        N2_green     201      94                 0
    ## 51      DAF2_green      LIPL4_blue     201     677               119
    ## 52      DAF2_green     LIPL4_brown     201     149                 7
    ## 53      DAF2_green RSKS1_turquoise     201     601               134
    ## 54      DAF2_green     RSKS1_brown     201     330                 0
    ## 55       DAF2_blue       N2_yellow     310     339               220
    ## 56       DAF2_blue        N2_brown     310     389                 3
    ## 57       DAF2_blue    N2_turquoise     310     872                11
    ## 58       DAF2_blue        N2_green     310      94                 0
    ## 59       DAF2_blue      LIPL4_blue     310     677               239
    ## 60       DAF2_blue     LIPL4_brown     310     149                12
    ## 61       DAF2_blue RSKS1_turquoise     310     601               233
    ## 62       DAF2_blue     RSKS1_brown     310     330                 0
    ## 63      DAF2_brown       N2_yellow     243     339                 1
    ## 64      DAF2_brown        N2_brown     243     389                 5
    ## 65      DAF2_brown    N2_turquoise     243     872                52
    ## 66      DAF2_brown        N2_green     243      94                13
    ## 67      DAF2_brown      LIPL4_blue     243     677                10
    ## 68      DAF2_brown     LIPL4_brown     243     149                 2
    ## 69      DAF2_brown RSKS1_turquoise     243     601                 4
    ## 70      DAF2_brown     RSKS1_brown     243     330                 5
    ## 71 RSKS1_turquoise       N2_yellow     601     339               267
    ## 72 RSKS1_turquoise        N2_brown     601     389                45
    ## 73 RSKS1_turquoise    N2_turquoise     601     872                53
    ## 74 RSKS1_turquoise        N2_green     601      94                 0
    ## 75 RSKS1_turquoise      LIPL4_blue     601     677               364
    ## 76 RSKS1_turquoise     LIPL4_brown     601     149                23
    ## 77 RSKS1_turquoise      DAF2_green     601     201               134
    ## 78 RSKS1_turquoise       DAF2_blue     601     310               233
    ## 79 RSKS1_turquoise      DAF2_brown     601     243                 4
    ## 80     RSKS1_brown       N2_yellow     330     339                 0
    ## 81     RSKS1_brown        N2_brown     330     389                13
    ## 82     RSKS1_brown    N2_turquoise     330     872                98
    ## 83     RSKS1_brown        N2_green     330      94                 6
    ## 84     RSKS1_brown      LIPL4_blue     330     677                13
    ## 85     RSKS1_brown     LIPL4_brown     330     149                 4
    ## 86     RSKS1_brown      DAF2_green     330     201                 0
    ## 87     RSKS1_brown       DAF2_blue     330     310                 0
    ## 88     RSKS1_brown      DAF2_brown     330     243                 5
    ##             pval
    ## 1  9.803770e-179
    ## 2   2.849248e-01
    ## 3   5.094222e-36
    ## 4  1.780775e-203
    ## 5   1.000000e+00
    ## 6  1.227035e-183
    ## 7   1.000000e+00
    ## 8   9.999985e-01
    ## 9   2.960469e-23
    ## 10  9.797214e-01
    ## 11  1.000000e+00
    ## 12  9.999990e-01
    ## 13  9.383061e-01
    ## 14  9.999354e-01
    ## 15  1.000000e+00
    ## 16  1.000000e+00
    ## 17  9.997597e-01
    ## 18  1.000000e+00
    ## 19  3.550733e-01
    ## 20  1.000000e+00
    ## 21  1.490869e-05
    ## 22  1.000000e+00
    ## 23  1.000000e+00
    ## 24  1.000000e+00
    ## 25  1.000000e+00
    ## 26  2.167173e-03
    ## 27  1.000000e+00
    ## 28  7.398978e-01
    ## 29 9.803770e-179
    ## 30  9.999985e-01
    ## 31  1.000000e+00
    ## 32  1.000000e+00
    ## 33  1.213525e-47
    ## 34 7.626527e-143
    ## 35  1.000000e+00
    ## 36 2.959245e-174
    ## 37  1.000000e+00
    ## 38  2.849248e-01
    ## 39  2.960469e-23
    ## 40  1.000000e+00
    ## 41  1.000000e+00
    ## 42  5.517590e-01
    ## 43  3.889256e-01
    ## 44  9.984997e-01
    ## 45  3.367095e-01
    ## 46  9.976080e-01
    ## 47  5.094222e-36
    ## 48  9.797214e-01
    ## 49  9.997597e-01
    ## 50  1.000000e+00
    ## 51  1.213525e-47
    ## 52  5.517590e-01
    ## 53  2.154763e-70
    ## 54  1.000000e+00
    ## 55 1.780775e-203
    ## 56  1.000000e+00
    ## 57  1.000000e+00
    ## 58  1.000000e+00
    ## 59 7.626527e-143
    ## 60  3.889256e-01
    ## 61 1.093832e-148
    ## 62  1.000000e+00
    ## 63  1.000000e+00
    ## 64  9.999990e-01
    ## 65  3.550733e-01
    ## 66  2.167173e-03
    ## 67  1.000000e+00
    ## 68  9.984997e-01
    ## 69  1.000000e+00
    ## 70  9.999787e-01
    ## 71 1.227035e-183
    ## 72  9.383061e-01
    ## 73  1.000000e+00
    ## 74  1.000000e+00
    ## 75 2.959245e-174
    ## 76  3.367095e-01
    ## 77  2.154763e-70
    ## 78 1.093832e-148
    ## 79  1.000000e+00
    ## 80  1.000000e+00
    ## 81  9.999354e-01
    ## 82  1.490869e-05
    ## 83  7.398978e-01
    ## 84  1.000000e+00
    ## 85  9.976080e-01
    ## 86  1.000000e+00
    ## 87  1.000000e+00
    ## 88  9.999787e-01

# Circos Plots

``` r
# circos.plot <- function (tissue) {
#   
#   tissue_dir <- file.path(path_to_RData, tissue)
#   
#   overlap_file_name <- list.files(tissue_dir, pattern="overlap")
#   overlap <- read.csv(file.path(tissue_dir, overlap_file_name))
#   overlap <- overlap[,-1]
#   
#   CP_N2_file_name <- list.files(tissue_dir, pattern="CP_N2")
#   CP_LIPL4_file_name <- list.files(tissue_dir, pattern="CP_LIPL4")
#   CP_DAF2_file_name <- list.files(tissue_dir, pattern="CP_DAF2")
#   CP_RSKS1_file_name <- list.files(tissue_dir, pattern="CP_RSKS1")
#   
#   CP_N2 <- readRDS(file= file.path(tissue_dir, CP_N2_file_name))
#   CP_LIPL4 <- readRDS(file= file.path(tissue_dir, CP_LIPL4_file_name))
#   CP_DAF2 <- readRDS(file= file.path(tissue_dir, CP_DAF2_file_name))
#   CP_RSKS1 <- readRDS(file= file.path(tissue_dir, CP_RSKS1_file_name))
#   
#   CP_list <- list(CP_N2 , CP_LIPL4, CP_DAF2, CP_RSKS1)
#   names(CP_list) <- c("N2","LIPL4","DAF2","RSKS1")
#   
#   ##############################################################################
#   
#   # Keep overlaps with pvalue < 0.05 
#   filtered <- overlap[overlap$pval<0.05,]
#   
#   # Remove duplicate combinations (I did this based on pvalues. Duplicate combinations 
#   # should have pvalues that are identical up to 3 significant figures)
#   filtered <- filtered[!duplicated(signif(filtered$pval,3)), ] 
#   
#   connections <- filtered[,c(1,2,6)]
#   
#   ##############################################################################
#   
#   # Extract the sizes of each module. 
#   geno_sizes <- overlap[,c(1,3)]
#   geno_sizes2 <- overlap[,c(2,4)]
#   
#   # Duplicates are removed
#   geno_sizes <- geno_sizes[!duplicated(geno_sizes),]
#   geno_sizes2 <- geno_sizes2[!duplicated(geno_sizes2),]
#   names(geno_sizes) <- c("module", "size")
#   names(geno_sizes2) <- c("module", "size")
#   geno_sizes <- rbind(geno_sizes,geno_sizes2)
#   geno_sizes <- geno_sizes[!duplicated(geno_sizes),]
#   
#   # Separate the geno name and module names
#   geno_module_sizes <- separate(data= geno_sizes, col= module, into= c("genotype", "module"), sep = "\\_")
#   
#   # Calculate the total number of genes for each genotype
#   ngenes_genotype <- aggregate(geno_module_sizes$size, 
#                                by=list(Category=geno_module_sizes$genotype), 
#                                FUN=sum)
#   names(ngenes_genotype) <- c("genotype", "size")
#   
#   # Calculate the genotype proportions, this will be used to set the sector sizes.
#   ngenes_genotype$proportion <- proportions(ngenes_genotype$size)
#   
#   ##############################################################################
#   
#   # Extract total amount of genes per genotype and store it in the geno_module_sizes df
#   geno_module_sizes$total_geno_size <- NA
#   genoz <- ngenes_genotype$genotype
#   
#   for (i in 1:length(genoz)) {
#     geno_module_sizes$total_geno_size[geno_module_sizes$genotype==genoz[i]] <- ngenes_genotype$size[ngenes_genotype$genotype==genoz[i]]
#   }
#   
#   # Calculate proportions of each module by genotype
#   geno_module_sizes$proportion <- NA
#   geno_module_sizes$proportion <- geno_module_sizes$size/geno_module_sizes$total_geno_size
#   
#   # Convert geno_module_sizes to a list 
#   split_tibble <- function(tibble, col = 'col') tibble %>% split(., .[, col])
#   mod_list <- split_tibble(geno_module_sizes, 'genotype')
#   
#   # Sort proportions from greatest to least
#   for (i in 1:length(genoz)) {
#     mod_list[[genoz[i]]] <- mod_list[[genoz[i]]][order(mod_list[[genoz[i]]]$proportion, decreasing = TRUE), ]
#   }
#   
#   # Each module will be represented by a box. To create these boxes, line segments are placed on each track where the right side of the box will be. To calculate the end positions, the ordered proportions are added sequentially. The end position of the first module will just be the value of its proportion. The end position of the second module will be the sum of its proportion with the end position of the first module. The end position of the third module will be the sum of its proportion with the end position of the second module. This continues until the end position of the last module is equal to 1.
#   for (i in 1:length(genoz)) {
#     mod_list[[genoz[i]]]$end_position[1] <- mod_list[[genoz[i]]]$proportion[1]
#     n_modules <- nrow(mod_list[[genoz[i]]])
#     x <- 2
#     while (x <= n_modules) {
#       mod_list[[genoz[i]]]$end_position[x] <- mod_list[[genoz[i]]]$proportion[x] + mod_list[[genoz[i]]]$end_position[x-1]
#       x = x + 1
#     }
#   }
#   
#   # Combine geno name with module name
#   for (i in 1:length(genoz)) {
#     mod_list[[genoz[i]]]$geno_mod <- paste(mod_list[[genoz[i]]]$genotype, 
#                                            mod_list[[genoz[i]]]$module, 
#                                            sep= "_")
#   }
#   
#   # Extract correlation values from CP_list and store in mod_list
#   for (i in 1:length(genoz)) {
#     for (j in 1:nrow(mod_list[[genoz[i]]])) {
#       mod <- mod_list[[genoz[i]]]$geno_mod[j]
#       mod_list[[genoz[i]]]$cor[j] <- CP_list[[genoz[i]]][mod,"cor"]
#     }
#   }
# 
#   ##############################################################################
#                                   # Circos Plot
#   ##############################################################################
#   
#   circos_file_name <- paste(tissue, sep= "_", 'circos_plot_reversed.pdf')
#   
#   pdf(file.path(plots_path, circos_file_name), width = 8, height = 8)
#   
#   geno_color <- c(2,3,4,6)
#   circos.par("track.height" = 0.05, gap.after = c(5,5,5,5), cell.padding = c(0, 0, 0, 0))
#   circos.initialize(genoz, xlim = c(0, 1),sector.width = ngenes_genotype$proportion)
#   circos.track(ylim = c(0,1))
#   circos.track(ylim = c(0,1),track.height =0.08)
#   
#   title(tissue)
#   
#   # Add line segments that split the modules of each genotype
#   for (i in 1:length(genoz)) {
#     for (j in 1:nrow(mod_list[[genoz[i]]])) {
#       x <- mod_list[[genoz[i]]]$end_position[j]
#       circos.segments(x,0,x,1, 
#                       sector.index = genoz[i], 
#                       track.index = 2)
#     }
#   }
#   
#   # Color the modules based on their correlation with aging
#   for (i in 1:length(genoz)) {
#     
#     for (j in 1:nrow(mod_list[[genoz[i]]])) {
#       xi <- mod_list[[genoz[i]]]$end_position[j-1]
#       x <- mod_list[[genoz[i]]]$end_position[j]
#       cor <- mod_list[[genoz[i]]]$cor[j]
#       
#       if (cor > 0) {color = "green"}
#       if (cor < 0) {color = "blue"}
#       
#       if (j == 1) {
#         circos.rect(0, 0, x, 1, 
#                     col= color, 
#                     sector.index = genoz[i], 
#                     track.index = 2)
#       }
#       else {
#         circos.rect(xi, 0, x, 1, 
#                     col= color, 
#                     sector.index = genoz[i], 
#                     track.index = 2)
#       }
#     }
#   }
#   
#   # Add the links between modules
#   transp <- 0.8
#   
#   skip_these <- data.frame(matrix(nrow=1,ncol=2))
#   names(skip_these) <- c("M1","M2")
#   
#   for (geno in genoz) {
#     
#     modulez <- connections[grep(geno, connections$M1) ,]
#     modulez <- rbind(modulez, connections[grep(geno, connections$M2) ,])
#     modulez[grep(geno, modulez$M2), c("M2","M1")] <- modulez[grep(geno, modulez$M2), c("M1","M2")]
#     
#     
#     t <- table(modulez$M1)
#     dup_1 <- names(which(t[]>1)) 
#   
#     if (!identical(dup_1, character(0))) {
#       for (i in 1:length(dup_1)) {
#         loc <- grep(dup_1[i], modulez$M1)
#         modulez <- separate(data= modulez, col= M2, into= c("genotype", "module"), sep = "\\_")
#         dummy_mod_1 <- data.frame(matrix(nrow= nrow(modulez), ncol= ncol(modulez)))
#         dummy_mod_1[loc,] <- modulez[loc,]
#         t <- table(dummy_mod_1[,2])
#         dup_g <- names(which(t[]>1)) 
#         
#         if (!identical(dup_g, character(0))) {
#           for (j in 1:length(dup_g)) {
#             loc <- which(dummy_mod_1[,2] == dup_g[j])
#             dummy_mod_2 <- c(rep(NA,nrow(modulez)))
#             dummy_mod_2[loc] <- modulez$pval[loc]
#             loc_min <- which.min(dummy_mod_2)
#             loc_remove <- loc[loc != loc_min]
#             modulez <- modulez[-loc_remove,]
#           }
#         }
#         modulez <- modulez %>% unite("M2", genotype:module, sep= "_", remove= TRUE)
#       }
#     }
#     
#     
#     t <- table(modulez$M2)
#     dup_2 <- names(which(t[]>1)) 
#   
#     if (!identical(dup_2, character(0))) {
#       for (i in 1:length(dup_2)) {
#         loc <- grep(dup_2[i], modulez$M2)
#         dummy_mod <- c(rep(NA,nrow(modulez)))
#         dummy_mod[loc] <- modulez$pval[loc]
#         loc_min <- which.min(dummy_mod)
#         loc_remove <- loc[loc != loc_min]
#         modulez <- modulez[-loc_remove,]
#       }}
#     
#     
#     
#     for (i in 1:nrow(modulez)) {
#       M1 <- modulez$M1[i]
#       M2 <- modulez$M2[i]
#       
#       test1 <- which(skip_these$M1 == M1 & skip_these$M2 == M2)
#       test2 <- which(skip_these$M1 == M2 & skip_these$M2 == M1)
#       
#       if (!identical(test1,integer(0)) | !identical(test2,integer(0))) next
#       
#       geno1 <- gsub("_.*","",M1)
#       geno2 <- gsub("_.*","",M2)
#       
#       loc_a <- which(mod_list[[geno1]]$geno_mod == M1)
#       loc_b <- which(mod_list[[geno2]]$geno_mod == M2)
#       
#       if (loc_a == 1) {
#         M1_i <- 0
#         M1_f <- mod_list[[geno1]][mod_list[[geno1]]$geno_mod == M1,"end_position"]}
#       
#       if (loc_a > 1) {
#         M1_i <- mod_list[[geno1]][loc_a-1,"end_position"]
#         M1_f <- mod_list[[geno1]][mod_list[[geno1]]$geno_mod == M1,"end_position"]}
#       
#       if (loc_b == 1) {
#         M2_i <- 0
#         M2_f <- mod_list[[geno2]][mod_list[[geno2]]$geno_mod == M2,"end_position"]}
#       
#       if (loc_b > 1) {
#         M2_i <- mod_list[[geno2]][loc_b-1,"end_position"]
#         M2_f <- mod_list[[geno2]][mod_list[[geno2]]$geno_mod == M2,"end_position"]}
#       
#       M1_cor <- mod_list[[geno1]][mod_list[[geno1]]$geno_mod == M1,"cor"]
#       M2_cor <- mod_list[[geno2]][mod_list[[geno2]]$geno_mod == M2,"cor"]
#       
#       if (M1_cor > 0 & M2_cor > 0) {color = "green"}
#       if (M1_cor < 0 & M2_cor < 0) {color = "blue"}
#       if (M1_cor > 0 & M2_cor < 0) {color = "red"}
#       if (M1_cor < 0 & M2_cor > 0) {color = "red"}
#         
#       circos.link(geno1, c(M1_i, M1_f), 
#                   geno2, c(M2_i, M2_f), 
#                   col=add_transparency(color, transp)) 
#     }
#     
#     skip_these <- rbind(skip_these, modulez[,1:2])
#     skip_these <- na.omit(skip_these)
#   }
#   
#   
#   # Highlight the genotypes
#   for(i in 1:length(genoz)) {
#     highlight.sector(sector.index = genoz[i], 
#                      track.index = 1, 
#                      col = geno_color[i], 
#                      text = genoz[i], 
#                      text.vjust = -2, 
#                      niceFacing = TRUE)
#   }
#   
#   circos.clear()
#   
#   dev.off()
# 
# }
```

``` r
# for (i in 1:length(tissuez)) {
#   circos.plot(tissuez[i])
# }
```

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
    ##  [1] tidyr_1.3.0           circlize_0.4.15       writexl_1.4.2        
    ##  [4] dplyr_1.1.4           WGCNA_1.72-1          fastcluster_1.2.3    
    ##  [7] dynamicTreeCut_1.63-1 GeneOverlap_1.32.0    SeuratObject_4.1.3   
    ## [10] Seurat_4.3.0.1        hdWGCNA_0.2.03       
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
    ##  [46] leiden_0.4.3           shape_1.4.6            future.apply_1.11.0   
    ##  [49] BiocGenerics_0.42.0    abind_1.4-5            scales_1.2.1          
    ##  [52] DBI_1.1.3              spatstat.random_3.1-5  miniUI_0.1.1.1        
    ##  [55] Rcpp_1.0.11            htmlTable_2.4.1        viridisLite_0.4.2     
    ##  [58] xtable_1.8-4           reticulate_1.35.0      foreign_0.8-85        
    ##  [61] bit_4.0.5              preprocessCore_1.58.0  Formula_1.2-5         
    ##  [64] stats4_4.2.1           htmlwidgets_1.6.2      httr_1.4.7            
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
    ## [115] lifecycle_1.0.4        GlobalOptions_0.1.2    spatstat.geom_3.2-4   
    ## [118] lmtest_0.9-40          RcppAnnoy_0.0.21       data.table_1.14.8     
    ## [121] cowplot_1.1.1          bitops_1.0-7           irlba_2.3.5.1         
    ## [124] httpuv_1.6.11          patchwork_1.1.3        R6_2.5.1              
    ## [127] promises_1.2.1         KernSmooth_2.23-20     gridExtra_2.3         
    ## [130] IRanges_2.30.1         parallelly_1.36.0      codetools_0.2-19      
    ## [133] MASS_7.3-60            gtools_3.9.4           sctransform_0.3.5     
    ## [136] GenomeInfoDbData_1.2.8 S4Vectors_0.34.0       parallel_4.2.1        
    ## [139] rpart_4.1.19           grid_4.2.1             rmarkdown_2.24        
    ## [142] Rtsne_0.16             spatstat.explore_3.2-1 base64enc_0.1-3       
    ## [145] Biobase_2.56.0         shiny_1.7.5
