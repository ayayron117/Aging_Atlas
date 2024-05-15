Bootstraping
================
Aaron Mohammed

``` r
# Bootsrapping function that takes a character vector of tissue names as input

bootstrap.cells <- function (tissue) {

  # Create a folder for each tissue
  dir.create(file.path(bootstrap_path,tissue))
  tissue_path <- file.path(bootstrap_path, tissue)
  
  # Subset the tissue 
  Idents(seurat_data)<-"tissue"
  subset(seurat_data,idents = tissue)->tissue_subset
  Idents(tissue_subset)<-"orig.ident"
  
  ##############################################################################
  ################################ Wild Type ###################################
  ##############################################################################
  
  #### N2 Day 1 ####
  
  # Subset day 1 data for N2
  tmp <- subset(tissue_subset, idents =c("N2D1","N2D1R","N2OPD1"))
  
  # Extract count data from the N2 day 1 subset
  counts <- tmp@assays$RNA@counts 
  metadata <- tmp@meta.data
  sce <- SingleCellExperiment(assays = list(counts = counts), 
                              colData = metadata)
  cells <- colData(sce)[, c("genotype")]
  
  # Only evaluate if there are at least 50 cells
  if (length(cells) >= 50) { 
    
    # Create a character vector that will be used for labeling cells
    # This vector has the same length as the number of cells
    labels <- rep('neg',length(cells))
    
    # Randomly change 15 labels in the vector from neg to pos
    'pos' -> labels[sample(c(1:length(cells)),size=15, replace=FALSE)]
    
    # Aggregate the count matrix according to the labels
    sam <- Matrix.utils:::aggregate.Matrix(t(counts(sce)), 
                           groupings = labels, fun = "sum") 
    
    # Extract the aggregated counts for the cells that were labeled as pos and
    # store in a data frame called Day
    Day <- data.frame(sam=sam['pos',])
    
    for (i in 1:100){ 
      
      # Repeat the same steps 100 more times and each time append the result to
      # the Day data frame
      labels <- rep('neg',length(cells))
      'pos' -> labels[sample(c(1:length(cells)),size=15, replace=FALSE)]
      sam <- Matrix.utils:::aggregate.Matrix(t(counts(sce)), 
                             groupings = labels, fun = "sum") 
      Day <- cbind(Day,sam['pos',])
      
    }
    
    # Save the Day data frame as a csv file in the tissue directory
    write.csv(Day,file.path(tissue_path,"N2_Day_1.csv"))
    
  }
  
  #### N2 Day 6 ###
  
  # Subset day 6 data for N2
  tmp <- subset(tissue_subset,idents =c("N2D6"))
  
  # Extract count data from the N2 day 6 subset
  counts <- tmp@assays$RNA@counts 
  metadata <- tmp@meta.data
  sce <- SingleCellExperiment(assays = list(counts = counts), 
                              colData = metadata)
  cells <- colData(sce)[, c("genotype")]
  
  # Only evaluate if there are at least 50 cells
  if (length(cells) >= 50) { 
    
    # Create a character vector that will be used for labeling cells
    # This vector has the same length as the number of cells
    labels <- rep('neg',length(cells))
    
    # Randomly change 15 labels in the vector from neg to pos
    'pos' -> labels[sample(c(1:length(cells)),size=15, replace=FALSE)]
    
    # Aggregate the count matrix according to the labels
    sam <- Matrix.utils:::aggregate.Matrix(t(counts(sce)), 
                           groupings = labels, fun = "sum") 
    
    # Extract the aggregated counts for the cells that were labeled as pos and
    # store in a data frame called Day
    Day <- data.frame(sam=sam['pos',])
    
    for (i in 1:100){ 
      
      # Repeat the same steps 100 more times and each time append the result to
      # the Day data frame
      labels <- rep('neg',length(cells))
      'pos' -> labels[sample(c(1:length(cells)),size=15, replace=FALSE)]
      sam <- Matrix.utils:::aggregate.Matrix(t(counts(sce)), 
                             groupings = labels, fun = "sum") 
      Day <- cbind(Day,sam['pos',])
      
    }
    
    # Save the Day data frame as a csv file in the tissue directory
    write.csv(Day,file.path(tissue_path,"N2_Day_6.csv"))
    
  }
  
  ### N2 Day 12 ###
  
  # Subset day 12 data for N2
  tmp <- subset(tissue_subset,idents =c("N2D12"))
  
  # Extract count data from the N2 day 12 subset
  counts <- tmp@assays$RNA@counts 
  metadata <- tmp@meta.data
  sce <- SingleCellExperiment(assays = list(counts = counts), 
                              colData = metadata)
  cells <- colData(sce)[, c("genotype")]
  
  # Only evaluate if there are at least 50 cells
  if (length(cells) >= 50) { 
    
    # Create a character vector that will be used for labeling cells
    # This vector has the same length as the number of cells
    labels <- rep('neg',length(cells))
    
    # Randomly change 15 labels in the vector from neg to pos
    'pos' -> labels[sample(c(1:length(cells)),size=15, replace=FALSE)]
    
    # Aggregate the count matrix according to the labels
    sam <- Matrix.utils:::aggregate.Matrix(t(counts(sce)), 
                           groupings = labels, fun = "sum") 
    
    # Extract the aggregated counts for the cells that were labeled as pos and
    # store in a data frame called Day
    Day <- data.frame(sam=sam['pos',])
    
    for (i in 1:100){ 
      
      # Repeat the same steps 100 more times and each time append the result to
      # the Day data frame
      labels <- rep('neg',length(cells))
      'pos' -> labels[sample(c(1:length(cells)),size=15, replace=FALSE)]
      sam <- Matrix.utils:::aggregate.Matrix(t(counts(sce)), 
                             groupings = labels, fun = "sum") 
      Day <- cbind(Day,sam['pos',])
      
    }
    
    # Save the Day data frame as a csv file in the tissue directory
    write.csv(Day,file.path(tissue_path,"N2_Day_12.csv"))
    
  }
  
  ### N2 Day 14 ###
  
  # Subset day 14 data for N2
  tmp <- subset(tissue_subset,idents =c("N2D14"))
  
  # Extract count data from the N2 day 14 subset
  counts <- tmp@assays$RNA@counts 
  metadata <- tmp@meta.data
  sce <- SingleCellExperiment(assays = list(counts = counts), 
                              colData = metadata)
  cells <- colData(sce)[, c("genotype")]
  
  # Only evaluate if there are at least 50 cells
  if (length(cells) >= 50) { 
    
    # Create a character vector that will be used for labeling cells
    # This vector has the same length as the number of cells
    labels <- rep('neg',length(cells))
    
    # Randomly change 15 labels in the vector from neg to pos
    'pos' -> labels[sample(c(1:length(cells)),size=15, replace=FALSE)]
    
    # Aggregate the count matrix according to the labels
    sam <- Matrix.utils:::aggregate.Matrix(t(counts(sce)), 
                           groupings = labels, fun = "sum") 
    
    # Extract the aggregated counts for the cells that were labeled as pos and
    # store in a data frame called Day
    Day <- data.frame(sam=sam['pos',])
    
    for (i in 1:100){ 
      
      # Repeat the same steps 100 more times and each time append the result to
      # the Day data frame
      labels <- rep('neg',length(cells))
      'pos' -> labels[sample(c(1:length(cells)),size=15, replace=FALSE)]
      sam <- Matrix.utils:::aggregate.Matrix(t(counts(sce)), 
                             groupings = labels, fun = "sum") 
      Day <- cbind(Day,sam['pos',])
      
    }
    
    # Save the Day data frame as a csv file in the tissue directory
    write.csv(Day,file.path(tissue_path,"N2_Day_14.csv"))
    
  }
  
  
  ##############################################################################
  ########################### Longevity Genotypes ##############################
  ##############################################################################
  
  ################################# lipl-4 Tg ##################################
  
  ### LIPL4 Day 6 ###
  
  # Subset day 6 data for LIPL4
  tmp <- subset(tissue_subset,idents =c("LIPL4D6"))
  
  # Extract count data from the LIPL4 day 6 subset
  counts <- tmp@assays$RNA@counts 
  metadata <- tmp@meta.data
  sce <- SingleCellExperiment(assays = list(counts = counts), 
                              colData = metadata)
  cells <- colData(sce)[, c("genotype")]
  
  # Only evaluate if there are at least 50 cells
  if (length(cells) >= 50) { 
    
    # Create a character vector that will be used for labeling cells
    # This vector has the same length as the number of cells
    labels <- rep('neg',length(cells))
    
    # Randomly change 15 labels in the vector from neg to pos
    'pos' -> labels[sample(c(1:length(cells)),size=15, replace=FALSE)]
    
    # Aggregate the count matrix according to the labels
    sam <- Matrix.utils:::aggregate.Matrix(t(counts(sce)), 
                           groupings = labels, fun = "sum") 
    
    # Extract the aggregated counts for the cells that were labeled as pos and
    # store in a data frame called Day
    Day <- data.frame(sam=sam['pos',])
    
    for (i in 1:100){ 
      
      # Repeat the same steps 100 more times and each time append the result to
      # the Day data frame
      labels <- rep('neg',length(cells))
      'pos' -> labels[sample(c(1:length(cells)),size=15, replace=FALSE)]
      sam <- Matrix.utils:::aggregate.Matrix(t(counts(sce)), 
                             groupings = labels, fun = "sum") 
      Day <- cbind(Day,sam['pos',])
      
    }
    
    # Save the Day data frame as a csv file in the tissue directory
    write.csv(Day,file.path(tissue_path,"LIPL4_Day_6.csv"))
    
  }
  
  ################################# daf-2(lf) ##################################
  
  ### DAF2 Day 6 ###
  
  # Subset day 6 data for DAF2
  tmp <- subset(tissue_subset,idents =c("DAF2D6"))
  
  # Extract count data from the DAF2 day 6 subset
  counts <- tmp@assays$RNA@counts 
  metadata <- tmp@meta.data
  sce <- SingleCellExperiment(assays = list(counts = counts), 
                              colData = metadata)
  cells <- colData(sce)[, c("genotype")]
  
  # Only evaluate if there are at least 50 cells
  if (length(cells) >= 50) { 
    
    # Create a character vector that will be used for labeling cells
    # This vector has the same length as the number of cells
    labels <- rep('neg',length(cells))
    
    # Randomly change 15 labels in the vector from neg to pos
    'pos' -> labels[sample(c(1:length(cells)),size=15, replace=FALSE)]
    
    # Aggregate the count matrix according to the labels
    sam <- Matrix.utils:::aggregate.Matrix(t(counts(sce)), 
                           groupings = labels, fun = "sum") 
    
    # Extract the aggregated counts for the cells that were labeled as pos and
    # store in a data frame called Day
    Day <- data.frame(sam=sam['pos',])
    
    for (i in 1:100){ 
      
      # Repeat the same steps 100 more times and each time append the result to
      # the Day data frame
      labels <- rep('neg',length(cells))
      'pos' -> labels[sample(c(1:length(cells)),size=15, replace=FALSE)]
      sam <- Matrix.utils:::aggregate.Matrix(t(counts(sce)), 
                             groupings = labels, fun = "sum") 
      Day <- cbind(Day,sam['pos',])
      
    }
    
    # Save the Day data frame as a csv file in the tissue directory
    write.csv(Day,file.path(tissue_path,"DAF2_Day_6.csv"))
    
  }
  
  ################################# rsks-1(lf) #################################
  
  ### RSKS1 Day 6 ###
  
  # Subset day 6 data for RSKS1
  tmp <- subset(tissue_subset,idents =c("RSKS1D6"))
  
  # Extract count data from the RSKS1 day 6 subset
  counts <- tmp@assays$RNA@counts 
  metadata <- tmp@meta.data
  sce <- SingleCellExperiment(assays = list(counts = counts), 
                              colData = metadata)
  cells <- colData(sce)[, c("genotype")]
  
  # Only evaluate if there are at least 50 cells
  if (length(cells) >= 50) { 
    
    # Create a character vector that will be used for labeling cells
    # This vector has the same length as the number of cells
    labels <- rep('neg',length(cells))
    
    # Randomly change 15 labels in the vector from neg to pos
    'pos' -> labels[sample(c(1:length(cells)),size=15, replace=FALSE)]
    
    # Aggregate the count matrix according to the labels
    sam <- Matrix.utils:::aggregate.Matrix(t(counts(sce)), 
                           groupings = labels, fun = "sum") 
    
    # Extract the aggregated counts for the cells that were labeled as pos and
    # store in a data frame called Day
    Day <- data.frame(sam=sam['pos',])
    
    for (i in 1:100){ 
      
      # Repeat the same steps 100 more times and each time append the result to
      # the Day data frame
      labels <- rep('neg',length(cells))
      'pos' -> labels[sample(c(1:length(cells)),size=15, replace=FALSE)]
      sam <- Matrix.utils:::aggregate.Matrix(t(counts(sce)), 
                             groupings = labels, fun = "sum") 
      Day <- cbind(Day,sam['pos',])
      
    }
    
    # Save the Day data frame as a csv file in the tissue directory
    write.csv(Day,file.path(tissue_path,"RSKS1_Day_6.csv"))
    
  }

}
```

``` r
# Run the bootstrap function
for(i in 1:length(tissuez)) {
bootstrap.cells(tissuez[i])
}
```

``` r
# Here is an example of what one of the data frames looks like. What's displayed
# is the first 10 genes and first 5 columns of aggregated counts. Each data 
# frame has 26208 rows and 101 columns

tissue_path <- file.path(bootstrap_path, "Neuron")

N2_day_1 <- read.csv(file.path(tissue_path, "N2_Day_1.csv"), row.names = 1)
N2_day_1[1:10, 1:5]
```

    ##        sam sam..pos.... sam..pos.....1 sam..pos.....2 sam..pos.....3
    ## nduo-6  13           13              3              7             13
    ## ndfl-4   3            3              0              0              7
    ## MTCE.7   3            5              1              3             11
    ## nduo-1   6            4              0              0             10
    ## atp-6   13           15              8             13             32
    ## nduo-2   7            1              3              0              4
    ## ctb-1    1            9              3              4             20
    ## ctc-3   15           12              4              8             54
    ## nduo-4   4            4              1              5              7
    ## ctc-1   26           16             17             19             68

``` r
# Not all tissues contained all time points for N2. The ones that did have
# day 1, day 6, day 12, day 14 data were copied to a new directory using a 
# function I wrote called filter.by.N2

# Create a directory where the files for the tissues with all 4 time points for 
# N2 will be copied to 
all_timepoints_dir <- file.path(bootstrap_path,c("tissues_with_all_N2_timepoints"))
dir.create(all_timepoints_dir)

# Filtering function
filter.by.N2 <- function (tissue) {
  
  # List the files in that contain 'N2' in their name within the tissue 
  # directory being evaluated
  contents <- dir(file.path(bootstrap_path, tissue),pattern = "N2")
  
  # Determine the number of files that contain 'N2' in their name and store that
  # number in a variable called l 
  l <- length(contents)

  # If l is equal to 4, copy the tissue directory being evaluated to the
  # 'tissues_with_all_N2_timepoints' directory
  if (l == 4) {
    new_tissue_dir <- file.path(all_timepoints_dir,tissue)
    copyDirectory(file.path(bootstrap_path, tissue), new_tissue_dir)
  }
  
}
  
# Run the filtering function
for (i in 1:length(tissuez)) {
  filter.by.N2(tissuez[i])
}

# These were the tissues with all 4 time points for N2
tissuez <- list.dirs(all_timepoints_dir,full.names=FALSE, recursive= FALSE)
tissuez
```

    ## [1] "Embroynic_cells" "Germline"        "Hypodermis"      "Intestine"      
    ## [5] "Muscle"          "Neuron"          "Pharynx"         "Vulva_uterus"
