Bootstraping
================
Aaron Mohammed

``` r
library(Seurat)
library(SingleCellExperiment)
library(R.utils)

readRDS(raw_data_path) -> seurat_data

Idents(seurat_data)<-"tissue"
tissuez <- names(table(seurat_data[[]]$tissue))

genotypez <- list(c("N2D1","N2D1R","N2OPD1"),
                  "N2D6","N2D12","N2D14",
                  "LIPL4D6","DAF2D6","RSKS1D6")
```

``` r
for (t in tissuez) {
  
  # Create a folder for each tissue
  dir.create(file.path(bootstrap_path,t))
  tissue_path <- file.path(bootstrap_path, t)
  
  # Subset the tissue 
  Idents(seurat_data) <- "tissue"
  subset(seurat_data,idents = t) -> tissue_subset
  Idents(tissue_subset) <- "orig.ident"
  
  for (g in genotypez) {
    # Subset genotype
    tmp <- subset(tissue_subset, idents = g)
    
    # Extract count data from genotype subset
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
      'pos' -> labels[sample(c(1:length(cells)),
                             size=15, 
                             replace=FALSE)]
      
      # Aggregate the count matrix according to the labels
      labeled <- Matrix.utils:::aggregate.Matrix(t(counts(sce)), 
                                                 groupings = labels,
                                                 fun = "sum") 
      
      # Extract the aggregated counts for the cells that were labeled as pos and
      # store in a data frame called Day
      Day <- data.frame(labeled=labeled['pos',])
      
      # Repeat the same steps 100 more times and each time append the result to
      # the Day data frame
      for (i in 1:100){ 
        
        labels <- rep('neg',length(cells))
        'pos' -> labels[sample(c(1:length(cells)),
                               size=15,
                               replace=FALSE)]
        labeled <- Matrix.utils:::aggregate.Matrix(t(counts(sce)), 
                                                   groupings = labels, 
                                                   fun = "sum") 
        Day <- cbind(Day,labeled['pos',])
      }
      
        # Save the Day data frame as a csv file in the tissue directory
        if (length(g) == 3) {
          write.csv(Day, file.path(tissue_path, "N2D1.csv"))
        } else {
          write.csv(Day, file.path(tissue_path, 
                                   paste0(g, ".csv")))
        }
    }
  }
}
```

``` r
# Here is an example of what one of the data frames looks like. What's displayed
# is the first 10 genes and first 5 columns of aggregated counts. Each data 
# frame has 26208 rows and 101 columns

tissue_path <- file.path(bootstrap_path, "Neuron")

N2D1 <- read.csv(file.path(tissue_path, "N2D1.csv"), row.names = 1)
N2D1[1:10, 1:5]
```

    ##        labeled labeled..pos.... labeled..pos.....1 labeled..pos.....2
    ## nduo-6      13               13                  3                  7
    ## ndfl-4       3                3                  0                  0
    ## MTCE.7       3                5                  1                  3
    ## nduo-1       6                4                  0                  0
    ## atp-6       13               15                  8                 13
    ## nduo-2       7                1                  3                  0
    ## ctb-1        1                9                  3                  4
    ## ctc-3       15               12                  4                  8
    ## nduo-4       4                4                  1                  5
    ## ctc-1       26               16                 17                 19
    ##        labeled..pos.....3
    ## nduo-6                 13
    ## ndfl-4                  7
    ## MTCE.7                 11
    ## nduo-1                 10
    ## atp-6                  32
    ## nduo-2                  4
    ## ctb-1                  20
    ## ctc-3                  54
    ## nduo-4                  7
    ## ctc-1                  68

``` r
# Not all tissues contained all time points for N2. The ones that did have
# day 1, day 6, day 12, day 14 data were copied to a new directory using a 
# function I wrote called filter.by.N2

# Create a directory where the files for the tissues with all 4 time points for 
# N2 will be copied to 
all_timepoints_dir <- file.path(bootstrap_path,c("tissues_w_all_N2_timepoints"))
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
    new_tissue_dir <- file.path(all_timepoints_dir, 
                                tissue)
    
    copyDirectory(file.path(bootstrap_path, tissue), 
                  new_tissue_dir)
  }
}
```

``` r
# Run the filtering function
for (i in 1:length(tissuez)) {
  filter.by.N2(tissuez[i])
}

# These were the tissues with all 4 time points for N2
tissuez <- list.dirs(all_timepoints_dir,
                     full.names=FALSE, 
                     recursive= FALSE)

tissuez
```

    ## [1] "Embroynic_cells" "Germline"        "Hypodermis"      "Intestine"      
    ## [5] "Muscle"          "Neuron"          "Pharynx"         "Vulva_uterus"
