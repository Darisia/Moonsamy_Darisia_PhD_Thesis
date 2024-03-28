# This R script contains the code to perform all downstream processing steps
# .. on the scRNA-seq datasets used in chapter three

# First we need to load all required packages
library(tidyverse)
library(Seurat)

####### Functions created in R to facilitate this process

# General R function to load an RData file using a specified object name
loadRData <- function(filename) {
  load(filename)
  get(ls()[ls() != "filename"])
}

# R function to replace the ENSEMBL IDs with hgnc gene names within a Seurat object
ensg_to_hgnc <- function(
  ENSG_gname38_path = "working_directory",
  seurat_object,
  sample_ID
){
  # Set the path to the ENSg_gname38.tsv file
  if (ENSG_gname38_path == "working_directory") {
    ENSG_gname38_path = "./ENSG_gname38.tsv"
  } else {
    ENSG_gname38_path = ENSG_gname38_path
  }
  # read in the ENSG_gname38.tsv file
  ENSG_gname38 <- read_tsv(ENSG_gname38_path)
  # Read in the saved seurat object
  if (endsWith(seurat_object, ".rds") == TRUE | endsWith(seurat_object, "RDS") == TRUE) {
    sample_ID_srt <- readRDS(file = seurat_object)
  } else {
    sample_ID_srt <- loadRData(seurat_object)
  }
  # Update Seurat object if saved as an old seurat object
  sample_ID_srt <- UpdateSeuratObject(sample_ID_srt)
  # Generate the dataframe that will be used to match and replace the ENSEMBL IDs
  pmatch_output <- data.frame(pmatch(ENSG_gname38$ENSEMBL_ID, rownames(sample_ID_srt@assays[["RNA"]]@counts)))
  colnames(pmatch_output) <- "gene_match"
  pmatch_output$match_index <- rownames(pmatch_output)
  pmatch_output <- pmatch_output[order(pmatch_output$gene_match), ]
  pmatch_output <- filter(pmatch_output, !is.na(gene_match))
  # Match and replace the ENSEMBL IDs in the relevant slots in the seurat object with the hgnc gene symbol
  rownames(sample_ID_srt@assays[["RNA"]]@counts) <- ENSG_gname38$gene_name[as.double(pmatch_output$match_index)]
  # Return the Seurat object with hgnc symbols as rownames
  return(sample_ID_srt)
}

# Read in and check the seurat object
read_check_srt <- function(
  seurat_object,
  sample_ID,
  ENSG_gname38_path = "working_directory"
){
  # Read in the saved seurat object
  if (endsWith(seurat_object, ".rds") == TRUE | endsWith(seurat_object, "RDS") == TRUE) {
    sample_ID_srt <- readRDS(file = seurat_object)
  } else {
    sample_ID_srt <- loadRData(seurat_object)
  }
  # Update Seurat object if saved as an old seurat object
  sample_ID_srt <- UpdateSeuratObject(sample_ID_srt)
  # Check if rownames use hgnc symbols or ENSEMBL IDs - change to hgnc symbol if ENSEMBL IDs are used
  if (startsWith(rownames(sample_ID_srt@assays[["RNA"]]@counts)[1], "ENSG") == TRUE) {
    sample_ID_srt <- ensg_to_hgnc(ENSG_gname38_path, seurat_object, sample_ID)
  }
  # Remove cell barcodes with counts of 0
  sample_ID_srt <- subset(sample_ID_srt, subset = nCount_RNA > 0)
  # Check if the mito.percent column is in the meta.data slot of the Seurat object
  if ("mito.percent" %in% colnames(sample_ID_srt@meta.data)) {
    sample_ID_srt <- sample_ID_srt
  } else {
    sample_ID_srt[['mito.percent']] <- PercentageFeatureSet(sample_ID_srt, pattern = '^MT-')
  }
  return(sample_ID_srt)
}

############################# Function for editing plots ##########
# Function for editing plots
edit_plots <- function(
  ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 14),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 14),
          axis.text.x = element_text(size = 12,
                                     colour = 'black'),
          axis.text.y = element_text(size = 12,
                                     colour = 'black'),
          legend.text = element_text(size = 12,
                                     colour = 'black'),
          legend.title = element_text(size = 14,
                                      colour = 'black'),
          strip.text = element_text(
            size = 14,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

########################
# small function to perform cell-cycle regression
########################

basic_ccreg <- function(
    seurat_object_path,
    sample_ID,
    output_directory,
    set_ident = TRUE
) {
    # Read in the Seurat object --- seurat object should contain the hgnc gene symbols at this point
    sample_ID_srt <- read_check_srt(
        seurat_object_path,
        smaple_ID
    )

    # Set up the cell-cycle information
    s.genes <- cc.genes.updated.2019$s.genes
    g2m.genes <- cc.genes.updated.2019$g2m.genes
    sample_ID_srt <- CellCycleScoring(
        sample_ID_srt,
        s.features = s.genes,
        g2m.features = g2m.genes,
        set.ident = set.ident
    )

    # Set up the vasr.to.regress variable
    vars.to.regress <- c("S.Score", "G2M.Score")

    # Run the normalisation function to regress out cell-cycle variation
    sample_ID_srt <- SCTransform(
        sample_ID_srt,
        vars.to.regress = vars.to.regress
    )

    # Save the cell-cycle regressed object to file
    save_to_dir <- sprintf(
      '%s/%s_ccreg.rds',
      output_directory,
      sample_ID
    )
    saveRDS(
      object = sample_ID_srt,
      file = save_to_dir
    )

    return(sample_ID_srt)

}
########################
# steps for merging the datasets
########################
# First we need the paths to the seurat objects we want to merge
sample_paths <- c()

# We need to read in the seurat objects and store them in a list
samples_list <- list()
for (path in sample_paths) {
    srt_object <- readRDS(path)
    samples_list <- c(
        samples_list,
        list(srt_object)
    )
}

# Now we run the merge function to merge all of the seurat objects
# .. Can rename the output if I need to
merged_samples <- merge(
    x = samples_list[[1]],
    y = c(
        samples_list[[2]],
        samples_list[[3]],
        samples_list[[4]]
    ),
    add.cell.ids = c( # Add the cell names - keep it the same as orig.ident
        "samplename1",
        "samplename2",
        "samplename3",
        "samplename4"
    ),
    merge.data = TRUE, # Think this is actually default
    project = "projectname" # name when performing these steps
)

########################
# small function for integrating datasets
########################

basic_combine <- function(
    list_of_srt_objects, # List of seurat objects
    normalization_method, # String
    sample_ID,
    output_directory,
    reference = NULL
) {
    # select repeatedly variable features
    common_features <- SelectIntegrationFeatures(
        object.list = list_of_srt_objects
    )

    # if normalisation method = SCT, then run PrepSCTIntegration
    if (normalization_method == "SCT") {
      list_of_srt_objects <- PrepSCTIntegration(
        list_of_srt_objects,
        anchor.features = common_features,
        assay = "SCT"
      )
    }

    # Identify anchor cells that will be used for integration
    integration_anchors <- FindIntegrationAnchors(
        object.list = list_of_srt_objects,
        anchor.features = common_features,
        normalization.method = normalization_method,
        reference = reference
    )

    # Combine the datasets using the anchors and produce new combined seurat object
    combined_datasets_srt <- IntegrateData(
        anchorset = integration_anchors,
        normalization.method = normalization_method
    )

    # Save the combined datasets to file as well
    save_to_dir <- sprintf(
      '%s/%s_integrated_datasets.rds',
      output_directory,
      sample_ID
  )
    saveRDS(
      object = combined_datasets_srt,
      file = save_to_dir
    )

    return(combined_datasets_srt)
}

########################
# Small function for reducing the dimensions and visualising PCs
########################

basic_reduction <- function(
    seurat_object,
    sample_ID,
    output_directory
) {
    # Run the PCA on the dataset to find principal components
    seurat_object_reduced <- RunPCA(
        object = seurat_object
    ) # No other parameters needed

    # Extracting standard deviation for 50 principal components to create a pretty elbow plot
    extracted_std_devs <- seurat_object_reduced@reductions$pca@stdev
    extracted_std_dev_df <- data.frame(extracted_std_devs)
    colnames(extracted_std_dev_df) <- "standard_deviation"

    # Creating the actual plot
    std_dev_plot <- ggplot(
        extracted_std_dev_df,
        aes(
            x = as.double(rownames(extracted_std_dev_df)),
            y = standard_deviation
        )
    ) +
    geom_point(
        colour = "#007fff"
    ) +
    labs(
        x = "Principal Components",
        y = "Standard deviation"
    ) +
    theme(
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")
    )

    std_dev_plot <- edit_plots(std_dev_plot)

    # Save the plot to file
    plot_save_to_dir = sprintf(
        "%s/%s_stdev_pc.png",
        output_directory,
        sample_ID
    )
    ggsave(
        file = plot_save_to_dir,
        std_dev_plot,
        height = 8,
        width = 11,
        dpi = 900
    )

    # Save the seurat object with reduced dimensions aswell
    save_to_dir = sprintf(
        "%s/%s_dim_reduced.rds",
        output_directory,
        sample_ID
    )

    saveRDS(
        object = seurat_object_reduced,
        file  = save_to_dir
    )

    # Return the seurat object so other analyses can be easily performed
    return(seurat_object_reduced)
}

########################
# Small function for clustering and visualising the clusters
########################

basic_cluster <- function(
    seurat_object,
    sample_ID,
    output_directory,
    dimensions, # Given in the form of dimensions - eg -> 1:15
    resolution
) {
    # Finding the neighbours
    seurat_object <- FindNeighbors(
        seurat_object,
        dims = dimensions
    )

    # Finding the clusters
    seurat_object <- FindClusters(
        seurat_object,
        resolution = resolution
    )

    # Visualisation of clusters using UMAP
    seurat_object <- RunUMAP(
        seurat_object,
        dims = dimensions
    )
    # Basic UMAP visualisation (everything combined)
    basic_UMAP_visualisation <- DimPlot(
        seurat_object,
        reduction = "umap",
        label = TRUE,
        repel = TRUE
    )
    # detailed UMAP visualisation showing clusters for the different datasets
    detailed_UMAP_visualisation <- DimPlot(
        seurat_object,
        reduction = "umap",
        split.by = "orig.ident",
        ncol = 2,
        label = TRUE,
        repel = TRUE
    )

    # Saving the basic UMAP visualisation
    basic_UMAP_save_to_dir <- sprintf(
        "%s/%s_basic_clusters_UMAP.png",
        output_directory,
        sample_ID
    )

    ggsave(
        file = basic_UMAP_save_to_dir,
        basic_UMAP_visualisation,
        width = 8,
        height = 8,
        dpi = 900
    )

    # Saving the detailed UMAP visualisation
    detailed_UMAP_save_to_dir <- sprintf(
        "%s/%s_detailed_clusters_UMAP.png",
        output_directory,
        sample_ID
    )

    ggsave(
        file = detailed_UMAP_save_to_dir,
        detailed_UMAP_visualisation,
        width = 8,
        height = 8,
        dpi = 900
    )

    # Save the clustered seurat object to file
    save_to_dir <- sprintf(
        "%s/%s_clustered.rds",
        output_directory,
        sample_ID
    )

    saveRDS(
        object = seurat_object,
        file = save_to_dir
    )

    return(seurat_object)
}

########################
# Small function for visualising expression of the known celltype markers
########################

known_marker_visualisations <- function(
    seurat_object, # Seurat object already read into R
    sample_ID,
    split_by_col, # The column to split the visualisations by (Generally the orig.ident column)
    output_directory,
    ncol, # Not sure if I even use this parameter - but should just set to 2 when using the function
    colours # Need to include colours to represent the different groups if more than 1 group is visualised
){
    # Create list of the PBMC cell type markers
    PBMC_celltype_markers <- list(
        "CD14_Monocytes_Classical" = c(
         "CD68",
         "CD14",
         "CCL2"
        ),
        "CD16_Monocytes_Non_classical" = c(
            "CD68",
            "FCGR3A",
            "MS4A7"
        ),
        "B_cells" = c(
            "MS4A1",
            "CD79A"
        ),
        "T_cells" = c(
            "CD3D",
            "CD3G"
        ),
        "Cytotoxic_cells" = c(
            "GNLY",
            "NKG7",
            "GZMA"
        ),
        "NK_Cytotoxic_cells" = c(
            "GNLY",
            "NKG7",
            "KLRF1"
        ),
        "CD8_T_cells" = c(
            "CD8A",
            "CD8B",
            "GZMA",
            "CCL5"
        ),
        "Blood_platelets" = c(
            "PPBP"
        ),
        "All_Monocytes" = c(
            "CD68",
            "LYZ"
        ),
        "Dendritic_cells" = c(
            "LYZ",
            "FCER1A",
            "HLA-DQA1",
            "GPR183"
        ),
        "Naive_T_cells" = c(
            "PTPRC",
            "SELL",
            "IL7R",
            "CCR7"
        ),
        "Central_Memory_T_cells" = c(
            "CCR7",
            "CD44",
            "S100A4"
        ),
        "Effector_Memory_T_cells" = c(
            "PTPRC"
        )
    )
    # Create PBMC celltype folders
    PBMC_celltype_folders <- c(
        "Monocytes",
        "Tcells",
        "Bcells",
        "DCs",
        "CytotoxicCells",
        "Platelets"
    )
    # Extract the cell-type names for each celltype in the list
    known_celltypes <- names(PBMC_celltype_markers)
    # Vector of patterns to subset the PBMC_celltype_markers list
    PBMC_celltype_patterns <- c(
        "Monocytes",
        "T_cells",
        "B_cells",
        "Dendritic_cells",
        "Cytotoxic_cells",
        "platelets"
    )
    # Subset the list according to the patterns
    # Create the empty list to store the subset lists
    PBMC_subset_by_celltype <- list()
    # Create the loop to subset the main celltype list
    for ( pattern in PBMC_celltype_patterns ) {
        # Subset the list according to the pattern
        subset_list <- PBMC_celltype_markers[grep(pattern, names(PBMC_celltype_markers))]
        # Store the subset list in the list create outside the loop
        PBMC_subset_by_celltype <- c(
            PBMC_subset_by_celltype,
            list(subset_list)
        )
    }
    # Set the names of the subset lists within the combined lists
    PBMC_subset_by_celltype <- setNames(
        PBMC_subset_by_celltype,
        PBMC_celltype_folders
    )
    # Create vector of major markers
    PBMC_major_markers <- c(
        "CD3D",
        "CD14",
        "FCGR3A",
        "CD68",
        "GZMA",
        "GNLY",
        "MS4A1",
        "FCER1A",
        "CCL5",
        "CCR7",
        "CD8A",
        "CD8B",
        "PPBP"
    )

    # Set the default assay of the seurat object to "SCT"
    DefaultAssay(seurat_object) <- "RNA"
    # DefaultAssay(seurat_object) <- "SCT"
    # For loop creating featureplot visualisations for all the celltype markers - should create a new folder to store these images
    for ( i in seq(1, length(PBMC_subset_by_celltype)) ) {
        for ( marker_set in names(PBMC_subset_by_celltype[[i]]) ) {
            for ( gene_marker in PBMC_subset_by_celltype[[i]][[marker_set]] ) {
                markers_feature <- FeaturePlot(
                    object = seurat_object,
                    features = gene_marker
                )
                # Create the filename for the featureplot
                markers_feature_name <- sprintf(
                    "%s/%s/%s_%s.png",
                    output_directory,
                    PBMC_celltype_folders[i],
                    marker_set,
                    gene_marker
                )
                # Save the plot to destination specified
                ggsave(
                    filename = markers_feature_name,
                    plot = markers_feature,
                    width = 8,
                    height = 7,
                    dpi = 900
                )
            }
        }
    }

    # Create the dotplot showing all the major cell type markers - splitting by dataset
    markers_dotplot <- DotPlot(
        object = seurat_object,
        features = PBMC_major_markers,
        dot.scale = 8,
        split.by = split_by_col,
        cols = colours
    ) +
    RotatedAxis() +
    theme_bw()
    # Create the filename for saving the dotplot
    markers_dotplot_name <- sprintf(
        "%s/%s_major_markers_dotplot.png",
        output_directory,
        sample_ID
    )
    # Save the dotplot to file
    ggsave(
        filename = markers_dotplot_name,
        plot = markers_dotplot,
        width = 8 ,
        height = 12,
        dpi = 900
    )
    # Create the dotplot showing all the major cell type markers - not splitting by dataset
    markers_combined_dotplot <- DotPlot(
        object = seurat_object,
        features = PBMC_major_markers,
        dot.scale = 8,
        cols = c("lightgrey", "#ff007f")
    ) +
    RotatedAxis() +
    theme_bw()
    # Create the filename for saving the dotplot
    markers_combined_dotplot_name <- sprintf(
        "%s/%s_major_markers_combined_dotplot.png",
        output_directory,
        sample_ID
    )
    # Save the dotplot to file
    ggsave(
        filename = markers_combined_dotplot_name,
        plot = markers_combined_dotplot,
        width = 8,
        height = 9,
        dpi = 900
    )
}

#####################################################
# Running the functions above for the IFNStim datasets
#####################################################

# Reminder: Code names for each of the datasets for record
# IFNAStim - 227291
# IFNLStim - 227292
# Unstim1 - 227293
# Unstim2 - 227294
# IFNBStim - 227295

## 14 September 2022 - Running downstream steps (PRJNA597786 and PRJNA381100) - 2291401 - standard integration workflow
# NB - new data processing code for integrated dataset: 2291401
#Create vector of paths the datasets to be integrated
sample_paths <- c(
        "1.Preprocessing/IFNStimMain/227291/IFNA_227291_logNorm.rds", # Stim with IFN-alpha (healthy)
        "1.Preprocessing/IFNStimMain/227292/IFNL_227292_logNorm.rds", # Stim with IFN-lambda (healthy)
        "1.Preprocessing/IFNStimMain/227293/UnS_227293_logNorm.rds", # Unstim (Healthy)
        "1.Preprocessing/IFNStimMain/227294/UnS_227294_logNorm.rds", # Stim with IFN-beta (SLE)
        "1.Preprocessing/IFNStimMain/227295/IFNB_227295_logNorm.rds" # Unstim (SLE)
)
#  Read in the samples and store them in a list
samples_list <- list()
for (path in sample_paths) {
    srt_object <- readRDS(path)
    samples_list <- c(
        samples_list,
        list(srt_object)
    )
}
# Integrate the samples using the basic integrate function
# Set the sample ID
sample_ID = "IFNInt_2291401" # Integrating all IFN datasets
# Set the reference and then call the function
#reference = c(1, 2)
integrated_seurat <- basic_combine(
        list_of_srt_objects = samples_list,
        normalization_method = "LogNormalize", # But change accordingly
        sample_ID = sample_ID,
        output_directory = "2.Downstream_processing/2.1.Integration/IFNStim/2291401"#,
        #reference = reference
)
# Now - with the integrated samples - we run the rest of the downstream processing workflow
#Read in the merged dataset
integrated_seurat <- readRDS("2.Downstream_processing/2.1.Integration/project/sample_ID_merged.rds")
#Perform Dimensional reduction using the basic function
#        LogNormalised data
integrated_seurat <- ScaleData(integrated_seurat)
        ## sctransform normalised and LogNormalised data
integrated_seurat <- basic_reduction(
        seurat_object = integrated_seurat,
        sample_ID = sample_ID,
        output_directory = "./2.Downstream_processing/2.2.Dimension_reduction/IFNStim/2291401"
) # Assess PCA components before moving to next step
# Perform clustering steps using the basic function
#integrated_seurat <- readRDS("./2.Downstream_processing/2.2.Dimension_reduction/project/sample_ID_integrated_dim_reduced.rds")
integrated_seurat <- basic_cluster(
        seurat_object = integrated_seurat,
        sample_ID = sample_ID,
        output_directory = "2.Downstream_processing/2.3.Clustering/IFNStim/2291401",
        dimensions = 1:30,
        resolution = 0.5
)
# Create the conditions column prior to cluster markers
srt_clustered <- readRDS("2.Downstream_processing/2.3.Clustering/project/sample_ID_integrated_clustered.rds")
# Create the conditions - if needed
condition_one <- c(
        "UnS_227293",
        "UnS_227294"
)
condition_two <- c(
        "IFNA_227291",
        "IFNL_227292",
        "IFNB_227295"
)
seurat_object@meta.data <- mutate(
        seurat_object@meta.data,
        condition = case_when(
                orig.ident %in% condition_one ~ "UnStimulated",
                orig.ident %in% condition_two ~ "IFNStimulated"
        )
)
# Run known marker visualisations for cluster annottaions
known_marker_visualisations(
        seurat_object = integrated_seurat,
        sample_ID = sample_ID,
        split_by_col = "orig.ident",
        output_directory = "./2.Downstream_processing/2.4.Clustering_annotation/IFNStim/2291401",
        ncol = 3,
        colours = c(
                "#83d338",
                "#545cb7"
                # "#00b3ce",
                # "goldenrod1",
                # "#da429a"
        )
) # Assess cluster markers and annotate accordingly
# Annotate the clusters and save the annotated seurat object
integrated_seurat <- readRDS("2.Downstream_processing/2.3.Clustering/IFNStim/2291401/IFNInt_2291401_clustered.rds")
integrated_seurat@meta.data <- mutate(
        integrated_seurat@meta.data,
        celltype_major = case_when(
                seurat_clusters %in% c(
                        0, 1
                ) ~ "CD4T",
                seurat_clusters %in% c(
                        5, 7, 12
                ) ~ "CD8T",
                seurat_clusters %in% c(
                        8
                ) ~ "NK",
                seurat_clusters %in% c(
                        14, 15
                ) ~ "DCs",
                seurat_clusters %in% c(
                        6, 11
                ) ~ "Bcells",
                seurat_clusters %in% c(
                        2, 10
                ) ~ "Monocytes",
                seurat_clusters %in% c(
                        13, 16
                ) ~ "Platelets",
                seurat_clusters %in% c(
                        17, 18
                ) ~ "Unclassified",
                seurat_clusters %in% c(
                        3, 4
                ) ~ "Other",
                seurat_clusters %in% c(
                  9
                ) ~ "TSle"
        )
)
saveRDS(
        object = integrated_seurat,
        file = "2.Downstream_processing/2.4.Clustering_annotation/IFNStim/2291401/IFNInt_2291401_clustered_annotated.rds"
) # Saved: integrated, clustered and annotated
