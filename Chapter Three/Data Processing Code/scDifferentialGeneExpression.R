# This R script contains the code to perform differential gene expression
# .. on the scRNA-seq datasets used in chapter three

# First we need to load all required packages
library(tidyverse)
library(Seurat)

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

#######################################################################
# D I F F E R E N T I A L    G E N E     E X P R E S S I O N    A N A L Y S I S
#######################################################################

#####################################################################################
# Function to find differentially expressed genes for each cell type
#####################################################################################
find_DE_ISGs <- function(
    seurat_object,
    sample_ID,
    #antiviral_csv,
    cell_type_column,
    group_1_values, # Which values (sample names) in orig.ident form the first group
    group_2_values, # Which values (sample names) in orig.ident form the second group
    output_directory,
    group_1_name = "Healthy", # The name of the first group (comparison group)
    group_2_name = "Stimulated" # The name of the second group (comparison group)
){
    # Set the default assay to "SCT"
    # DefaultAssay(seurat_object) <- "SCT"
    DefaultAssay(seurat_object) <- "RNA"
    # Create new column called condition
    seurat_object@meta.data <- mutate(
        seurat_object@meta.data,
        condition = case_when(
            orig.ident %in% group_1_values ~ group_1_name,
            orig.ident %in% group_2_values ~ group_2_name
        )
    )
    # Create another new column called type_condition
    # .. combining the condition and cell type columns
    seurat_object$type_condition <- paste(
        as.character(seurat_object@meta.data[[cell_type_column]]),
        seurat_object@meta.data[,'condition'],
        sep = "_"
    )
    # Set the identity to the newly created type_ocndition column
    Idents(seurat_object) <- 'type_condition'
    # Read in the isg information from csv
    # isg_information <- read_csv2(
    #     file = antiviral_csv
    # )
    # Extract only the gene names from this isg table
    #antiviral_isgs <- isg_information$`Gene Symbol`
    # Now we want to run the DEG tests
    # Create vector with the cell type or cluster information
    all_celltypes <- levels(as.factor(seurat_object@meta.data[,cell_type_column]))
    # Make sure that the celltypes we test have cells in both conditions
    celltypes <- c()
    for (celltype in all_celltypes ) {
        cond1_val <- sprintf(
            "%s_%s",
            celltype,
            group_1_name
        )
        cond2_val <- sprintf(
            "%s_%s",
            celltype,
            group_2_name
        )
        if ( cond1_val %in% seurat_object$type_condition && cond2_val %in% seurat_object$type_condition ) {
            celltypes <- c(
                celltypes,
                celltype
            )
        }
    }
    # Empty list to store the output dataframes after each FindMarkers run
    DEG_cell_specific_dfs <- list()
    # Empty vector to store names for each of the DEG dfs returned
    DEG_cell_specific_names <- c()
    # For loop to run the FindMarkers function
    for ( celltype in celltypes ) {
        ident_cond1 <- sprintf(
            "%s_%s",
            celltype,
            group_1_name
        )
        ident_cond2 <- sprintf(
            "%s_%s",
            celltype,
            group_2_name
        )
        # FindMarkers code
        DEG_df <- FindMarkers(
            seurat_object,
            ident.1 = ident_cond1,
            ident.2 = ident_cond2,
            logfc.threshold = 0,
            min.cells.group = 0,
            min.cells.feature = 0,
            min.pct = 0
        )
        # Create a column with the gene names which are currently stores as rownames
        DEG_df$genes <- rownames(DEG_df)
        # Find the ISGs within the dataframe (if there are any)
        # DEG_df$is_isg <- ifelse(
        #     DEG_df$genes %in% antiviral_isgs,
        #     TRUE,
        #     FALSE
        # )
        # Create the names for each of the DEG dataframes
        df_name <- sprintf(
            "%s_%s_vs_%s",
            celltype,
            group_1_name,
            group_2_name
        )
        # Add a column with the information about the comparison
        DEG_df$comparison <- df_name
        # Add a column with the differential expression classification
        DEG_df <- mutate(
            DEG_df,
            DEG_direction = case_when(
                avg_log2FC > 0.5 & p_val_adj < 0.05 ~ "Upregulated",
                avg_log2FC < -0.5 & p_val_adj < 0.05 ~ "Downregulated"
            )
        )
        DEG_df$classification <- ifelse(
            is.na(DEG_df$DEG_direction),
            "non_DE",
            DEG_df$DEG_direction
        )
        # Store in the DEG list
        DEG_cell_specific_dfs <- c(
            DEG_cell_specific_dfs,
            list(DEG_df)
        )
        # Store the dataframe names in the vector
        DEG_cell_specific_names <- c(
            DEG_cell_specific_names,
            df_name
        )
    }
    # Set the DEG df names for each of the dfs in the DEG list
    DEG_cell_specific_dfs <- setNames(
        DEG_cell_specific_dfs,
        DEG_cell_specific_names
    )
    # Create file path for the list of dfs so it can be saved as an R object
    DEG_save_name <- sprintf(
        "%s/%s_DEGs.rds",
        output_directory,
        sample_ID
    )
    # Save the R object list to file for easy access
    saveRDS(
        object = DEG_cell_specific_dfs,
        file = DEG_save_name
    )
    # Since we added some important columns to the seurat object meta data slot - we should save the seurat object
    # Create the file path for the new seurat object so it can be saved as an R object
    updated_seurat_object_name <- sprintf(
        "%s/%s_clustered_seurat.rds",
        output_directory,
        sample_ID
    )
    # Save the R object to file
    saveRDS(
        object = seurat_object,
        file = updated_seurat_object_name
    )
    ## One last thing - we want to find DE genes between all healthy vs stimulated PBMCs
    # we need to set the identity of the cells to the condition column
    # Idents(seurat_object) <- 'condition'
    # # FindMarkers code
    # all_pbmc_DEG <- FindMarkers(
    #     seurat_object,
    #     ident.1 = group_1_name,
    #     ident.2 = group_2_name,
    #     logfc.threshold = 0, # We need to include all genes in the output df so we can make the volcano plot
    #     min.pct = 0
    # )
    # Create a column for the gene names
    # all_pbmc_DEG$genes <- rownames(all_pbmc_DEG)
    # # Find the ISGs within the dataframe (if there are any)
    # all_pbmc_DEG$is_isg <- ifelse(
    #     all_pbmc_DEG$genes %in% antiviral_isgs,
    #     TRUE,
    #     FALSE
    # )
    # # Classify the genes based on the significance thresholds
    # all_pbmc_DEG <- mutate(
    #     all_pbmc_DEG,
    #     classification = case_when(
    #             avg_log2FC > 0 & p_val_adj < 0.5 ~ "Upregulated",
    #             avg_log2FC < 0 & p_val_adj < 0.5 ~ "Downregulated",
    #             p_val_adj >= 0.5 ~ "non_DE"
    #         )
    # )
    # Create file path for the pbmc df so it can be saved as an R object
    # all_pbmc_DEG_name <- sprintf(
    #     "%s/%s_all_pbmc_DEG_df.rds",
    #     output_directory,
    #     sample_ID
    # )
    # # Save the R object to file
    # saveRDS(
    #     object = all_pbmc_DEG,
    #     file = all_pbmc_DEG_name
    # )
}


##########################################
# Run the differential gene expression analysis on the specified pairs
##########################################

##########################################
# DGEA on IFNA vs Unstim1
##########################################
# First we need to subset the integrated seurat object (since it contains information for literally everything)
Idents(integrated_seurat) <- integrated_seurat@meta.data$orig.ident
IFNAvUnSintSeurat <- subset(
  integrated_seurat,
  idents = c(
    "IFNA_227291",
    "UnS_227293"
  )
)
Idents(IFNAvUnSintSeurat) <- IFNAvUnSintSeurat@meta.data$seurat_clusters
# Now we only have 2 conditions and we can run the DEG function
find_DE_ISGs(
  seurat_object = IFNAvUnSintSeurat,
  sample_ID = "IFNAvUnS_2291402",
  cell_type_column = "celltype_major",
  group_1_values = c("IFNA_227291"),
  group_2_values = c("UnS_227293"),
  output_directory = "4.DEG_analysis/2291402/",
  group_1_name = "IFNaStimulated",
  group_2_name = "UnStimulated"
)
# Output stored in 2291402

##########################################
# DGEA on IFNL vs Unstim1
##########################################
# First we need to subset the integrated seurat object (since it contains information for literally everything)
Idents(integrated_seurat) <- integrated_seurat@meta.data$orig.ident
IFNLvUnSintSeurat <- subset(
  integrated_seurat,
  idents = c(
    "IFNL_227292",
    "UnS_227293"
  )
)
Idents(IFNLvUnSintSeurat) <- IFNLvUnSintSeurat@meta.data$seurat_clusters
# Now we only have 2 conditions and we can run the DEG function
find_DE_ISGs(
  seurat_object = IFNLvUnSintSeurat,
  sample_ID = "IFNLvUnS_2291403",
  cell_type_column = "celltype_major",
  group_1_values = c("IFNL_227292"),
  group_2_values = c("UnS_227293"),
  output_directory = "4.DEG_analysis/2291403/",
  group_1_name = "IFNlStimulated",
  group_2_name = "UnStimulated"
)
# Output stored in 2291403

##########################################
# DGEA on IFNB vs Unstim2
##########################################
# First we need to subset the integrated seurat object (since it contains information for literally everything)
Idents(integrated_seurat) <- integrated_seurat@meta.data$orig.ident
IFNBvUnSintSeurat <- subset(
  integrated_seurat,
  idents = c(
    "IFNB_227295",
    "UnS_227294"
  )
)
Idents(IFNBvUnSintSeurat) <- IFNBvUnSintSeurat@meta.data$seurat_clusters
# Now we only have 2 conditions and we can run the DEG function
find_DE_ISGs(
  seurat_object = IFNBvUnSintSeurat,
  sample_ID = "IFNBvUnS_2291404",
  cell_type_column = "celltype_major",
  group_1_values = c("IFNB_227295"),
  group_2_values = c("UnS_227294"),
  output_directory = "4.DEG_analysis/2291404/",
  group_1_name = "IFNbStimulated",
  group_2_name = "UnStimulated"
)
# Output stored in 2291404