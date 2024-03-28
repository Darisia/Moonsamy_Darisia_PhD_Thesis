# This R script contains the code for the classification of differentially 
# .. expressed genes from IFNa, IFNl, IFNb comparisons

# libraries I need to load
library(Seurat)
library(tidyverse)
library(circlize)
library(ComplexHeatmap)
library(ggupset)
library(ggvenn)
#library(UpSetR)

# Functions we need that we have previously created
# Edit plots function
edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          #plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}
# The function that is required for the upset plot
Convert_up <- function(TB){
  TB_tmp = TB
  for(i in colnames(TB)){
    TB_tmp[i] = i
  }
  TB_tmp[TB == 0] = ""
  TB_t <- data.frame(t(TB_tmp), stringsAsFactors = F)
  TB_tmp$upset = as.list(TB_t)
  return(TB_tmp)
}

##### Read in the datasets that we need #####
# Cell types we want to work with
KeepCellTypes <- c(
  "Bcells",
  "CD4T",
  "CD8T",
  "Monocytes",
  "NK"
) # As listed when we annotated (as originally stored in celltype_major)
# But now what about the re-analysis using the proper doublet count?
# IFN-alpha DEGs
IFNaDEG <- readRDS(
  "./4.DEG_analysis/2291402/IFNAvUnS_2291402_DEGs.rds"
)
# IFN-beta DEGs
IFNbDEG <- readRDS(
  "./4.DEG_analysis/2291404/IFNBvUnS_2291404_DEGs.rds"
)
# We also need the IFNlDEGs since we are originally working with both
IFNlDEG <- readRDS(
  "./4.DEG_analysis/2291403/IFNLvUnS_2291403_DEGs.rds"
)
# Sort out the IFNaDEG so it's easier to work with
CtNamesExtracted <- c()
for ( Ct_i in seq(1, length(IFNaDEG)) ) {
  DEGName <- names(IFNaDEG)[Ct_i]
  # extract the celltype name only from the comparison name
  CellTypeOnly <- unlist(
    str_split(
      DEGName,
      pattern = "_"
    )
  )[1]
  # Now add it to the open vector
  CtNamesExtracted <- c(
    CtNamesExtracted,
    CellTypeOnly
  )
} #lookfor
# re name the celltypes in the DEGs list
names(IFNaDEG) <- CtNamesExtracted
# Now, using the cell type names - extract Ct of interest
IFNaDEG <- IFNaDEG[
  which(
    names(IFNaDEG) %in% KeepCellTypes
  )
]
# Add celltype info to the dataframes within the DEG list
for ( i in KeepCellTypes ) {
  IFNaDEG[[i]]$CellType <- i
}
# Sort out the IFNbDEG so it's easier to work with
CtNamesExtracted <- c()
for ( Ct_i in seq(1, length(IFNbDEG)) ) {
  DEGName <- names(IFNbDEG)[Ct_i]
  # extract the celltype name only from the comparison name
  CellTypeOnly <- unlist(
    str_split(
      DEGName,
      pattern = "_"
    )
  )[1]
  # Now add it to the open vector
  CtNamesExtracted <- c(
    CtNamesExtracted,
    CellTypeOnly
  )
} #lookfor
# re name the celltypes in the DEGs list
names(IFNbDEG) <- CtNamesExtracted
# Now, using the cell type names - extract Ct of interest
IFNbDEG <- IFNbDEG[
  which(
    names(IFNbDEG) %in% KeepCellTypes
  )
]
# Add celltype info to the dataframes within the DEG list
for ( i in KeepCellTypes ) {
  IFNbDEG[[i]]$CellType <- i
}
# Sort out the IFNlDEG so it's easier to work with
CtNamesExtracted <- c()
for ( Ct_i in seq(1, length(IFNlDEG)) ) {
  DEGName <- names(IFNlDEG)[Ct_i]
  # extract the celltype name only from the comparison name
  CellTypeOnly <- unlist(
    str_split(
      DEGName,
      pattern = "_"
    )
  )[1]
  # Now add it to the open vector
  CtNamesExtracted <- c(
    CtNamesExtracted,
    CellTypeOnly
  )
} #lookfor
# re name the celltypes in the DEGs list
names(IFNlDEG) <- CtNamesExtracted
# Now, using the cell type names - extract Ct of interest
IFNlDEG <- IFNlDEG[
  which(
    names(IFNlDEG) %in% KeepCellTypes
  )
]
# Add celltype info to the dataframes within the DEG list
for ( i in KeepCellTypes ) {
  IFNlDEG[[i]]$CellType <- i
} # Note to self: The dataset is listed 7 times for this 'sorting' section of code

###########################################################
# C L A S S I F Y    D E G S
###########################################################

##### Classify DEGs according to filtering parameters #####
# What we can do before running the loop is put all datasets into one large list 
# .. and actually just save it like that
IFNResponseList <- list(
  "IFN-α" = IFNaDEG,
  "IFN-β" = IFNbDEG,
  "IFN-λ" = IFNlDEG
)
rm(
  IFNaDEG,
  IFNbDEG,
  IFNlDEG
) # So we free up some space
for (IFNi in seq(1, length(IFNResponseList))) {
  for (Ct in KeepCellTypes) {
    IFNResponseList[[IFNi]][[Ct]] <- filter(
      IFNResponseList[[IFNi]][[Ct]],
      !(pct.1 < 0.1 & pct.2 < 0.1) # Percentage of cells that detect genes > 10% in at least one conditions
    )
    IFNResponseList[[IFNi]][[Ct]] <- filter(
      IFNResponseList[[IFNi]][[Ct]],
      abs(avg_log2FC) > 0.5 # The absolute log2 fold change value > 0.5
    )
    IFNResponseList[[IFNi]][[Ct]] <- filter(
      IFNResponseList[[IFNi]][[Ct]],
      p_val_adj < 0.05 # Adjusted p-value < 0.05
    )
  }
} # If genes satisfy all three thresholds - they will be classified as differentially expressed
# Save the overall Response List
saveRDS(
  object = IFNResponseList,
  file = "./A24011001/IFNResponseList.rds" # A24011001 = Analysis Code where we can find these files stored (for records)
) # Then we can use this every time we want to work with this information

###########################################################
# V I S U A L I S E    D E G    N U M B E R S
###########################################################

##### Subset the IFNResponseList and create the dataframe #####
IFNResponseList <- IFNResponseList[1:3]

# Combine the dataframes within the list
# Create new list with all the cell types combined for each treatment separately
IFNablCombinedDEGs <- list()
for ( i in seq(1, length(IFNResponseList)) ) {
  IFNCombined <- bind_rows(
    IFNResponseList[[i]]
  )
  # Add this combined df to the list opened before the loop
  IFNablCombinedDEGs <- c(
    IFNablCombinedDEGs,
    list(IFNCombined)
  )
}
names(IFNablCombinedDEGs) <- c("IFN-α", "IFN-β", "IFN-λ")
for ( i in seq(1, length(IFNablCombinedDEGs)) ) {
  IFNablCombinedDEGs[[i]]$Treatment <- names(IFNablCombinedDEGs)[i]
}
# Now create one extra large dataframe that holds all the information that we need for plotting
IFNablCombinedDEGsDf <- bind_rows(
  IFNablCombinedDEGs
)
# Now use the dataframe above to group_by and count to make the bar plot of responsive genes
IFNablDECount <- IFNablCombinedDEGsDf %>%
  group_by(
    Treatment,
    CellType
  ) %>%
  summarise(
    TotalDE = n()
  )
# Now create the bar plot
IFNablDECountBar <- ggplot(
  IFNablDECount,
  aes(
    x = CellType,
    y = TotalDE,
    fill = Treatment
  )
) +
  geom_bar(
    stat = "identity",
    position = "dodge",
    colour = "black"
  ) +
  labs(
    x = "",
    y = "Number of differentially\nexpressed genes",
    fill = ""
  ) +
  scale_y_continuous(
    expand = c(0,0),
    limits = c(0, 1000)
  ) +
  scale_x_discrete(
    breaks = KeepCellTypes,
    labels = c(
      "B cells",
      bquote(CD4^'+'~ 'T cells'),
      bquote(CD8^'+'~ 'T cells'),
      "Monocytes",
      "NK cells"
    )
  ) +
  scale_fill_manual(
    values = c(
      "#371ea3",
      "#b3eebe",
      "grey30"
    ),
    labels = c(
      "IFN-α",
      "IFN-β",
      "IFN-λ"
    )
  )
IFNablDECountBar <- edit_plots(IFNablDECountBar) +
  theme(
    axis.text.x = element_text(
      size = 16,
      colour = "black",
      angle = 90
    ),
    axis.title.y = element_text(
      size = 16,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 14,
      colour = "black"
    ),
    legend.text = element_text(
      size = 16
    )
  )
IFNablDECountBar
# Now save this plot
ggsave(
  filename = "./A24020703/IFNablDECountBar.png",
  plot = IFNablDECountBar,
  width = 9,
  height = 7,
  dpi = 900
)

###########################################################
# S H A R E D/U N I Q U E    T Y P E I    D E G S    A C R O S S    C E L L T Y P E S
###########################################################
# Fist we plot the upSet plot
# We need to create AlphaVecs and BetaVecs
# AlphaVecs List
AlphaVecs <- list()
# BetaVecs List
BetaVecs <- list()
# Open a for loop so we can extract the vectors of genes for all cell types
for (Ct in KeepCellTypes) {
  # Start with Alpha
  Alpha <- IFNResponseList[[1]][[Ct]]$genes
  AlphaVecs <- c(
    AlphaVecs,
    list(Alpha)
  )
  # And then Beta
  Beta <- IFNResponseList[[2]][[Ct]]$genes
  BetaVecs <- c(
    BetaVecs,
    list(Beta)
  )
}
# Add cell type names to the vector list
names(AlphaVecs) <- KeepCellTypes
names(BetaVecs) <- KeepCellTypes
# Now we need to convert the list to a matrix
# Note: need complex heatmap for this
# Alpha
AlphaVecsMat <- list_to_matrix(AlphaVecs)
# Beta
BetaVecsMat <- list_to_matrix(BetaVecs)
# From the matrix - we need to create a dataframe for each
# Alpha
AlphaVecsMatDf <- as.data.frame(AlphaVecsMat)
# Beta
BetaVecsMatDf <- as.data.frame(BetaVecsMat)
# Use the Convert_up function to convert the Mat dataframe into an upset-compatable
# .. dataframe
# Alpha
AlphaVecsConvert <- Convert_up(AlphaVecsMatDf)
# Beta
BetaVecsConver <- Convert_up(BetaVecsMatDf)
# Add identifiers to each of the above dataframes
# Stimulated with IFN-alpha
AlphaVecsConvert$Stimulant <- "IFN-α"
# Stimulated with IFN-beta
BetaVecsConver$Stimulant <- "IFN-β"
# Now, since they have identifies we can combine these two dataframes into one
AlphaBetaVecsConvert <- rbind(
  AlphaVecsConvert,
  BetaVecsConver
)
# Now that we have this dataframe - we can create the plot
AlphaBetaResponseUpset <- ggplot(
  AlphaBetaVecsConvert,
  aes(
    x = upset
    ,
    fill = Stimulant
  )
) +
  geom_bar(
    # position = "dodge",
    position = position_dodge2(width = 0.9, preserve = "single")
    ,
    colour = "black"
    # ,
    # fill = "#b3eebe"
  ) +
  labs(
    x = "",
    y = "Number of genes",
    fill = ""
  ) +
  #geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) +
  scale_x_upset() +
  scale_y_continuous(
    expand = c(0,0)
  ) +
  geom_hline(
    yintercept = 0
  ) +
  scale_fill_manual(
    values = c(
      "#371ea3",
      "#b3eebe"
    )
  ) #+
# facet_wrap(
#   ~Regulation,
#   nrow = 2
# )
AlphaBetaResponseUpset
AlphaBetaResponseUpset <- edit_plots(
  AlphaBetaResponseUpset
)
AlphaBetaResponseUpset
# Save the upset plot
ggsave(
  filename = "./A24011001/AlphaBetaResponseUpset.png",
  plot = AlphaBetaResponseUpset,
  width = 12,
  height = 5,
  dpi = 900
) # UpSet Plot illustrating shared and unique gene sets for IFN-alpha and IFN-beta

###########################################################
# C O R E    T Y P E I    D E G S    A C R O S S    C E L L T Y P E S
###########################################################
# First the core IFNa response genes
CoreIFNa <- IFNResponseList[[1]][["Bcells"]]$genes[
  IFNResponseList[[1]][["Bcells"]]$genes %in% 
    IFNResponseList[[1]][["CD4T"]]$genes
]
CoreIFNa <- CoreIFNa[
  CoreIFNa %in%
    IFNResponseList[[1]][["CD8T"]]$genes
]
CoreIFNa <- CoreIFNa[
  CoreIFNa %in%
    IFNResponseList[[1]][["Monocytes"]]$genes
]
CoreIFNa <- CoreIFNa[
  CoreIFNa %in%
    IFNResponseList[[1]][["NK"]]$genes
] # There are 129 genes in this vector

# Now the core IFNb response genes
CoreIFNb <- IFNResponseList[[2]][["Bcells"]]$genes[
  IFNResponseList[[2]][["Bcells"]]$genes %in% 
    IFNResponseList[[2]][["CD4T"]]$genes
]
CoreIFNb <- CoreIFNb[
  CoreIFNb %in%
    IFNResponseList[[2]][["CD8T"]]$genes
]
CoreIFNb <- CoreIFNb[
  CoreIFNb %in%
    IFNResponseList[[2]][["Monocytes"]]$genes
]
CoreIFNb <- CoreIFNb[
  CoreIFNb %in%
    IFNResponseList[[2]][["NK"]]$genes
] # There are 103 genes in this vector

# Now for those that are core across IFNa and IFNb conditions
CoreIFNab <- CoreIFNa[
  CoreIFNa %in%
    CoreIFNb
] # Only 86 genes are found here


###########################################################
# M O N O C Y T E - S P E C I F I C    T Y P E I    D E G S
###########################################################
# First in IFNa conditions
MonoSpIFNa <- IFNResponseList[[1]][["Monocytes"]]$genes[
  !(IFNResponseList[[1]][["Monocytes"]]$genes %in% 
    IFNResponseList[[1]][["Bcells"]]$genes)
]
MonoSpIFNa <- MonoSpIFNa[
  !(MonoSpIFNa %in% 
      IFNResponseList[[1]][["CD4T"]]$genes)
]
MonoSpIFNa <- MonoSpIFNa[
  !(MonoSpIFNa %in% 
      IFNResponseList[[1]][["CD8T"]]$genes)
]
MonoSpIFNa <- MonoSpIFNa[
  !(MonoSpIFNa %in% 
      IFNResponseList[[1]][["NK"]]$genes)
] # 543 monocyte specific genes instead of 1K and some change
# Now the IFNb condition
# First in IFNa conditions
MonoSpIFNb <- IFNResponseList[[2]][["Monocytes"]]$genes[
  !(IFNResponseList[[2]][["Monocytes"]]$genes %in% 
      IFNResponseList[[2]][["Bcells"]]$genes)
]
MonoSpIFNb <- MonoSpIFNb[
  !(MonoSpIFNb %in% 
      IFNResponseList[[2]][["CD4T"]]$genes)
]
MonoSpIFNb <- MonoSpIFNb[
  !(MonoSpIFNb %in% 
      IFNResponseList[[2]][["CD8T"]]$genes)
]
MonoSpIFNb <- MonoSpIFNb[
  !(MonoSpIFNb %in% 
      IFNResponseList[[2]][["NK"]]$genes)
] # 494 monocyte-specific genes in the IFNb condition
# Now let's check the overlap
MonoSpIFNab <- MonoSpIFNa[
  MonoSpIFNa %in% 
    MonoSpIFNb
] # Only 235 genes overlap
# Note that 35 of these genes are ribosomal protein genes and will be removed prior to visualisation


##### Create lists for the different Analysis Results #####
# This step is so I can save everything in a away that makes it 
# more easily accessible when I need to access the data for figures etc
# First it's the core signatures
IFNabCoreList <- list(
  "IFN-α" = CoreIFNa, # 129
  "IFN-β" = CoreIFNb, # 103
  "IFN-α/β" = CoreIFNab # 86
)
# Now let's save this
saveRDS(
  object = IFNabCoreList,
  file = "./A24011001/IFNabCoreList.rds"
)
# Now for the monocyte-specific signatures
MonoSpIFNabList <- list(
  "IFN-α" = MonoSpIFNa, # 543
  "IFN-β" = MonoSpIFNb, # 494
  "IFN-α/β" = MonoSpIFNab # 235 - 35 RP = 200
)
# Now let's save this too
saveRDS(
  object = MonoSpIFNabList,
  file = "./A24011001/MonoSpIFNabList.rds"
)