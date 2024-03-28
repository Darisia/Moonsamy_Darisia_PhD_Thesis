# This R script contains all of the code for the DESeq2
# .. analsyis of data derived from bulk RNA-sequencing

# The project accessions for the bulk RNA-seq datasets

# PRJNA694583
# PRJNA739760
# PRJNA916231
    # Note that this dataset encompasses the Huh7.5 wild type analysis
    #                                    the Huh7.5 STAT1 knockout analysis
    #                                    the Huh7.5 STAT2 Knockout analysis
    #                                    the Huh7.5 IRF9 Knockout analysis

# libraries I need to load
library(tidyverse)
library(Seurat) # I don't know why but it is just a habit
library(DESeq2)
library(tximport)
library(ggrepel)

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
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}


# Note: All of the bulk RNA-seq datasets underwent the same processing in order to obtain 
# .. differentially expressed genes. So the code below demonstrates the process for the 
# .. IFN-alpha stimulated THP-1 dataset (PRJNA694583)

##### Create the transcripts to genes mapping dataframe #####
# Create the tximport mapping file from the transcripts to genes mapping file
# .. that we already have
# first we need to read in the files
transcriptsToGenesfile <- read_tsv(
  "./transcripts_to_genes.txt",
  col_names = FALSE
) # This is the same transcripts_to_genes mapping file we used for scRNA-seq
# Convert to a data frame (just to make sure it's a dataframe type object)
transcriptsToGenesfile <- as.data.frame(transcriptsToGenesfile)
# The transcriptsToGenes mapping file we have contains three columns with the 
# ENST, ENSG and HGNC symbols - for this purpose, we just need the first two
transcriptsToGenesfile <- transcriptsToGenesfile[ , c(1,3)]
# Rename the columns - but we actually don't have to do it
colnames(transcriptsToGenesfile) <- c("TXNAME", "GENEID")

##### The TxImport function to summarise transcript counts into gene counts #####
# Providing the directory to each of the kallisto files
import = c(
  "./6.0.BulkRNASeq/24011103/24011105/abundance.h5", # THP-1 control 1
  "./6.0.BulkRNASeq/24011103/24011106/abundance.h5", # THP-1 control 2
  "./6.0.BulkRNASeq/24011103/24011107/abundance.h5", # THP-1 control 3
  "./6.0.BulkRNASeq/24011103/24011108/abundance.h5", # THP-1 IFNA-stimulated 1
  "./6.0.BulkRNASeq/24011103/24011109/abundance.h5", # THP-1 IFNA-stimulated 2
  "./6.0.BulkRNASeq/24011103/24011110/abundance.h5" # THP-1 IFNA-stimulated 3
)

# importing the files (Index prepared using --gencode and tx2gene using gtf so not need ignoreTxVersion/ignoreAfterBar)
txi.kallisto <- tximport(
  import, # Vector of file paths to the kallisto stored output
  type = "kallisto", # The type of alignment software used to generate the counts
  ignoreTxVersion = F, # whether or not to ignore what comes after the . in the ENSG names
  ignoreAfterBar = T, # whether to split the transcript ID to facilitate matching
  txOut = F, # Whether transcript level should be output
  tx2gene = transcriptsToGenesfile # The mapping file to convert transcript to gene IDs and aggregate the counts
)   
# 18143 transcripts missing from tx2gene when trying to summarise to gene count

##### Creating the meta data table in preparation for running DESeq2 #####
# Creating a table containing the sample information for each file (meta information)
MetaData <- data.frame(
  condition = factor(
    rep(
      c("Control", "IFNAStimulated"),
      each = 3
    )
  ),
  Sample = c(
    "R1",
    "R2",
    "R3",
    "R1",
    "R2",
    "R3"
  )
)       
# Naming for DESeq2 to read
rownames(MetaData) <- colnames(txi.kallisto$counts)
# Import into DESeq2
dds <- DESeqDataSetFromTximport(
  txi.kallisto,
  MetaData,
  ~condition
)

##### PCA plot #####
PCAVis <- plotPCA(
  vst(dds)
) +
  geom_point(
    size = 1
  ) +
  labs(
    colour = ""
  ) +
  scale_colour_manual(
    values = c(
      "#371ea3",
      "mediumaquamarine"
    )
  ) #+
# geom_text_repel(
#   aes(label = MetaData$Sample),
#   colour = "black"
# )
PCAVis <- edit_plots(PCAVis) +
  theme(
    legend.position = "top"
  )
PCAVis # It isn't visualising properly - not sure why - but it is geom_text_repel that is the problem
# Save this PCA plot - not saving with the replicate labels
ggsave(
  filename = "./6.0.BulkRNASeq/24011103/PCAVis.png",
  plot = PCAVis,
  width = 7,
  height = 5,
  dpi = 900
)
##### Running DeSeq2 #####
res <- DESeq(dds)
extractRes <- results(
  res,
  contrast = c(
    "condition",
    "IFNAStimulated",
    "Control"
  )
)
##### Extract the different data slots so that we can save them to file #####
ResDataFrame <- as.data.frame(extractRes@listData)
rownames(ResDataFrame) <- extractRes@rownames
ResDataFrame$genes <- rownames(ResDataFrame)
# Save the dataframe and the ddsObject
saveRDS(
  object = dds,
  file = "./6.0.BulkRNASeq/24011103/dds.rds"
)
saveRDS(
  object = res,
  file = "./6.0.BulkRNASeq/24011103/res.rds"
)
saveRDS(
  object = extractRes,
  file = "./6.0.BulkRNASeq/24011103/ExtractRes.rds"
)
saveRDS(
  object = ResDataFrame,
  file = "./6.0.BulkRNASeq/24011103/ResDataFrame.rds"
)

ResDataFramewCounts <- merge(
  ResDataFrame,
  as.data.frame(counts(res, normalized = TRUE)),
  by = "row.names"
)
# Save the Df with the normalised counts merged 
saveRDS(
  object = ResDataFramewCounts,
  file = "./6.0.BulkRNASeq/24011103/ResDataFramewCounts.rds"
)

##### Plot a preliminary Volcano plot #####
PrelimVolcano <- ggplot(
  ResDataFrame,
  aes(
    x = log2FoldChange,
    y = -log10(padj)
  )
) +
  geom_point(
    alpha = 0.5,
    size = 4
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    colour = "#da429a",
    linewidth = 1
  )
PrelimVolcano <- edit_plots(PrelimVolcano)
PrelimVolcano

# All bulk RNA-seq datasets from chapter were processed exactly as above

