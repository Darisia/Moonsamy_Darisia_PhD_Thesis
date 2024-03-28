# This Python script contains the code to perform single-cell RNA sequencing data preprocessing with SASCRiP
# This code includes processing of:
#   Unstimulated1 PBMC
#   Unstimulated2 PBMC
#   IFN-alpha
#   IFN-lamda 
#   IFN-beta

# Prior to the processing of the data, the packages that will use needed to be imported
import SASCRiP
from SASCRiP import sascrip_functions
import os

#####################################################
# Assess the quality of the sequences obtained
#####################################################

# Where are the fastq files stored?
fastq_dirs = [
"../Raw/PRJNA597786/IFNAStimFastq/fastqFiles/",
"../Raw/PRJNA597786/IFNLStimFastq/fastqFiles/",
"../Raw/PRJNA597786/UnTreatedFastq/fastqFiles/",
"../Raw/PRJNA597786/UnTreated2Fastq/fastqFiles/",
"../Raw/PRJNA597786/IFNBStimFastq/fastqFiles/"
]
# Where should the fastQC output be stored
fastqc_output = [
"./0.fastqc/PRJNA597786/227291", # Data processing number for IFN-alpha
"./0.fastqc/PRJNA597786/227292", # Data processing number for IFN-lambda
"./0.fastqc/PRJNA597786/227293", # Data processing number for Unstimulated1
"./0.fastqc/PRJNA597786/227294", # Data processing numbeer for Unstimulated2
"./0.fastqc/PRJNA597786/227295" # Data processing number for IFN-beta
]
# run fastqc to assess the quality of the sequences
sascrip_functions.run_fastqc(
fastq_input_list = fastq_dirs,
fastqc_output_list = fastqc_output,
fastqc_options_list = ["--extract"]
)

#####################################################
# Demultiplexing and Pseudoalignemnt 
#####################################################
########### kallisto_bustools_count #################

# Create directories for the outputs of each sample
# Note: output from this function will have the same name - so if new folders for each sample are not created - each new sample will override the previous
sample_names = [
"./1.Preprocessing/IFNStimMain/227291/Counts_227291",
"./1.Preprocessing/IFNStimMain/227292/Counts_227292",
"./1.Preprocessing/IFNStimMain/227293/Counts_227293",
"./1.Preprocessing/IFNStimMain/227294/Counts_227294",
"./1.Preprocessing/IFNStimMain/227295/Counts_227295"
]
# Run the multi_kallisto_bustools_count function to pseduoalign and demultiplex sequencing reads
# Note: Not all samples are the same chemistry. The chemistries are as follows
# IFN-alpha - 10xv3
# IFN-lambda - 10xv3
# Unstimulated1 - 10xv3
# Unstimulated2 - 10xv1
# IFN-beta - 10xv1

sascrip_functions.kallisto_bustools_count(
    list_of_fastqs = fastq_dirs,
    single_cell_technology = ["10xv3", "10xv3", "10xv3", "10xv1", "10xv1"],
    output_directory_path = sample_names,
    species_index = "./1.Preprocessing/kallisto_index.idx", # Path to where the index is stored
    species_t2g = "./1.Preprocessing/transcripts_to_genes.txt", # Path to where the t2g file is stored
    main_output_directory = "./1.Preprocessing/IFNStimMain/"
    input_directories = True,
    read_separators = [
        ["R1", "R2"],
        ["R1", "R2"],
        ["R1", "R2"],
        ["R1", "R2", "R3"],
        ["R1", "R2", "R3"]
    ]
)
#

#####################################################
# Convert the counts into seurat_compatible counts
#####################################################
########### multi_seurat_matrix #################

# Create a list of file paths to the files needed for seurat_matrix
# list of the paths to the matrix file
matrix_files = [
    "./1.Preprocessing/IFNStimMain/227291/Counts_227291/Count_analysis/filtered_counts/filtered_counts.mtx",
    "./1.Preprocessing/IFNStimMain/227292/Counts_227292/Count_analysis/filtered_counts/filtered_counts.mtx",
    "./1.Preprocessing/IFNStimMain/227293/Counts_227293/Count_analysis/filtered_counts/filtered_counts.mtx",
    "./1.Preprocessing/IFNStimMain/227294/Counts_227294/Count_analysis/filtered_counts/filtered_counts.mtx",
    "./1.Preprocessing/IFNStimMain/227295/Counts_227295/Count_analysis/filtered_counts/filtered_counts.mtx"
]
# list of the paths to the gene file
gene_indices = [
    "./1.Preprocessing/IFNStimMain/227291/Counts_227291/Count_analysis/filtered_counts/filtered_counts.genes.txt",
    "./1.Preprocessing/IFNStimMain/227292/Counts_227292/Count_analysis/filtered_counts/filtered_counts.genes.txt",
    "./1.Preprocessing/IFNStimMain/227293/Counts_227293/Count_analysis/filtered_counts/filtered_counts.genes.txt",
    "./1.Preprocessing/IFNStimMain/227294/Counts_227294/Count_analysis/filtered_counts/filtered_counts.genes.txt",
    "./1.Preprocessing/IFNStimMain/227295/Counts_227295/Count_analysis/filtered_counts/filtered_counts.genes.txt"
]
# list of the paths to the barcode file
barcode_indices = [
    "./1.Preprocessing/IFNStimMain/227291/Counts_227291/Count_analysis/filtered_counts/filtered_counts.barcodes.txt",
    "./1.Preprocessing/IFNStimMain/227292/Counts_227292/Count_analysis/filtered_counts/filtered_counts.barcodes.txt",
    "./1.Preprocessing/IFNStimMain/227293/Counts_227293/Count_analysis/filtered_counts/filtered_counts.barcodes.txt",
    "./1.Preprocessing/IFNStimMain/227294/Counts_227294/Count_analysis/filtered_counts/filtered_counts.barcodes.txt",
    "./1.Preprocessing/IFNStimMain/227295/Counts_227295/Count_analysis/filtered_counts/filtered_counts.barcodes.txt"

]
# list of paths for directory to store the output (same folder as the current matrix files)
output_dirs = [
"./1.Preprocessing/IFNStimMain/227291/Counts_227291/Count_analysis/filtered_counts/",
"./1.Preprocessing/IFNStimMain/227292/Counts_227292/Count_analysis/filtered_counts/",
"./1.Preprocessing/IFNStimMain/227293/Counts_227293/Count_analysis/filtered_counts/",
"./1.Preprocessing/IFNStimMain/227294/Counts_227294/Count_analysis/filtered_counts/",
"./1.Preprocessing/IFNStimMain/227295/Counts_227295/Count_analysis/filtered_counts/"
]
# Now we call the function
sascrip_functions.multi_seurat_matrix(
    matrix_files = matrix_files,
    gene_indices = gene_indices,
    barcode_indices = barcode_indices,
    output_directories = output_dirs,
    t2g_file = "1.Preprocessing/transcripts_to_genes.txt",
    add_hgnc = True
)
## Output files will be stored as:


#####################################################
# Perform cell quality control and view the overall features of the data 
#####################################################
########### multi_run_cqc #################
#
# First we create the sample IDs (these will prefix all file names and will serve as their origin identity names in the seurat object)
sample_IDs = [
"IFNA_227291",
"IFNL_227292",
"UnS_227293", # Unstimulated1
"UnS_227294", # Unstimulated2
"IFNB_227295"
]
# Paths where the output of the cqc function will be stored
output_paths = [
"./1.Preprocessing/IFNStimMain/227291",
"./1.Preprocessing/IFNStimMain/227292",
"./1.Preprocessing/IFNStimMain/227293",
"./1.Preprocessing/IFNStimMain/227294",
"./1.Preprocessing/IFNStimMain/227295"
]

sascrip_functions.multi_run_cqc(
    input_files_or_folders = output_dirs, # the input here is the output from the previous function
    all_samples_ID = sample_IDs,
    output_directory_paths = output_paths,
    main_output_directory = "./1.Preprocessing/IFNStimMain/",
    output_prefix = "IFNStimAll",
    gene_column = 2
)
#
# Run the normalisation step to normalise the raw gene counts from good-quality cells
# for i in range(0, len(dataset_sample_IDs)):
#     filename = dataset_sample_IDs[i] + "_subset_seurat.rds"
#     seurat_object_path = os.path.join(
#     "1.Preprocessing/PRJNA646704/Cell_Quality_Control",
#     filename
#     )
#     sascrip_functions.sctransform_normalize(
#     seurat_object = seurat_object_path,
#     sample_ID = dataset_sample_IDs[i],
#     output_directory_path = "./1.Preprocessing/PRJNA646704/Normalised/"
# )