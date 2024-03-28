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
# The format of the output files that will likely be used for the next steps is as follows:
# .. "output_path/sample_ID_subset_seurat.rds"
# Furthermore, a preQC seurat object is also saved so users are able to perform their own QC and visualisation if need be


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

# Gonna use the new multi_log_normalise function to normalise the gene expression information
# Since it's still new - gonna add the full code here (in case it hasn't been fully updated)
#####################################################
# Function to create the multi_log_noramlise code 
#####################################################
# lognormalise function that Runs LogNormalize from within Seurat to normalise the single-cell data for seuqencing depth
def log_normalize(
    seurat_object,
    sample_ID,
    output_directory_path = "working_directory",
    output_log_matrix = False,
    output_count_matrix = False):

    '''
    LogNormalise function within Seurat to perform normalisation on the raw counts

    Parameters:
    seurat_object(str): Path to the saved seurat object
    sample_ID(str): Name of the sample
    output_directory_path(str): Path to the output directory where all generated files and folders will be saved
    ENSG_gname38_path(str): Path to the file called ENSG_gname38.tsv that will allow us to convert ENSG name into hgnc symbols

    '''

    # Sort out the output directory
    if output_directory_path == "working_directory":
        output_directory_path = "./"
    else:
        gen_func.mkdirs(output_directory_path)

    # Generate the directory where the output log mtx matrix would be saved
    if output_log_matrix == True:
        output_matrix_dir = os.path.join(output_directory_path, sample_ID + '_log_normalised_matrix')
        gen_func.mkdirs(output_matrix_dir)

    # Generate the directory where the output count mtx matrix would be saved
    if output_log_matrix == True:
        output_matrix_dir = os.path.join(output_directory_path, sample_ID + '_normalised_matrix')
        gen_func.mkdirs(output_matrix_dir)

    # if transcripts_to_genes supplied - create ensg_gname
    if transcripts_to_genes_file != None:
        # Path to the required R script
        Bash_file = pkg_resources.resource_filename('SASCRiP', 'cut_t2g.sh')

        command = 'sh {} {} {}'.format(
        Bash_file,
        transcripts_to_genes_file,
        output_directory_path
        )

        check_process = subprocess.run(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)

        # check to see if that variable has a returncode above 0 - if it does, print the error message
        if check_process.returncode != 0:
            print(check_process.stderr.decode())
            print(check_process.stdout.decode())

        # Print all standard output to log file
        logger.info(check_process.stdout.decode())

        ensg_gname_path = os.path.join(output_directory_path, "ensg_gname.tsv")
    else:
        logger.warning("No transcripts_to_genes_file given. If ESEMBL gene names are used, a Seurat object cannot be properly generated and an error will be returned")
        ensg_gname_path = "None"

    # all variables for the bash code
    # 1) path to seurat object
    # 2) sample_ID
    # 3) output_directory_path
    # 4) ENSG_gname_path

    # Path to the required R script
    R_file = pkg_resources.resource_filename('SASCRiP', 'logNorm_seurat.R')
    #R_file = "./logNorm_seurat.R" # This needs to be in the current working durectory
    # let's put together the bash code
    command = "Rscript {} {} {} {} {} {} {} {} {} {}".format(
        R_file,
        seurat_object,
        sample_ID,
        output_directory_path,
        ensg_gname_path,
        output_log_matrix,
        output_count_matrix)

    # Use subprocess to run the bash command through python
    check_process = subprocess.run(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)

    # Check to see if there was an error. If so, print the standard error
    if check_process.returncode != 0:
        print(check_process.stderr)

####################################################################################################

# Defines the function that runs sctransform_noramlize on multiple samples simultaneously
def multi_log_normalize(
    seurat_objects,
    all_sample_IDs,
    output_directory_paths,
    transcripts_to_genes_file = None,
    output_log_matrix = False,
    output_count_matrix = False

):
    for i in range(0, len(all_sample_IDs)):
        sascrip_functions.log_normalize(
            seurat_object = seurat_objects[i],
            sample_ID = all_sample_IDs[i],
            output_directory = output_directory_paths[i] if isinstance(output_directory_paths, list) == True else output_directory_paths,
            output_log_matrix = output_log_matrix[i] if isinstance(output_log_matrix, list) == True else output_log_matrix,
            output_count_matrix = output_count_matrix[i] if isinstance(output_count_matrix, list) == True else output_count_matrix,
            transcripts_to_genes_file = transcripts_to_genes_file
        )

#####################################################
# Normalise the gene counts relative to the sequencing depth of each cell
#####################################################


# Sample_IDs (as before - keep this constant)
sample_IDs = [
"IFNA_227291",
"IFNL_227292",
"UnS_227293", # Unstimulated1
"UnS_227294", # Unstimulated2
"IFNB_227295"
]
# Output_paths (as before - all stored together)
output_paths = [
"./1.Preprocessing/IFNStimMain/227291",
"./1.Preprocessing/IFNStimMain/227292",
"./1.Preprocessing/IFNStimMain/227293",
"./1.Preprocessing/IFNStimMain/227294",
"./1.Preprocessing/IFNStimMain/227295"
]

# Now we create a path to the filtered seurat objects from the previous functions
seurat_objects = [
    "./1.Preprocessing/IFNStimMain/227291/IFNA_227291_subset_seurat.rds",
    "./1.Preprocessing/IFNStimMain/227292/IFNL_227292_subset_seurat.rds",
    "./1.Preprocessing/IFNStimMain/227293/UnS_227293_subset_seurat.rds",
    "./1.Preprocessing/IFNStimMain/227294/UnS_227294_subset_seurat.rds",
    "./1.Preprocessing/IFNStimMain/227295/IFNB_227295_subset_seurat.rds"
]
# Now we call the function
multi_log_normalize(
    seurat_objects = seurat_objects,
    all_sample_IDs = sample_IDs,
    output_directory_paths = output_paths
)

# The final preprocessing output will look as follows:
# .. "output_path/sample_ID_logNorm.rds" - these files will be used for downstream analysis

#####################################################
# End of preprocessing - next steps in R
#####################################################