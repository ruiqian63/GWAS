library(EthSEQ)

# Set your VCF file path
#your_vcf <- "/neuro/labs/grantlab/research/enrique.mondragon/morton_lab/new_gwas/rawdata/merged_1.vcf"  # Replace with the actual path to your VCF file
your_vcf <- "/neuro/labs/grantlab/research/enrique.mondragon/morton_lab/new_gwas/rawdata/merged_folate.vcf"
#your_vcf <- "/Users/Rui/documents/gwas_all/filtered_point_duplicated.vcf"  # Replace with the actual path to your VCF file

# Get the built-in reference GDS file path
reference_gds <- system.file("extdata", "Reference.Gencode.Exome.10000SNPs.gds", package="EthSEQ")

# Check if the reference file exists
if (reference_gds == "") {
  stop("Built-in reference database not found. Please ensure the EthSEQ package is correctly installed.")
}

# Set the output directory
out_dir <- "/neuro/labs/grantlab/research/enrique.mondragon/morton_lab/new_gwas/EthSEQ/folate_result/"  # Replace with the directory path where you want to store the results
# Or use a temporary directory
# out_dir <- file.path(tempdir(), "EthSEQ_CustomAnalysis")

# Run the EthSEQ analysis
ethseq.Analysis(
  target.vcf = your_vcf,
  model.gds = reference_gds,
  out.dir = out_dir,
  verbose = TRUE,
  cores = 4,
  composite.model.call.rate = 1,
  space = "3D",
)

# Load and view the results
result_file <- file.path(out_dir, "Report_folate.txt")
if (!file.exists(result_file)) {
  stop("Analysis result file not found. Please check if the analysis completed successfully.")
}
ethseq_annotations <- read.delim(result_file, sep = "\t", as.is = TRUE, header = TRUE)
head(ethseq_annotations)

# Clean up temporary files (if using a temporary directory)
# unlink(out_dir, recursive = TRUE)
