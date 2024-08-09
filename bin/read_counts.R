library(data.table)
library(stringr)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# File path and parent folder name
file_path <- args[1]
parent_folder <- args[2]

# Check if the file path is valid
if (file_path == "") {
  stop("No filenames provided. Please provide filenames as arguments.")
}

# Read count results and extract the required columns
count_table <- fread(file_path)
count_table <- count_table[, .SD, .SDcols = c("Name", "NumReads")]
names(count_table) <- c("GeneID", parent_folder)

# Construct output filename
out_file <- file.path(dirname(file_path), paste0(parent_folder, ".tsv"))

# Write the combined columns to an output file
fwrite(x = count_table, file = out_file, 
       append = FALSE, quote = FALSE, 
       sep = '\t', row.names = FALSE, 
       col.names = TRUE)

# Print a message indicating successful processing
cat("Processed file:", file_path, "Output written to:", out_file, "\n")