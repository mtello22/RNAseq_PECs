# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if any arguments were passed
if (length(args) == 0) {
  stop("No filenames provided. Please provide filenames as arguments.")
}

# Print the filenames
cat("Filenames provided:\n")
print(args)