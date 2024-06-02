blocks <- read.table('/N/u/viyella/Quartz/Downloads/LD_blocks_ADNI_19.blocks.det', header = TRUE)

# Create a new directory for the .rds files
dir.create("rds_files")

# Loop over all the rows in the blocks dataframe
for (i in 1:nrow(blocks)) {
  # Extract the i-th block
  block <- blocks[i, ]
  
  # Split the SNPs string into a vector of SNP IDs
  snps <- strsplit(block$SNPS, "\\|")[[1]]
  
  # Extract the SNPs from the VCF file
  block_snps <- vcf[vcf@fix[, "ID"] %in% snps, ]
  
  samples <- data.frame(block_snps@gt)
  rsIDs <- data.frame(block_snps@fix)
  
  # Select the "ID", "REF", and "ALT" columns from the rsids data frame
  selected_columns <- rsIDs[, c("ID", "REF", "ALT")]
  
  combined_df <- cbind(samples, selected_columns)
  
  # Get the number of columns in the data frame
  num_columns <- ncol(combined_df)
  
  # Create a vector of column indices in the desired order
  new_order <- c((num_columns-2):num_columns, 1:(num_columns-3))
  
  # Rearrange the columns
  rearranged_df <- combined_df[, new_order]
  
  # Remove the fourth column
  Block_df <- rearranged_df[,-4]
  
  # Save the dataframe as an .rds file in the new directory
  saveRDS(Block_df, file = paste0("rds_files/Block", i, ".rds"))
}

