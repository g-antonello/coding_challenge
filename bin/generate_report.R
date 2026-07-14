#----------------------------------------------
# Pipeline report script to generate a csv file
# with required fields (see assignment)
#----------------------------------------------

#----------------------------------------------
# custom functions to filter output
# split FORMAT and per-sample genotype column together, keyed by name
parse_format_field <- function(format_col, sample_col, field) {
  format_keys <- strsplit(format_col, ":")
  sample_vals <- strsplit(sample_col, ":")
  mapply(function(keys, vals) {
    idx <- which(keys == field)
    if (length(idx) == 0) return(NA)
    vals[idx]
  }, format_keys, sample_vals)
}
#----------------------------------------------

library(optparse)

option.list <- list(
  make_option(c("-v", "--vcf_file"), type = "character", default = NULL),
  make_option(c("-r", "--ref_table"), type = "character", default = NULL),
  make_option(c("-o", "--output_filename"), type = "character", default = NULL),
  make_option(c("-d", "--min_read_depth"), type = "integer", default = 10),
  make_option(c("-f", "--min_allele_freq"), type = "double", default = 0.05)
)
opt_parser <- OptionParser(option_list = option.list)
opts <- parse_args(opt_parser)

#------------------------------------------
# sanity checks
if (is.null(opts$vcf_file)) {
  stop("no vcf file provided")
}
if (is.null(opts$ref_table)) {
  stop("no reference table provided")
}

# extract sample_id from vcf filename
sample_id <- basename(gsub(".vcf$", "", opts$vcf_file))
sample_id_makenames <- make.names(sample_id)

# default output filename based on sample_id
if (is.null(opts$output_filename)) {
  opts$output_filename <- sprintf("%s_report.csv", sample_id)
}

#------------------------------------------
# load vcf header separately
colnames <- grep("#CHROM", readLines(opts$vcf_file), value = TRUE)
colnames <- gsub("#", "", strsplit(colnames, "\t")[[1]], fixed = TRUE)
colnames.makenames <- make.names(colnames) 

# note: read.table automatically starts reading lines that do not start with #
vcf_full_noheaders <- read.table(file = opts$vcf_file, sep = "\t", header = FALSE, col.names = colnames.makenames)

# get Allele depth and frequency explicitly out of the sample and format fileds
vcf_full_noheaders$AD_raw <- parse_format_field(vcf_full_noheaders$FORMAT, vcf_full_noheaders[[sample_id_makenames]], "AD")
vcf_full_noheaders$AD <- sapply(strsplit(vcf_full_noheaders$AD_raw, ",", fixed = TRUE), "[", 2)
vcf_full_noheaders$AF     <- as.numeric(parse_format_field(vcf_full_noheaders$FORMAT, vcf_full_noheaders[[sample_id_makenames]], "AF"))
# get positions from genome reference table
positions_ref <- read.table(opts$ref_table, header = TRUE, sep = "\t")
positions_ref <- positions_ref[positions_ref$ID_in_fasta_reference_db %in% unique(vcf_full_noheaders$CHROM), ]

ref_pos <- c(positions_ref$A2058_PO, positions_ref$A2059_POS)

# initialize the report with fixed elements, relative to the vcf created anyway
report_1.0 <- data.frame(
  sample_id = rep(colnames[length(colnames)], 2),
  gene = rep(positions_ref$ID_in_fasta_reference_db, 2),
  ecoli_pos = c(2058, 2059),
  ref_pos = ref_pos,
  ref_base = c(positions_ref$A2058_REF, positions_ref$A2058_REF),
  stringsAsFactors = FALSE)

# subset the vcf with only the positions wanted from the reference rrl gene
vcf_subset <- vcf_full_noheaders[vcf_full_noheaders$POS %in% ref_pos,]

if(nrow(vcf_subset) == 0){
  # if mutation is not detected in the vcf, fill remaining fields
  # with NA
  report_2.0 <- report_1.0
  report_2.0$alt_base  <- NA
  report_2.0$alt_depth  <- NA
  report_2.0$alt_freq  <- NA
  report_2.0$mutation_detected <- NA
} else{
  # if at least one mutation is found, merge 1.0 with vcf_subset based on the 
  ## get coverage/seq depth in each position
  vcf_subset_for_merging <- data.frame(
    ref_pos = vcf_subset$POS,
    alt_base = vcf_subset$ALT,
    alt_depth = vcf_subset$AD,
    alt_freq =vcf_subset$AF,
    mutation_detected = paste0(vcf_subset$REF, vcf_subset$POS, vcf_subset$ALT)
  )
  
  report_2.0 <- merge(report_1.0, vcf_subset_for_merging, by = "ref_pos")
}

# final step is a decision on call. this is arbitrary and will be subject to change
report_2.0$decision <- ifelse(
  report_2.0$alt_depth > opts$min_read_depth & report_2.0$alt_freq  > opts$min_allele_freq,
  "PASS",
  ifelse(
    report_2.0$alt_depth <= opts$min_read_depth &
      report_2.0$alt_freq  <= opts$min_allele_freq,
    "FAIL",
    "WARN"
  )
)

# add notes about steps above.
# so far only the NA in decision is noted as "variant not found"
report_2.0$notes <- ifelse(is.na(report_2.0$decision), "VARIANT_NOT_FOUND", "")


# write report
write.csv(report_2.0, file = opts$output_filename, row.names = FALSE)