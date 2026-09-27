# twosamplemr_harmonise - official TwoSampleMR::read_exposure_data +
# read_outcome_data + harmonise_data on two separate summary-statistics
# files, emitting the merged snake_case TSV the `twosamplemr` node
# consumes on its input port.
#
# Adapted from the legacy Rust wrapper (nodes-io
# twosamplemr_harmonise_container.rs). The command stays ["Rscript"] and
# the interpreter executes this file directly; every parameter arrives as
# a TWOSAMPLEMR_-prefixed environment value. No decompression stanza on
# purpose: the legacy wrapper never decompressed inputs.

exp_path <- Sys.getenv("AUTONOMICS_INPUT0")
out_path <- Sys.getenv("AUTONOMICS_INPUT1")
harm_path <- Sys.getenv("AUTONOMICS_OUTPUT0")
log_path <- Sys.getenv("AUTONOMICS_OUTPUT1")
# Labels the legacy validate() rejected before the container started;
# the emptiness checks move into the script (coloc precedent).
if (trimws(Sys.getenv("TWOSAMPLEMR_ID_EXPOSURE")) == "") stop("id_exposure cannot be empty", call. = FALSE)
if (trimws(Sys.getenv("TWOSAMPLEMR_EXPOSURE")) == "") stop("exposure cannot be empty", call. = FALSE)
if (trimws(Sys.getenv("TWOSAMPLEMR_ID_OUTCOME")) == "") stop("id_outcome cannot be empty", call. = FALSE)
if (trimws(Sys.getenv("TWOSAMPLEMR_OUTCOME")) == "") stop("outcome cannot be empty", call. = FALSE)
sink(log_path, split = TRUE)
cat("TwoSampleMR:", as.character(packageVersion("TwoSampleMR")), "\n")
exp <- TwoSampleMR::read_exposure_data(exp_path)
out <- TwoSampleMR::read_outcome_data(out_path)
cat("exposure rows:", nrow(exp), " columns:", paste(names(exp), collapse = ","), "\n")
cat("outcome rows:", nrow(out), " columns:", paste(names(out), collapse = ","), "\n")
exp$id.exposure <- Sys.getenv("TWOSAMPLEMR_ID_EXPOSURE")
exp$exposure <- Sys.getenv("TWOSAMPLEMR_EXPOSURE")
# Optional free-text units labels: nzchar() gates the legacy
# Option<String> branches.
if (nzchar(Sys.getenv("TWOSAMPLEMR_UNITS_EXPOSURE"))) {
  exp$units.exposure <- Sys.getenv("TWOSAMPLEMR_UNITS_EXPOSURE")
}
out$id.outcome <- Sys.getenv("TWOSAMPLEMR_ID_OUTCOME")
out$outcome <- Sys.getenv("TWOSAMPLEMR_OUTCOME")
if (nzchar(Sys.getenv("TWOSAMPLEMR_UNITS_OUTCOME"))) {
  out$units.outcome <- Sys.getenv("TWOSAMPLEMR_UNITS_OUTCOME")
}
harm <- TwoSampleMR::harmonise_data(exp, out, action = as.numeric(Sys.getenv("TWOSAMPLEMR_HARMONISE_ACTION")))
cat("harmonised rows:", nrow(harm), " columns:", paste(names(harm), collapse = ","), "\n")
if (nrow(harm) == 0) stop("harmonise_data produced 0 rows; check SNP overlap between exposure and outcome")
# The downstream wrapper reads data$snp (lowercase) while harmonise_data
# emits SNP; the rewrite must cover both cases.
names(harm)[names(harm) == "SNP"] <- "snp"
names(harm) <- gsub("\\.", "_", names(harm))
write.table(harm, harm_path, sep = "\t", quote = FALSE, row.names = FALSE)
sink()
