# twosamplemr - official R TwoSampleMR with official PLINK2 offline
# instrument clumping, in one container run.
#
# Adapted from the legacy Rust wrapper (nodes-io twosamplemr_container.rs).
# Every parameter arrives as a TWOSAMPLEMR_-prefixed environment value;
# the optional `chr` renders empty for the all-autosomes default and the
# `clump` bool renders "true"/"false". No decompression stanza on
# purpose: the legacy wrapper never decompressed inputs (read.delim
# reads the table verbatim).

set -eu
: > "${AUTONOMICS_OUTPUT3}"

if [ "$TWOSAMPLEMR_CLUMP" = "true" ]; then
  # Cross-field rule the v0 param DSL cannot express (legacy validate()):
  # clump_p2 must be greater than or equal to clump_p1.
  if ! awk -v p1="$TWOSAMPLEMR_CLUMP_P1" -v p2="$TWOSAMPLEMR_CLUMP_P2" \
      'BEGIN { exit !(p2 >= p1) }'; then
    echo "clump_p2 must be greater than or equal to clump_p1" >&2
    exit 1
  fi

  # Exposure p-values when present, else the two-sided normal
  # approximation of beta/se - byte-identical to the legacy stage.
  Rscript --vanilla -e 'data <- read.delim(Sys.getenv("AUTONOMICS_INPUT0"), check.names = FALSE, stringsAsFactors = FALSE); p <- if (is.null(data$pval_exposure)) 2 * pnorm(-abs(data$beta_exposure / data$se_exposure)) else data$pval_exposure; write.table(data.frame(SNP = data$snp, P = p), "/work/clump_input.tsv", sep = "\t", quote = FALSE, row.names = FALSE)'

  mkdir -p /work/per_chr
  : > /work/all.clumps
  # Optional single-chromosome run: absence of TWOSAMPLEMR_CHR means
  # autosomes 1..=22 (legacy `seq 1 22` vs `seq <chr> <chr>`).
  CHR_SEQ="1 22"
  CHR_ARGS=""
  if [ -n "$TWOSAMPLEMR_CHR" ]; then
    CHR_SEQ="$TWOSAMPLEMR_CHR $TWOSAMPLEMR_CHR"
    CHR_ARGS="--chr $TWOSAMPLEMR_CHR"
  fi
  for chr in $(seq $CHR_SEQ); do
    out=/work/per_chr/chr${chr}
    plink2 \
       --bfile /panels/plink_ref/1000G.EUR.QC.${chr} \
       --clump /work/clump_input.tsv \
       --clump-id-field SNP --clump-p-field P \
       $CHR_ARGS \
       --clump-p1 "${TWOSAMPLEMR_CLUMP_P1}" \
       --clump-p2 "${TWOSAMPLEMR_CLUMP_P2}" \
       --clump-r2 "${TWOSAMPLEMR_CLUMP_R2}" \
       --clump-kb "${TWOSAMPLEMR_CLUMP_KB}" \
       --threads 1 \
       --out "$out" >>"${AUTONOMICS_OUTPUT3}" 2>&1
    if [ -f "$out.clumps" ]; then cat "$out.clumps" >> /work/all.clumps; fi
  done
fi

Rscript --vanilla - <<'RSCRIPT'
input <- Sys.getenv("AUTONOMICS_INPUT0")
table_path <- Sys.getenv("AUTONOMICS_OUTPUT0")
harmonised_path <- Sys.getenv("AUTONOMICS_OUTPUT1")
result_path <- Sys.getenv("AUTONOMICS_OUTPUT2")
log_path <- Sys.getenv("AUTONOMICS_OUTPUT3")
# Labels the legacy validate() rejected before the container started;
# the emptiness checks move into the script (coloc precedent).
if (trimws(Sys.getenv("TWOSAMPLEMR_ID_EXPOSURE")) == "") stop("id_exposure cannot be empty", call. = FALSE)
if (trimws(Sys.getenv("TWOSAMPLEMR_EXPOSURE")) == "") stop("exposure cannot be empty", call. = FALSE)
if (trimws(Sys.getenv("TWOSAMPLEMR_ID_OUTCOME")) == "") stop("id_outcome cannot be empty", call. = FALSE)
if (trimws(Sys.getenv("TWOSAMPLEMR_OUTCOME")) == "") stop("outcome cannot be empty", call. = FALSE)
data <- read.delim(input, check.names = FALSE, stringsAsFactors = FALSE)
pval <- if (is.null(data$pval_exposure)) 2 * pnorm(-abs(data$beta_exposure / data$se_exposure)) else data$pval_exposure
exposure_dat <- data.frame(SNP = data$snp, beta.exposure = data$beta_exposure, se.exposure = data$se_exposure, effect_allele.exposure = data$effect_allele_exposure, other_allele.exposure = data$other_allele_exposure, eaf.exposure = data$eaf_exposure, pval.exposure = pval, id.exposure = Sys.getenv("TWOSAMPLEMR_ID_EXPOSURE"), exposure = Sys.getenv("TWOSAMPLEMR_EXPOSURE"))
outcome_dat <- data.frame(SNP = data$snp, beta.outcome = data$beta_outcome, se.outcome = data$se_outcome, effect_allele.outcome = data$effect_allele_outcome, other_allele.outcome = data$other_allele_outcome, eaf.outcome = data$eaf_outcome, id.outcome = Sys.getenv("TWOSAMPLEMR_ID_OUTCOME"), outcome = Sys.getenv("TWOSAMPLEMR_OUTCOME"))
# Instrument selection: the clumped index SNPs when the PLINK2 stage
# ran, otherwise every input SNP (legacy clump_selection branches).
if (Sys.getenv("TWOSAMPLEMR_CLUMP") == "true") {
  clumps <- read.delim("/work/all.clumps", check.names = FALSE, stringsAsFactors = FALSE)
  selected <- unique(clumps$ID)
  exposure_dat <- exposure_dat[exposure_dat$SNP %in% selected, , drop = FALSE]
} else {
  selected <- unique(exposure_dat$SNP)
}
# method_list membership replaces the legacy Rust SUPPORTED_METHODS gate
# (the v0 param DSL has no enum type); the manifest minItems=1 covers the
# empty case.
method_list <- strsplit(Sys.getenv("TWOSAMPLEMR_METHOD_LIST"), " ", fixed = TRUE)[[1]]
supported <- c("mr_wald_ratio", "mr_two_sample_ml", "mr_egger_regression", "mr_egger_regression_bootstrap", "mr_simple_median", "mr_weighted_median", "mr_penalised_weighted_median", "mr_ivw", "mr_ivw_radial", "mr_ivw_mre", "mr_ivw_fe", "mr_simple_mode", "mr_weighted_mode", "mr_weighted_mode_nome", "mr_simple_mode_nome", "mr_sign", "mr_uwr", "mr_grip")
for (method in method_list) {
  if (!(method %in% supported)) stop(paste0("unsupported TwoSampleMR method `", method, "`"), call. = FALSE)
}
# Seed contract (WO-R-03): in v0.7.9 the stochastic method families all
# consume the global R RNG only inside mr() — mr_weighted_median and
# mr_simple_median via weighted_median_bootstrap (nboot = 1000), the mode
# family via its boot() (mr_mode.R), and mr_egger_regression_bootstrap via
# a vectorised rnorm matrix. harmonise_data reads no RNG. Therefore
# RNGkind + set.seed immediately before the mr() call below makes every
# bootstrap quantity deterministic for a given input, while deterministic
# methods (ivw, egger, wald ratio, mre/fe) are unaffected by the seed.
# The RNGkind triple mirrors the frozen host canonical runs recorded in the
# MRlap seed1 acceptance (Mersenne-Twister / Inversion / Rejection).
# An empty TWOSAMPLEMR_SEED reproduces the pre-fix unseeded behaviour.
seed_raw <- Sys.getenv("TWOSAMPLEMR_SEED")
seed_int <- NA_integer_
if (nzchar(seed_raw)) {
  seed_int <- suppressWarnings(as.integer(seed_raw))
  if (is.na(seed_int) || seed_int < 1L || seed_int > 2147483647L) {
    stop(paste0("seed must be an integer in 1..2147483647, got `", seed_raw, "`"), call. = FALSE)
  }
}
sink(log_path, split = TRUE, append = TRUE)
cat("TwoSampleMR:", as.character(packageVersion("TwoSampleMR")), "\n")
cat("Selected instruments:", length(unique(exposure_dat$SNP)), "\n")
cat("RNG seed:", if (nzchar(seed_raw)) seed_raw else "unset", "\n")
harmonised <- TwoSampleMR::harmonise_data(exposure_dat, outcome_dat, action = as.numeric(Sys.getenv("TWOSAMPLEMR_HARMONISE_ACTION")))
if (nzchar(seed_raw)) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed_int)
}
estimates <- TwoSampleMR::mr(harmonised, method_list = method_list)
print(estimates)
sink()
write.table(estimates, table_path, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(harmonised, harmonised_path, sep = "\t", quote = FALSE, row.names = FALSE)
result <- list(exposure_dat = exposure_dat, outcome_dat = outcome_dat, harmonised = harmonised, mr = estimates)
saveRDS(result, result_path)
RSCRIPT
