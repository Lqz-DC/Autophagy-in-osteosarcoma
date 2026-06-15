################################################################################
## drug_independent_and_anchor.R
## INPUT : drugPred.csv
##         external_risk_labels.csv  (from external_validation.R)
##         "The average drug sensitivity of the four cell lines of osteosarcoma.csv.csv"
## OUTPUT: drug_independent_results.csv  + console summary
##
## Run from the repository root, AFTER external_validation.R:
##   Rscript drug_independent_and_anchor.R
################################################################################

suppressMessages(library(tidyverse))

DRUGS <- c("Gemcitabine", "Cytarabine", "Staurosporine", "Obatoclax Mesylate",
           "AZD4547", "Vinorelbine", "Vincristine", "Vinblastine", "Sepantronium bromide")

## ---- (A) independent-cohort drug check -------------------------------------
dp  <- read.csv("drugPred.csv", check.names = FALSE)        # DRUGS x samples
rownames(dp) <- dp$DRUGS; dp$DRUGS <- NULL
lab <- read.csv("external_risk_labels.csv", stringsAsFactors = FALSE)
risk <- setNames(lab$risk, lab$id)

cohort_samples <- function(coh) {
  s <- colnames(dp)
  if (coh == "GSE21257") s[grepl("^GSM", s)] else s[grepl("TARGET", s)]
}

out <- list()
for (coh in c("GSE21257", "TARGET")) {
  cols <- intersect(cohort_samples(coh), names(risk))
  hi <- cols[risk[cols] == "high"]; lo <- cols[risk[cols] == "low"]
  cat(sprintf("\n=== Independent drug check in %s  (risk from the other cohort's locked model; high=%d, low=%d) ===\n",
              coh, length(hi), length(lo)))
  cat(sprintf("%-22s%12s%12s%10s   %s\n", "drug", "median_high", "median_low", "p(MWU)", "more_sensitive"))
  for (d in DRUGS) {
    if (!d %in% rownames(dp)) next
    a <- as.numeric(dp[d, hi]); b <- as.numeric(dp[d, lo])
    p <- wilcox.test(a, b)$p.value
    sens <- if (median(a) < median(b)) "HIGH-risk" else "low-risk"   # lower lnIC50 = more sensitive
    cat(sprintf("%-22s%12.3f%12.3f%10.1e   %s\n", d, median(a), median(b), p, sens))
    out[[length(out) + 1]] <- data.frame(cohort = coh, drug = d,
        median_high = median(a), median_low = median(b), p = p, more_sensitive = sens)
  }
}
write.csv(do.call(rbind, out), "drug_independent_results.csv", row.names = FALSE)

## ---- (B) measured anchoring in the four OS cell lines ----------------------
avg <- read.csv("The average drug sensitivity of the four cell lines of osteosarcoma.csv.csv",
                check.names = FALSE)
colnames(avg) <- c("drug", "mean_LN_IC50")
avg <- avg[is.na(suppressWarnings(as.numeric(avg$drug))), ]          # drop numeric-id rows
avg$pct <- rank(avg$mean_LN_IC50) / nrow(avg) * 100                  # low percentile = more sensitive
cat("\n=== Measured sensitivity percentile across the four OS lines (lower IC50 = more sensitive) ===\n")
for (d in c("Gemcitabine", "Cytarabine", "Staurosporine", "Obatoclax Mesylate", "AZD4547")) {
  r <- avg[tolower(avg$drug) == tolower(d), ]
  if (nrow(r)) cat(sprintf("  %-20s mean lnIC50 = %.2f  ->  %.0fth percentile (of %d drugs)\n",
                           d, r$mean_LN_IC50[1], r$pct[1], nrow(avg)))
}
cat("\nNOTE: confirm whether this file holds GDSC-measured or Precily-predicted IC50 for the four lines.\n")
cat("For the strongest 'measured' claim, substitute raw per-line GDSC values\n")
cat("(GDSC2_fitted_dose_response: LN_IC50 for HOS, MG63, U2OS, SAOS2; drug = gemcitabine, cytarabine).\n")
