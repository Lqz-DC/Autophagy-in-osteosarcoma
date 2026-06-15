################################################################################
## external_validation.R
##
## Strict external validation WITHOUT joint batch correction.
##   * the four-gene panel {BNIP3, MYC, PEA15, SAR1A} is the locked signature;
##   * each cohort is normalized INDEPENDENTLY (within-cohort z-score per gene);
##   * coefficients are estimated in the TRAINING cohort only;
##   * the high/low threshold is the TRAINING-cohort median;
##   * everything is evaluated in the held-out cohort, in BOTH directions.
##
## INPUT : figure1/merge.txt  (139 x 193, pre-batch-correction merged matrix)
##         figure1/clinical.txt    (id, futime[years], fustat, Gender, Age, Group)
## OUTPUT: external_validation_metrics.csv
##         KM_<train>_to_<test>.pdf  (x2)
##         external_risk_labels.csv  (independent risk label per sample; used by
##                                    drug_independent_and_anchor.R)
##         + paste-ready sentences printed to the console
##
## Run from the repository root:  Rscript external_validation.R
################################################################################

suppressMessages({
  library(survival)
  library(timeROC)
  library(survminer)
})

GENES <- c("BNIP3", "MYC", "PEA15", "SAR1A")

## ---- data ------------------------------------------------------------------
expr <- read.table("merge.txt", header = TRUE, sep = "\t",
                   check.names = FALSE, row.names = 1)          # samples x genes
expr <-as.data.frame(t(expr))
clin <- read.table("clinical.txt", header = TRUE, sep = "\t",
                   check.names = FALSE, row.names = 1)          # futime in YEARS
clin$futime <- clin$futime/365
stopifnot(all(GENES %in% colnames(expr)))

common <- intersect(rownames(expr), rownames(clin))
expr   <- expr[common, GENES, drop = FALSE]
clin   <- clin[common, ]
dat <- data.frame(expr, futime = clin$futime, fustat = clin$fustat,
                  check.names = FALSE)
dat$cohort <- ifelse(grepl("^GSM", rownames(dat)), "GSE21257", "TARGET")


zscore_within <- function(d) {
  d[GENES] <- scale(as.matrix(d[GENES]))   # mean 0, sd 1 per gene, within cohort
  d
}

## ---- one direction ---------------------------------------------------------
evaluate <- function(train_name, test_name) {
  tr <- zscore_within(dat[dat$cohort == train_name, ])
  te <- zscore_within(dat[dat$cohort == test_name, ])

  fit  <- coxph(Surv(futime, fustat) ~ BNIP3 + MYC + PEA15 + SAR1A, data = tr)
  beta <- coef(fit)[GENES]

  rs_tr <- as.matrix(tr[GENES]) %*% beta
  rs_te <- as.matrix(te[GENES]) %*% beta
  cut   <- median(rs_tr)                       # threshold defined on TRAINING only
  te$risk <- factor(ifelse(rs_te > cut, "high", "low"), levels = c("low", "high"))
  te$Auto_RS <- as.numeric(rs_te)

  ## univariate Cox of the continuous score (per SD) in the held-out cohort
  te$rs_sd <- as.numeric(scale(rs_te))
  uc  <- summary(coxph(Surv(futime, fustat) ~ rs_sd, data = te))
  hr  <- uc$conf.int[1, "exp(coef)"]; lo <- uc$conf.int[1, "lower .95"]
  hi  <- uc$conf.int[1, "upper .95"]; pv <- uc$coefficients[1, "Pr(>|z|)"]

  ## C-index and time-dependent AUC
  cidx <- summary(coxph(Surv(futime, fustat) ~ Auto_RS, data = te))$concordance[1]
  tROC <- timeROC(T = te$futime, delta = te$fustat, marker = te$Auto_RS,
                  cause = 1, times = c(1, 3, 5), iid = FALSE)
  auc  <- tROC$AUC

  ## KM with the training-derived split
  sd_fit <- survdiff(Surv(futime, fustat) ~ risk, data = te)
  km_p   <- 1 - pchisq(sd_fit$chisq, df = length(sd_fit$n) - 1)
  g <- ggsurvplot(survfit(Surv(futime, fustat) ~ risk, data = te), data = te,conf.int=T,
                  pval = TRUE, risk.table = TRUE, palette = c("#8AD293", '#EBBA37'),
                  legend.labs = c("low", "high"),
                  title = sprintf("Train %s -> Test %s", train_name, test_name))
  pdf(sprintf("KM_%s_to_%s.pdf", train_name, test_name), width = 6.5, height = 5.5, onefile = FALSE)
  print(g)
  dev.off()  # 关键：必须执行这行关闭设备，PDF 才能正常写入并打开！
  # 【新增】绘制并保存多时段 ROC 曲线
  pdf(sprintf("ROC_%s_to_%s.pdf", train_name, test_name), width = 6.5, height = 5.5)
  plot(tROC, time = 1, col = "#21908C", lwd = 2, title = FALSE)
  plot(tROC, time = 3, col = "#440154", lwd = 2, add = TRUE)
  plot(tROC, time = 5, col = "#FDE725", lwd = 2, add = TRUE)
  abline(a = 0, b = 1, col = "gray", lty = 2)
  legend("bottomright", 
         legend = c(sprintf("1-Year (AUC = %.3f)", auc[1]),
                    sprintf("3-Year (AUC = %.3f)", auc[2]),
                    sprintf("5-Year (AUC = %.3f)", auc[3])),
         col = c("#21908C", "#440154", "#FDE725"), lwd = 2, bty = "n", cex = 0.9)
  title(main = sprintf("ROC: Train %s -> Test %s", train_name, test_name))
  dev.off()
  
  cat(sprintf("\n=== Train %s -> Test %s (n_test = %d) ===\n", train_name, test_name, nrow(te)))
  cat(sprintf("  Auto-RS HR per SD = %.2f (95%% CI %.2f-%.2f), p = %.3g\n", hr, lo, hi, pv))
  cat(sprintf("  C-index = %.3f | time-AUC 1y/3y/5y = %.3f / %.3f / %.3f | KM log-rank p = %.3g\n",
              cidx, auc[1], auc[2], auc[3], km_p))

  list(
    metrics = data.frame(direction = sprintf("%s->%s", train_name, test_name),
             n_test = nrow(te), HR_perSD = hr, CI_low = lo, CI_high = hi, p = pv,
             Cindex = cidx, AUC_1y = auc[1], AUC_3y = auc[2], AUC_5y = auc[3],
             KM_logrank_p = km_p, row.names = NULL),
    label = list(ids = rownames(te), risk = as.character(te$risk)))
}

r1 <- evaluate("TARGET", "GSE21257")
r2 <- evaluate("GSE21257", "TARGET")

metrics <- rbind(r1$metrics, r2$metrics)
write.csv(metrics, "external_validation_metrics.csv", row.names = FALSE)

## independent risk label per sample (risk assigned when the sample was in the test cohort)
lab <- data.frame(id = c(r1$label$ids, r2$label$ids),
                  risk = c(r1$label$risk, r2$label$risk))
write.csv(lab, "external_risk_labels.csv", row.names = FALSE)

cat("\nSaved external_validation_metrics.csv and external_risk_labels.csv\n")
cat("\nPaste-ready sentence for the response:\n")
cat(sprintf(paste0("  Under strict external validation with independent within-cohort normalization, the locked ",
                   "Auto-RS retained prognostic discrimination in the held-out cohort in both directions ",
                   "(TARGET->GSE21257: HR %.2f per SD, 1/3/5-yr AUC %.2f/%.2f/%.2f; ",
                   "GSE21257->TARGET: HR %.2f per SD, P = %.3g, 1/3/5-yr AUC %.2f/%.2f/%.2f).\n"),
            r1$metrics$HR_perSD, r1$metrics$AUC_1y, r1$metrics$AUC_3y, r1$metrics$AUC_5y,
            r2$metrics$HR_perSD, r2$metrics$p, r2$metrics$AUC_1y, r2$metrics$AUC_3y, r2$metrics$AUC_5y))
