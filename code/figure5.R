library(impute)
library(pheatmap)
library(tidyverse)
library(ggplot2)
library(ggridges)
library(dplyr)
library(tidyr)
library(rstatix)
library(ggpubr)
#Read the drug prediction results (one drug per row, one cell line per column)
df1 <- read.csv("drugPred.csv", sep = ',')

#Convert to long format
df1_long <- pivot_longer(df1, cols = -DRUGS, names_to = "CellLine", values_to = "Prediction")

#Convert back to the wide format, using "DRUGS" as the line name and "CellLine" as the column name
Pred_wide <- df1_long %>%
  pivot_wider(names_from = CellLine, values_from = Prediction)

#Set the line name and delete the "DRUGS" column
Pred_wide <- as.data.frame(Pred_wide)
rownames(Pred_wide) <- Pred_wide$DRUGS

Predictions <- as.data.frame(Pred_wide[, -1])

# Load the mean and standard deviation table, with the row names being the drug names
sdmean <- read.csv("Drugs_mean_sd.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE, row.names = 1)

#Combine the mean/standard deviation and the prediction results
pred <- merge(sdmean, Predictions, by = "row.names")
rownames(pred) <- pred$Row.names

#Implement z-score standardization
Predictions <- sweep(pred[, 5:ncol(pred)], 1, pred$mean, "-") / pred$sd
rownames(Predictions) <- pred$Row.names

#Read the sample risk information
risk_info <- read.table("sample_risk.txt", header = TRUE, row.names = 1)


# Ensure the sequence is consistent
risk_info <- risk_info[colnames(mat), , drop = FALSE]  

# Set as a column comment
labels <- data.frame(Risk = as.factor(risk_info$risk))
rownames(labels) <- rownames(risk_info)


ann_colors <- list(
  Risk = c("high" = "#ebba37", "low" = "#00d293"))


#Create a color vector
my_colors <- c(
  "#000000",    
  "#00008B",
  "#5e60ce",    
  "#87CEFA",    
  "#FFFFFF",    
  "#FFFACD",    
  "#FFA07A",    
  "#F08080"    
)

my_palette <- colorRampPalette(my_colors)

gradient_colors <- my_palette(100)

res <- pheatmap(
  Predictions,
  fontsize_col = 8,
  col = gradient_colors,
  angle_col = 45,
  fontsize_row = 10,
  cluster_cols = F,
  cluster_rows = T,
  annotation_col = labels,
  annotation_colors = ann_colors,
  annotation_names_row = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  show_rownames = FALSE,
  show_colnames = F,
  fontsize = 8)
print(res)
dev.off()

#ridge plot
sdmean <- tibble::rownames_to_column(sdmean, var = "DRUGS")

pred <- df1 %>%
  left_join(sdmean, by = "DRUGS")


pred_cols <- setdiff(colnames(df1), "DRUGS")


pred_z <- pred %>%
  mutate(across(all_of(pred_cols), ~ (. - mean) / sd))


df_long <- pred_z %>%
  select(DRUGS, all_of(pred_cols)) %>%
  pivot_longer(cols = -DRUGS, names_to = "Sample", values_to = "zscore") %>%
  rename(Drug = DRUGS)

# Read the risk grouping information and ensure that the Sample column names are consistent
risk_df <- read.table("sample_risk.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  rename(Sample = id, Risk = risk)


df <- df_long %>%
  left_join(risk_df, by = "Sample")

selected_drugs <- c("Cytarabine","Staurosporine","Obatoclax Mesylate","Vinblastine",
                    "Sepantronium bromide",
                    "AZD4547","Gemcitabine",
                    "Vincristine",
                    "Vinorelbine"
)
df_sub <- df[df$Drug %in% selected_drugs, ]

pdf('Inter-group difference ridge map.pdf',width = 6,height = 5)
ggplot(df_sub, aes(x = zscore, y = Drug, fill = Risk)) +
  geom_density_ridges(alpha = 0.8, scale = 1.1, rel_min_height = 0.01,bandwidth = 0.1) +
  theme_minimal(base_size = 16) +
  scale_fill_manual(values = c("high" = "#E41A1C", "low" = "#377EB8")) +
  xlab("z-score IC50") + ylab("Drug") +
  theme(axis.text.y = element_text(size = 10))
dev.off()

#Statistical analysis
# Calculate the median
medians <- df_sub %>%
  group_by(Drug, Risk) %>%
  summarize(median_zscore = median(zscore, na.rm = TRUE), .groups = "drop")

print(medians)


stat_test <- df_sub %>%
  group_by(Drug) %>%
  wilcox_test(zscore ~ Risk) %>%      
  adjust_pvalue(method = "BH") %>%    
  add_significance("p.adj")

print(stat_test)

####correlation####

df_long <- read.csv("drugPred.csv") %>%
  pivot_longer(cols = -DRUGS, names_to = "Sample", values_to = "Predicted") %>%
  rename(Drug = DRUGS)


risk_df <- read.table("sample_risk.txt", header = TRUE, sep = "\t")
colnames(risk_df) <- c("Sample", "Risk")

df_long <- left_join(df_long, risk_df, by = "Sample")


real_ic50 <- read.csv("The average drug sensitivity of the four cell lines of osteosarcoma.csv", sep = ",")
colnames(real_ic50) <- c("Drug", "Real")

df_plot <- left_join(df_long, real_ic50, by = "Drug") %>%
  filter(!is.na(Real) & !is.na(Predicted) & !is.na(Risk))

#Calculate the average predicted value for each group
df_avg_grouped <- df_plot %>%
  group_by(Drug, Risk) %>%
  summarise(
    Predicted_mean = mean(Predicted, na.rm = TRUE),
    Real = first(Real), .groups = "drop"
  )

pdf('correlation.pdf',width = 8,height = 8)
ggscatter(df_avg_grouped,
          x = "Real", y = "Predicted_mean",
          color = "Risk", shape = "Risk",
          palette = c("high" = "#cd9820", "low" = "#4fcd8d"),
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "pearson", cor.coef.size = 5,
          size = 3, alpha = 0.7) +
  stat_cor(aes(color = Risk),
           label.x.npc = "left", label.y.npc = c("top", "top"),
           size = 4, show.legend = FALSE) +
  labs(title = "Correlation Between Mean Predicted and Observed LN IC50",
       subtitle = "Stratified by Risk Group (Each dot = one drug)",
       x = "Observed LN IC50 (mean of 4 osteosarcoma cell lines)",
       y = "Predicted LN IC50 (mean across group)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )
dev.off()


####IDWAS####
analyze_gene_drug <- function(gene_index, expr_matrix, drug_prediction) {
  if (gene_index > nrow(expr_matrix)) stop("gene_index out of bounds")
  gene_expr <- expr_matrix[gene_index, ]
  results <- list()
  
  for(drug in colnames(drug_prediction)) {
    if (length(drug_prediction[, drug]) != length(gene_expr)) {
      stop("Length mismatch between drug_prediction and gene_expr")
    }
    valid <- !is.na(gene_expr) & !is.na(drug_prediction[, drug])
    if (sum(valid) < 3) {
      results[[drug]] <- c(beta = NA, pvalue = NA)
      next
    }
    model <- try(lm(drug_prediction[valid, drug] ~ gene_expr[valid]), silent = TRUE)
    if (inherits(model, "try-error")) {
      results[[drug]] <- c(beta = NA, pvalue = NA)
      next
    }
    coef_summary <- coef(summary(model))
    if (nrow(coef_summary) < 2) {
      results[[drug]] <- c(beta = NA, pvalue = NA)
    } else {
      results[[drug]] <- c(
        beta = coef_summary[2, 1],
        pvalue = coef_summary[2, 4]
      )
    }
  }
  return(results)
}


expr_idwas <- function(drug_prediction, expr_matrix, n = 10, output_file = "gene_drug_associations.csv") {
  #incomplete parameter checking
  if (!is.data.frame(drug_prediction)) stop("drug_prediction must be data.frame")
  if (!is.matrix(expr_matrix)) stop("expr_matrix must be matrix")
  
  #Sample intersection
  common_samples <- intersect(rownames(drug_prediction), colnames(expr_matrix))
  drug_prediction <- drug_prediction[common_samples, , drop = FALSE]
  expr_matrix <- expr_matrix[, common_samples, drop = FALSE]
  
  # Filter out genes with insufficient expression samples
  valid_genes <- apply(expr_matrix, 1, function(x) sum(!is.na(x))) >= n
  expr_matrix <- expr_matrix[valid_genes, , drop = FALSE]
  
  #paralyzer
  library(parallel)
  cl <- makeCluster(detectCores())
  clusterExport(cl, varlist = c("expr_matrix", "drug_prediction", "analyze_gene_drug"), envir = environment())
  all_results <- parLapply(cl, 1:nrow(expr_matrix), function(i) {
    analyze_gene_drug(i, expr_matrix, drug_prediction)
  })
  stopCluster(cl)
  
  #output data.frame
  result_df <- do.call(rbind, lapply(1:length(all_results), function(i) {
    gene <- rownames(expr_matrix)[i]
    data.frame(
      gene = gene,
      drug = names(all_results[[i]]),
      beta = sapply(all_results[[i]], `[`, "beta"),
      pvalue = sapply(all_results[[i]], `[`, "pvalue"),
      stringsAsFactors = FALSE
    )
  }))
  
  #Multiple inspection correction
  result_df$fdr <- p.adjust(result_df$pvalue, method = "BH")
  
  #Derive significant results
  significant <- result_df[result_df$fdr < 0.05, ]
  write.csv(significant, file = output_file, row.names = FALSE)
  
  return(significant)
}
#input
drug_prediction <- read.csv("drugPred.csv", row.names = 1)  # Sample × Drug
drug_prediction <- t(drug_prediction)
drug_prediction <- as.data.frame(drug_prediction)
range(drug_prediction)
expr_matrix <- read.table("exp.txt", sep = "\t", header = TRUE, row.names = 1)  # Gene × Sample


# running analysis

results <- expr_idwas(
  drug_prediction = drug_prediction,
  expr_matrix = as.matrix(expr_matrix),
  n = 10,
  output_file = "gene_drug_associations.csv"
)

results$drug <- sub("_\\d+$", "", results2$drug)
results$drug <- gsub("\\.","-",results2$drug)
write.csv(results,"gene_drug_association.csv",row.names = F)
####Frequency distribution histogram####

results <- results %>% filter(fdr < 0.05) %>% filter(beta<0)
# Set the target drugs and target genes
target_drug <- "Staurosporine"#AZD4547,Cytarabine,Gemcitabine,Obatoclax Mesylate
target_gene <-c("LY86","PEA15","SLAMF8")# c("MYC","PEA15","BIN2")#PEA15,SAR1A,BNIP3,MYC

drug_subset <- subset(results, drug == target_drug)

target_value <- -log10(drug_subset$pvalue[drug_subset$gene %in% target_gene])


#Draw and add red lines for marking
pdf('Staurosporine.pdf',width = 4.5,height = 4.5)
ggplot(drug_subset, aes(x = -log10(pvalue))) +
  geom_histogram(binwidth = 1, fill = "lightblue", color = "black",size = 0.2) +
  geom_vline(xintercept = target_value, color = "red", linetype = "solid", size = 1) +
  annotate("text", x = target_value, y = Inf, label = target_gene,
           vjust = -0.5, hjust = 1.1, color = "red", angle = 90) +
  theme_minimal() +
  labs(
    title = paste("Drug", target_drug, "Significant distribution with all genes"),
    x = "-log10(pvalue)",
    y = "Frequency"
  )
dev.off()

#Scatter plot
target_drug <- "Staurosporine"
drug_subset <- subset(results, drug == target_drug)
top100 <- drug_subset %>%
  arrange(pvalue) %>%
  head(100)

# Find the rows corresponding to the smallest p-value and the smallest beta (with the largest negative value)
min_p_gene <- top100[which.min(top100$pvalue), ]
min_beta_gene <- top100[which.min(top100$beta), ]

pdf('Staurosporine_The top 100 scatter plots.pdf',width = 4.5,height = 4.5)
ggplot(top100, aes(x = -log10(pvalue), y = beta)) +
  geom_point(color = "lightblue", size =4) +
  # Mark the point with the smallest p value
  geom_text(data = min_p_gene, aes(label = gene), color = "red",
            vjust = -1, size =2) +
  
  # Mark the point with the smallest (most negative) beta_sensitive
  geom_text(data = min_beta_gene, aes(label = gene), color = "blue",
            vjust = 1.5, size = 2) +
  theme_minimal() +
  labs(
    title = paste("Drug", target_drug, "beta and significance"),
    x = "-log10(pvalue)",
    y = "Beta"
  )
dev.off()


