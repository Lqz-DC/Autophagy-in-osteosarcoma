################################################################################
## drug_ranking_plot.R  -  waterfall ranking of the nominated drugs
##
## Reuses the `avg` object built in drug_independent_and_anchor.R, part (B):
##   avg has columns: drug, mean_LN_IC50  (measured GDSC mean over the 4 OS lines)
## Run AFTER that block, or rebuild avg with the two lines marked below.
################################################################################

suppressMessages({ library(ggplot2); library(ggrepel); library(dplyr) })

## --- rebuild `avg` if running standalone (otherwise reuse the one in memory) ---
avg <- read.csv("The average drug sensitivity of the four cell lines of osteosarcoma.csv.csv",
                check.names = FALSE)
colnames(avg) <- c("drug","mean_LN_IC50")
avg <- avg[is.na(suppressWarnings(as.numeric(avg$drug))), ]

avg <- avg %>%
  arrange(mean_LN_IC50) %>%
  mutate(rank = row_number(),
         pct  = rank / n() * 100)

hi <- c("Gemcitabine","Cytarabine","Obatoclax Mesylate","AZD4547","Staurosporine")
pal <- c("Gemcitabine"="#C0392B","Cytarabine"="#C0392B",     # FDA agents, high-risk
         "Obatoclax Mesylate"="#2C7FB8","AZD4547"="#2C7FB8",  # other high-risk-leaning
         "Staurosporine"="#8E44AD")                            # low-risk preferential

ord <- function(n) { n <- round(n)
  suf <- ifelse(n %% 100 %in% 11:13, "th",
         ifelse(n %% 10 == 1, "st", ifelse(n %% 10 == 2, "nd",
         ifelse(n %% 10 == 3, "rd", "th")))); paste0(n, suf) }

sub <- avg %>% filter(tolower(drug) %in% tolower(hi)) %>%
  mutate(col = pal[drug],
         lab = sprintf("%s\n%s pct  (lnIC50 %.2f)", drug, ord(pct), mean_LN_IC50))

med <- median(avg$mean_LN_IC50)

p <- ggplot(avg, aes(rank, mean_LN_IC50)) +
  geom_area(fill = "#e9edf2") +
  geom_line(color = "#9aa5b1", linewidth = 0.5) +
  geom_point(color = "#c5ccd6", size = 0.7) +
  geom_hline(yintercept = med, linetype = "dashed", color = "#94a3b8", linewidth = 0.4) +
  annotate("text", x = nrow(avg), y = med + 0.2, hjust = 1, size = 3,
           color = "#64748b", label = sprintf("median (lnIC50=%.1f)", med)) +
  geom_point(data = sub, aes(color = drug), size = 3.2, show.legend = FALSE) +
  geom_label_repel(data = sub, aes(label = lab, color = drug),
                   size = 3, label.size = 0.4, box.padding = 0.8,
                   min.segment.length = 0, seed = 1, show.legend = FALSE) +
  scale_color_manual(values = pal) +
  labs(x = "Drug rank by mean sensitivity in 4 OS cell lines  (1 = most sensitive, of 257)",
       y = "Mean LN(IC50)  across HOS, MG63, U2OS, SAOS2",
       title = "Measured GDSC sensitivity ranking of nominated agents") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12))

ggsave("drug_sensitivity_ranking.pdf", p, width = 9.2, height = 4.8)
ggsave("drug_sensitivity_ranking.png", p, width = 9.2, height = 4.8, dpi = 300)
cat("saved drug_sensitivity_ranking.{pdf,png}\n")
