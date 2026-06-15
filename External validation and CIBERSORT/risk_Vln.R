library(tidyverse)
library(dplyr)
results <- read.table("CIBERSORT.filter.txt",header = T,sep = "\t",quote = "",row.names = 1)
risk_info <- read.table("risk.txt",header = T,sep = "\t",quote = "")
#data_process
data_long<-results%>%mutate(Sample=rownames(results))%>%gather(key="Cell_Type",value="Abundance",-Sample)
data_long<-merge(data_long,risk_info,by.x='Sample',by.y='id',all.x=T)

p_values<-data_long%>%group_by(Cell_Type)%>%summarise(p_value=wilcox.test(Abundance~risk)$p.value)

data_long<-data_long%>%left_join(p_values,by="Cell_Type")%>%
                mutate(Significance=case_when(
                      p_value<0.001~"***",
                      p_value<0.01~"**",
                      p_value<0.05~"*",
                      TRUE~"ns"))
# 1. 专门创建一个用于标注显著性的独立数据集
annotation_data <- data_long %>%
  group_by(Cell_Type) %>%
  summarise(
    Max_Abundance = max(Abundance), # 获取每种细胞各自的最大值
    Significance = first(Significance) # 每种细胞只保留一个显著性标记
  ) %>%
  ungroup()
# 2. 开始绘图
plot <- ggplot(data_long, aes(x=Cell_Type, y=Abundance, fill=risk, color=risk)) +
  geom_boxplot(fill=NA, outlier.shape=NA, size=1.2) +
  labs(y="Abundance", x="", title="Immune Cell Abundance", color='') +
  theme_bw() +
  scale_fill_manual(values=c("high" = '#EBBA37', "low" = "#8AD293")) +
  scale_color_manual(values=c("high" = '#EBBA37', "low" = "#8AD293")) +
  
  # 【核心修改】指定独立的数据集，不再使用全局的 max()
  geom_text(data = annotation_data, 
            aes(x = Cell_Type, y = Max_Abundance * 1.05, label = Significance), 
            size = 6, color = 'black', show.legend = FALSE, inherit.aes = FALSE) +
  
  theme(plot.title=element_text(hjust=0.5, size=30, face='bold'),
        legend.text=element_text(size=20),
        plot.margin=unit(c(3,3,3,3),'cm'),
        axis.ticks=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, size=20, color='black'),
        axis.text.y=element_text(size=20, color='black'))

plot
ggsave(plot=plot, filename='risk_Vln2.pdf', width=20, height=8.5)


stats_table <- data_long %>%
  group_by(Cell_Type) %>%
  summarise(
    n_high     = sum(risk == "high"),
    n_low      = sum(risk == "low"),
    mean_high  = mean(Abundance[risk == "high"]),
    mean_low   = mean(Abundance[risk == "low"]),
    median_high= median(Abundance[risk == "high"]),
    median_low = median(Abundance[risk == "low"]),
    p_value    = wilcox.test(Abundance ~ risk)$p.value,
    .groups    = "drop"
  ) %>%
  mutate(
    direction = ifelse(mean_high > mean_low, "high > low", "low > high"),
    fold      = pmax(mean_high, mean_low) / pmin(mean_high, mean_low),
    p_BH      = p.adjust(p_value, method = "BH"),     # 多重校正
    sig_raw   = case_when(p_value < .001 ~ "***",
                          p_value < .01  ~ "**",
                          p_value < .05  ~ "*",
                          TRUE           ~ "ns"),
    sig_BH    = case_when(p_BH < .001 ~ "***",
                          p_BH < .01  ~ "**",
                          p_BH < .05  ~ "*",
                          TRUE        ~ "ns")
  ) %>%
  arrange(p_value)

print(stats_table, n = Inf, width = Inf)

# Only look at those corrected by BH
print(stats_table %>% filter(p_BH < 0.05))

# save
write.csv(stats_table, "immune_celltype_stats.csv", row.names = FALSE)
