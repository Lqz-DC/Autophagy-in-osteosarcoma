library(Seurat)
library(monocle3)
library(harmony)
library(CellChat)
library(patchwork) 
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(DoubletFinder)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSEABase)
library(GSVA)
library(limma)
library(enrichplot)
library(tidyverse)
library(ComplexHeatmap)
library(nichenetr) 
library(tidyverse)
set.seed(1)


####sub CAFs cells####
CAFs <- subset(scRNA_harmony,ident='Fibroblasts (CAFs)')

#CreateSeuratObject
CAFs <- CreateSeuratObject(counts = GetAssayData(CAFs,assay = 'RNA',layer = 'counts'),meta.data = CAFs@meta.data)

#Complete the preprocessing
CAFs <- NormalizeData(CAFs)
CAFs <- FindVariableFeatures(CAFs)
CAFs <- ScaleData(CAFs)
CAFs <- RunPCA(CAFs)

ElbowPlot(CAFs, ndims = 50)
DimPlot(CAFs, group.by = "orig.ident", reduction = "pca")  
CAFs <- RunHarmony(CAFs, group.by.vars = "orig.ident")

DimPlot(CAFs, reduction = "harmony", group.by = "orig.ident")

CAFs <- RunUMAP(CAFs, reduction = "harmony", dims = 1:20)
CAFs <- RunTSNE(CAFs, reduction = "harmony", dims = 1:20)
CAFs <- FindNeighbors(CAFs, reduction = "harmony", dims = 1:20)

CAFs <- FindClusters(CAFs, resolution =seq(from=0.01,to=0.1,by=0.01))
CAFs <- FindClusters(CAFs, resolution =seq(from=0.1,to=0.1,by=1))


Idents(CAFs) <- "RNA_snn_res.0.1"


custom_colors <- c("myCAF"="#54A849",
                   "apCAF"="#BC86C4",
                   "vCAF"="#EB6D65",
                   "iCAF"="#69B5E8"
)


DimPlot(CAFs, reduction = "umap",group.by = "celltype",label = T,cols = custom_colors,pt.size = 5,alpha =0.5,stroke=0)

marker <- c("HLA-DRA", "CD74", "IGFBP3", "HLA-DPB1",#Antigen-presenting CAFs (apCAFs)
            "IL6",  "CD34", "PLA2G2A", "DPP4",  "C3", "PI16","CXCL14", "CXCL2", "IGF1", "ALDH1A1", "CCL2",#Inflammatory CAFs (iCAFs)
            "MCAM","NOTCH3", "COL18A1", "NR2F2", #Vessel-associated CAFs (vCAFs)
            "MMP11", "POSTN", "LRRC15", "FAP","VIM", "S100A4","LUM","LOXL1","COL6A3","FN1", "COL3A1", "SPON2", "COL5A1", "INHBA"#Matrix CAFs (mCAFs)
)

Idents(CAFs) <- "celltype"
DotPlot(CAFs, features = marker,group.by = "celltype",dot.scale = 7) +
  RotatedAxis() +
  scale_color_gradientn(colors = colors <- c(
    "#FDFEFF",
    "#B2E8EE",
    "#70D1DF",
    "#4CBCD0",
    "#6CCBBB",
    "#ADDAA7",
    "#FBEB96",
    "#F6C065", 
    "#F39A3C", 
    "#E31A1C" ))  
ggsave("CAF marker bubble plot.pdf",width = 10,height = 5.6)


CAFs$celltype <- recode(CAFs@meta.data$seurat_clusters,
                        "0"="mCAF","1"="apCAF","2"="mCAF",
                        "3"="apCAF","4"="apCAF","5"="iCAF")

#order
desired_order <- c("iCAF","apCAF","vCAF","myCAF")
CAFs$celltype <- factor(CAFs$celltype, levels = desired_order)

genes <- c("THBS1", "THBS2",  "CDH11", "PTX3", "ITGB2", "PTN", "CYP1B1", "FLNB", "CXCL12", "GAS6","MDK","TNC","IBSP","APP","HLA-F","CSF1","JAM3","MPZL1","CD99",# apCAF
           "CCL2", "CXCL2", "IGF1","HAS1","FBLN5", "THBS4", "VEGFB", "SCG2", "RARRES2","LIF","FGF7","MMP3","PDPN","CCL11","VEGFD","FGF18",#iCAF 
           "MMP9", "ANGPT2", "MCAM", "PTP4A3", "TGFB3","PDGFA", "CSPG4","MMP11","PGF","PTK2","NOTCH1","COL4A1","CD46","JAG1","SPP1","THY1",# vCAF 
           "BIRC5", "INHBA", "TGFBI", "SRPX2", "IL11", "WNT2", "PTTG1", "ITGA5", "FSTL3","PTHLH","TGFB1","MIF"#myCAF
)

custom_colors <- c("myCAF"="#54A849",
                   "apCAF"="#BC86C4",
                   "vCAF"="#EB6D65",
                   "iCAF"="#69B5E8"
)
CAFs$celltype <- factor(CAFs$celltype,levels = c("myCAF","vCAF","apCAF","iCAF"))

Idents(CAFs) <- "celltype"
DotPlot(CAFs, features = genes,group.by = "celltype",dot.scale = 7) +
  RotatedAxis() +
  #coord_flip()+
  ggtitle("DE Marker of CAFs")+
  scale_color_gradientn(colors = c(
    "#FDFEFF",
    "#B2E8EE",
    "#70D1DF",
    "#4CBCD0",
    "#6CCBBB",
    "#ADDAA7",
    "#FBEB96", 
    "#F6C065", 
    "#F39A3C", 
    "#E31A1C" ))  
ggsave("DE_CAFs_marke.pdf",width =15,height = 5)

#calculate RS
gene_expr <- FetchData(CAFs, vars = c("BNIP3", "MYC", "SAR1A", "PEA15"))

CAFs$RiskScore <- with(gene_expr,
                       0.4775 * BNIP3 +
                         0.6158 * MYC +
                         0.4827 * SAR1A -
                         0.6128 * PEA15)

risk_cluster_CAFs <- data.frame(
  Celltype =CAFs$celltype,
  RiskScore = CAFs$RiskScore
)

write.csv(risk_cluster_CAFs,"risk_cluster_CAFs.csv",row.names = F,quote = F)


VlnPlot(CAFs, features = "RiskScore", group.by = "Celltype", pt.size = 0.1) +
  ggtitle("CAFs Cox Model-based PrognCAFstic Risk Score per Cluster")

#Classify high and low risks based on the median
CAFs$risk_group <- ifelse(CAFs$RiskScore > median(CAFs$RiskScore), "High", "Low")


#Calculate the average RiskScore of each cluster
avg_risk_by_cluster <- CAFs@meta.data %>%
  group_by(cluster = Celltype) %>%
  summarise(mean_risk = mean(RiskScore, na.rm = TRUE),
            sd_risk = sd(RiskScore, na.rm = TRUE),
            n = n(),
            se = sd_risk / sqrt(n),
            lower = mean_risk - 1.96 * se,
            upper = mean_risk + 1.96 * se)

print(avg_risk_by_cluster)

#Violin picture display
desired_order <- c("apCAF","iCAF","vCAF","myCAF")
CAFs$celltype <- factor(CAFs$celltype, levels = desired_order)

Idents(CAFs) <- "celltype"

ggplot(risk_cluster_CAFs, aes(x = Celltype, y = RiskScore, fill = Celltype)) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "red") +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2, fill = "red") +
  scale_fill_manual(values = custom_colors) +
  labs(x = NULL, y = "Risk Score") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.pCAFsition = "none"
  )
ggsave("CAFs_riskscore.pdf",width = 6,height = 4.5)




##Draw a stacked bar chart####
cellnum <- table(CAFs$celltype, CAFs$orig.ident)
cell.prop <- as.data.frame(prop.table(cellnum))
table(CAFs$celltype)
colnames(cell.prop) <- c("Celltype", "Group", "Proportion")

p.bar <- ggplot(cell.prop, aes(x = Group, y = Proportion, fill = Celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = custom_colors) +
  ggtitle("Cell Type Proportion by Group") +
  theme_bw() +
  theme(axis.ticks.length = unit(0.5, 'cm')) +
  guides(fill = guide_legend(title = NULL))

print(p.bar)
ggsave("CAF Cell proportion bar chart.pdf",width = 4.5,height = 5.6)

####GSVA analysis####

CAFs$celltype <- factor(CAFs$celltype,levels = c("myCAF","vCAF","apCAF","iCAF"))
Idents(CAFs) <- CAFs$celltype

table(CAFs$celltype)

expr_mat <- GetAssayData(CAFs, layer = "data")  

cluster_info <- CAFs$celltype

#average expression per celltype
library(Matrix)
cluster_avg_expr <- sapply(levels(cluster_info), function(clust) {
  cells <- colnames(CAFs)[cluster_info == clust]
  rowMeans(expr_mat[, cells])
})

# Prepare the GMT file
gmt_file <- "c5.go.bp.v2025.1.Hs.symbols.gmt"

gene_sets <- getGmt(gmt_file)

####Run GSVA####
param <- gsvaParam(cluster_avg_expr, gene_sets, kcdf = "Gaussian")

gsva_res <- gsva(param)

rownames(gsva_res) <- gsub("^GOBP_", "", rownames(gsva_res))

# Z-score 
gsva_z <- t(scale(t(gsva_res)))
write.csv(gsva_z,"GSVA_c5BP_results.csv",row.names = T)

pdf("CAF_GSVA_gobp pathway.pdf",width =30,height =450)
head(gsva_z)
dim(gsva_z)
pheatmap(gsva_z,
         cluster_rows = TRUE,
         cluster_cols = F,
         color = colorRampPalette(c("#410a5e", "#589088", "#eaf24e"))(100),
         fontsize_row = 8,
         fontsize_col = 10,
         main = " GSVA Pathway Enrichment in CAFs Clusters")
dev.off()

####sub Myeloid_sc cells####

Myeloid_sc <- subset(scRNA_harmony,ident='Myeloid cells')

Myeloid_sc <- CreateSeuratObject(counts = GetAssayData(Myeloid_sc,assay = 'RNA',layer = 'counts'),meta.data = Myeloid_sc@meta.data)


Myeloid_sc <- NormalizeData(Myeloid_sc)
Myeloid_sc <- FindVariableFeatures(Myeloid_sc)
Myeloid_sc <- ScaleData(Myeloid_sc)
Myeloid_sc <- RunPCA(Myeloid_sc)

ElbowPlot(Myeloid_sc, ndims = 50)
DimPlot(Myeloid_sc, group.by = "orig.ident", reduction = "pca")  
Myeloid_sc <- RunHarmony(Myeloid_sc, group.by.vars = "orig.ident")

DimPlot(Myeloid_sc, reduction = "harmony", group.by = "orig.ident")

Myeloid_sc <- RunUMAP(Myeloid_sc, reduction = "harmony", dims = 1:20)
Myeloid_sc <- RunTSNE(Myeloid_sc, reduction = "harmony", dims = 1:20)
Myeloid_sc <- FindNeighbors(Myeloid_sc, reduction = "harmony", dims = 1:20)

Myeloid_sc <- FindClusters(Myeloid_sc, resolution =seq(from=0.01,to=0.1,by=0.01))
Myeloid_sc <- FindClusters(Myeloid_sc, resolution =seq(from=0.1,to=1,by=0.1))


Myeloid_sc$seurat_clusters <- Idents(Myeloid_sc)
table(Myeloid_sc$RNA_snn_res.0.3)
levels(Idents(Myeloid_sc))

DotPlot(Myeloid_sc, features = marker,dot.scale = 8) +
  RotatedAxis() +
  #coord_flip()+
  scale_color_gradientn(colors = c(
    "#FDFEFF",
    "#B2E8EE",
    "#70D1DF",
    "#4CBCD0",
    "#6CCBBB",
    "#ADDAA7",
    "#FBEB96", 
    "#F6C065", 
    "#F39A3C", 
    "#E31A1C"))  


custom_colors <- c(
  "0" = "#54A849",  
  "1" = "#BC86C4",  
  "2" = "#EB6D65",  
  "3" = "#E572AF",  
  "4" = "#7084CB",  
  "5" = "#69B5E8",  
  "6" = "#BF811D",  
  "7" = "#57B094",  
  "8" = "#B3D465",  
  "9" = "#F6C667"
)

Myeloid_sc$seurat_clusters <- Idents(Myeloid_sc)

Myeloid_sc$celltype <- recode(Myeloid_sc@meta.data$seurat_clusters,
                              "0"="TAMs","1"="TAMs","2"="TAMs",
                              "3"="Macrophages(IFNs response)","4"="Macrophages(Proliferating)","5"="Monocytes","6"="Macrophages(Stress-like)","7"="DCs","8"="DCs","9"="Mast cells")

custom_colors <- c(
  "TAMs" = "#efa6cd",  
  "Macrophages(IFNs response)" = "#bf811d",  
  "Macrophages(Proliferating)" = "#7084cb",  
  "Monocytes" = "#bc86c4",  
  "Macrophages(Stress-like)" = "#f5863a",  
  "DCs" = "#6ab5e8",  
  "Mast cells" = "#54a849"
)


DimPlot(Myeloid_sc, reduction = "umap",group.by = "seurat_clusters",label = T,pt.size = 2,cols = custom_colors,alpha =0.6,stroke=0)
ggsave("Myeloid_umap0.3.pdf",width = 10,height = 8)

Idents(Myeloid_sc) <- Myeloid_sc$celltype

DotPlot(Myeloid_sc, features=marker,dot.scale = 6) +
  RotatedAxis() +
  scale_color_gradientn(colors = c(
    "#FDFEFF",
    "#B2E8EE",
    "#70D1DF",
    "#4CBCD0",
    "#6CCBBB",
    "#ADDAA7",
    "#FBEB96",
    "#F6C065", 
    "#F39A3C", 
    "#E31A1C"))  
ggsave("Myeloid marker gene bubble plot.pdf",width = 16,height = 6)



gene_expr <- FetchData(Myeloid_sc, vars = c("BNIP3", "MYC", "SAR1A", "PEA15"))


# Cell stacking bar chart####
Idents(Myeloid_sc) <- "seurat_clusters"
cellnum <- table(Myeloid_sc$seurat_clusters, Myeloid_sc$orig.ident)
cell.prop <- as.data.frame(prop.table(cellnum))


colnames(cell.prop) <- c("Celltype", "Group", "Proportion")


p.bar <- ggplot(cell.prop, aes(x = Group, y = Proportion, fill = Celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = custom_colors) + 
  ggtitle("Cell Type Proportion by Group") +
  theme_bw() +
  theme(axis.ticks.length = unit(0.5, 'cm')) +
  guides(fill = guide_legend(title = NULL))

print(p.bar)
ggsave("Bar chart of myeloid cell proportion.pdf",width = 8,height = 9.84)


# calculate RS####
Myeloid_sc$RiskScore <- with(gene_expr,
                             0.4775 * BNIP3 +
                               0.6158 * MYC +
                               0.4827 * SAR1A -
                               0.6128 * PEA15)

risk_cluster_Myeloid_sc <- data.frame(
  Celltype =Myeloid_sc$celltype,
  RiskScore = Myeloid_sc$RiskScore
)
#Classify high and low risks based on the median
Myeloid_sc$risk_group <- ifelse(Myeloid_sc$RiskScore > median(Myeloid_sc$RiskScore), "High", "Low")


#Calculate the average RiskScore of each cluster
avg_risk_by_cluster <- Myeloid_sc@meta.data %>%
  group_by(cluster = celltype) %>%
  summarise(mean_risk = mean(RiskScore, na.rm = TRUE),
            sd_risk = sd(RiskScore, na.rm = TRUE),
            n = n(),
            se = sd_risk / sqrt(n),
            lower = mean_risk - 1.96 * se,
            upper = mean_risk + 1.96 * se)

print(avg_risk_by_cluster)

desired_order <- c("DCs","Monocytes","Macrophages(IFNs response)","Macrophages(Proliferating)","Macrophages(Stress-like)","TAMs","Mast cell")
risk_cluster_Myeloid_sc$Celltype <- factor(Myeloid_sc$celltype, levels = desired_order)


Idents(Myeloid_sc) <- "celltype"

ggplot(risk_cluster_Myeloid_sc, aes(x = Celltype, y = RiskScore, fill = Celltype)) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "red") +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2, fill = "red") +
  scale_fill_manual(values = custom_colors) +
  labs(x = NULL, y = "Risk Score") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
ggsave("Myeloid_sc_riskscore.pdf",width = 6,height = 4.5)


####cellchat####
load("combined.Rdata")

#Check the risk scores of each group
gene_expr <- FetchData(combined, vars = c("BNIP3", "MYC", "SAR1A", "PEA15"))

#Calculate the weighted risk score
combined$RiskScore <- with(gene_expr,
                           0.4775 * BNIP3 +
                             0.6158 * MYC +
                             0.4827 * SAR1A -
                             0.6128 * PEA15)
table(combined$celltype)
risk_cluster_combined <- data.frame(
  Celltype =combined$celltype,
  RiskScore = combined$RiskScore
)
# Classify high and low risks based on the median
combined$risk_group <- ifelse(combined$RiskScore > median(combined$RiskScore), "High", "Low")

#Calculate the average RiskScore of each cluster
avg_risk_by_cluster <- combined@meta.data %>%
  group_by(cluster = combined$celltype) %>%
  summarise(mean_risk = mean(RiskScore, na.rm = TRUE),
            sd_risk = sd(RiskScore, na.rm = TRUE),
            n = n(),
            se = sd_risk / sqrt(n),
            lower = mean_risk - 1.96 * se,
            upper = mean_risk + 1.96 * se)

print(avg_risk_by_cluster)


avg_risk_by_cluster <- avg_risk_by_cluster %>%
  arrange(desc(mean_risk))

print(avg_risk_by_cluster)

# Extract the expression matrix and meta information
data.input <- GetAssayData(combined, slot = "data")
meta <- combined@meta.data

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
cellchat@DB <- CellChatDB.human  
#Run the CellChat analysis
cellchat <- subsetData(cellchat)  # Extract receptor-ligand expression
future::plan("multisession", workers = 1)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat,population.size = FALSE )
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

####The heatmap displays the number of interactions in order of risk score.####
order_vec <- c(
  "ROS0/2/3",
  "myCAF",
  "vCAF",
  "iCAF",
  "ROS5/1/4",
  "apCAF",
  "Mast cells",
  "NK cells",
  "TAMs",
  "T cells",
  "Macrophages(Stress-like)",
  "Macrophages(Proliferating)",
  "Macrophages(IFNs response)",
  "Monocytes",
  "DCs"
)
mat <- cellchat@net$count

# Rearrange the rows and columns in order
mat <- mat[order_vec, order_vec]
head(mat)
mat <- mat[, rev(colnames(mat))]

colors_vec <- c(
  "#FFFFFF", "#FFF9DD", "#FFF3BB", "#FFEB99", "#FFE176",
  "#FFD054", "#FEB24C", "#FD9842", "#FC7E38", "#F9612E",
  "#F24524", "#E31A1C", "#C4121D", "#A00D1E", "#800026"
)
pdf("interaction_number2.pdf", width = 6.7, height = 5)
pheatmap::pheatmap(mat, border_color = "black", color = colors_vec,#colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                   cluster_cols = F, fontsize = 10, cluster_rows = F,
                   display_numbers = T,number_color="black",angle_col = 45,number_format = "%.0f")
dev.off()


#Show all the intercellular ligand-receptor interactions
p <- netVisual_bubble(cellchat, sources.use = 1:length(levels(cellchat@idents)), 
                      targets.use = 1:length(levels(cellchat@idents)), remove.isolate = FALSE)
ggsave("CCI_all.pdf", p, width =20, height = 25, limitsize = F)


pdf("interaction_strength.pdf",width =7,height = 5)
netVisual_heatmap(cellchat, measure = "weight",color.heatmap =c( "#FFFFFF",  "#A020F0"))
dev.off()


#Nichnet

organism <- "human"

lr_network <- readRDS("lr_network_human_21122021.rds")
ligand_target_matrix <- readRDS("ligand_target_matrix_nsga2r_final.rds")
weighted_networks <- readRDS("weighted_networks_nsga2r_final.rds")


lr_network <- lr_network %>% distinct(from, to)
head(lr_network)
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
head(weighted_networks$gr) # interactions and their weights in the gene regulatory network

table(combined$celltype)
####1. Define a set of potential ligands for sender agnosticism and sender focusing methods####
receiver = "ROS0/2/3"
expressed_genes_receiver <- get_expressed_genes(receiver, combined, pct = 0.05)

all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)

potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()

sender_celltypes <- c( "ROS0/2/3",
                       "myCAF",
                       "vCAF",
                       "iCAF",
                       "ROS5/1/4",
                       "apCAF",
                       "Mast cells",
                       "NK cells",
                       "TAMs",
                       "T cells",
                       "Macrophages(Stress-like)",
                       "Macrophages(Proliferating)",
                       "Macrophages(IFNs response)",
                       "Monocytes",
                       "DCs")
# lapply is used to obtain the expressed genes of each sending cell type respectively
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, combined, 0.05)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()
potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 

length(expressed_genes_sender)
length(potential_ligands)
length(potential_ligands_focused)
####2. Define the set of genes of interest####
condition_oi <-  "ROS0/2/3"
condition_reference <- "ROS5/1/4"

seurat_obj_receiver <- combined

DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
                                  ident.1 = condition_oi, ident.2 = condition_reference,
                                  group.by = "celltype",
                                  min.pct = 0.05) %>% rownames_to_column("gene")

geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

####3.Define background gene####
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

length(background_expressed_genes)
##[1] 10604
length(geneset_oi)
## [1] 3310

####4. Carry out the ligand activity analysis of NicheNet####
ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
ligand_activities

p_hist_lig_activity <- ggplot(ligand_activities, aes(x=aupr_corrected)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(10, aupr_corrected) %>% pull(aupr_corrected))),
             color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()

p_hist_lig_activity

str(ligand_activities)

best_upstream_ligands <- ligand_activities %>% top_n(10, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)

vis_ligand_aupr <- ligand_activities %>%
  filter(test_ligand %in% best_upstream_ligands) %>%
  arrange(aupr_corrected) %>%
  column_to_rownames(var = "test_ligand") %>%
  dplyr::select(aupr_corrected) %>%  
  as.matrix(ncol=1)

(make_heatmap_ggplot(vis_ligand_aupr,
                     "Prioritized ligands", "Ligand activity", 
                     legend_title = "AUPR", color = "darkorange") + 
    theme(axis.text.x.top = element_blank()))  

####5. Infer the target genes and receptors of the top ligand####
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()

nrow(active_ligand_target_links_df)

head(active_ligand_target_links_df)

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33) 
nrow(active_ligand_target_links)
## [1] 232
head(active_ligand_target_links)
order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                    color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")


ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

(make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                     y_name = "Ligands", x_name = "Receptors",  
                     color = "mediumvioletred", legend_title = "Prior interaction potential"))
####6. Infer the target genes and receptors of the top ligand####

ligand_activities_all <- ligand_activities 
best_upstream_ligands_all <- best_upstream_ligands

ligand_activities <- ligand_activities %>% filter(test_ligand %in% potential_ligands_focused)
best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>%
  pull(test_ligand) %>% unique()

ligand_aupr_matrix <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>%  dplyr::select(aupr_corrected)  %>% arrange(aupr_corrected)
vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1) 

p_ligand_aupr <- make_heatmap_ggplot(vis_ligand_aupr,
                                     "Prioritized ligands", "Ligand activity", 
                                     legend_title = "AUPR", color = "darkorange") + 
  theme(axis.text.x.top = element_blank())

p_ligand_aupr

# Target gene plot
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33) 

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

p_ligand_target <- make_heatmap_ggplot(vis_ligand_target[,1:50], "Prioritized ligands", "Predicted target genes",
                                       color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")

p_ligand_target

# Receptor plot
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

p_ligand_receptor <- make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                                         y_name = "Ligands", x_name = "Receptors",  
                                         color = "mediumvioletred", legend_title = "Prior interaction potential")

p_ligand_receptor
best_upstream_ligands_all %in% rownames(combined) %>% table()

#order
desired_order <- c( "ROS0/2/3",
                    "myCAF",
                    "vCAF",
                    "iCAF",
                    "ROS5/1/4",
                    "apCAF",
                    "Mast cells",
                    "NK cells",
                    "TAMs",
                    "T cells",
                    "Macrophages(Stress-like)",
                    "Macrophages(Proliferating)",
                    "Macrophages(IFNs response)",
                    "Monocytes",
                    "DCs")
combined$celltype <- factor(combined$celltype, levels = desired_order)
levels(combined$celltype)
Idents(combined) <- combined$celltype
# Dotplot of sender-focused approach
p_dotplot <- DotPlot(subset(combined, celltype %in% sender_celltypes),
                     features = rev(best_upstream_ligands),dot.scale = 4)+
  scale_color_gradientn(colors =c(
    "#FDFEFF",
    "#B2E8EE",
    "#70D1DF",
    "#4CBCD0",
    "#6CCBBB",
    "#ADDAA7",
    "#FBEB96", # 黄
    "#F6C065", # 黄橙
    "#F39A3C", # 橙
    "#E31A1C" )) + 
  coord_flip() +
  scale_y_discrete(position = "right")
p_dotplot


(make_line_plot(ligand_activities = ligand_activities_all,
                potential_ligands = potential_ligands_focused) +
    theme(plot.title = element_text(size=11, hjust=0.1, margin=margin(0, 0, -5, 0))))

####7. Visual summary of the NicheNet analysis####

figures_without_legend <- cowplot::plot_grid(
  p_ligand_aupr + theme(legend.position = "none"),
  p_dotplot + theme(legend.position = "none",
                    axis.ticks = element_blank(),
                    axis.title.y = element_blank(),
                    axis.title.x = element_text(size = 12),
                    axis.text.y = element_text(size = 9),
                    axis.text.x = element_text(size = 9,  angle = 90, hjust = 0)) +
    ylab("Expression in Sender"),
  p_ligand_target + theme(legend.position = "none",
                          axis.title.y = element_blank()),
  align = "hv",
  nrow = 1,
  rel_widths = c(1.5,1.5,4))

legends <- cowplot::plot_grid(
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_dotplot)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target)),
  nrow = 1,
  align = "h", rel_widths = c(3,1,1.5))

combined_plot <-  cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
combined_plot



#When the receiver is an immune cell####

table(combined$celltype)
receiver = "T cells"
expressed_genes_receiver <- get_expressed_genes(receiver, combined, pct = 0.05)

all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)

potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()

sender_celltypes <- c( "ROS0/2/3",
                       "myCAF",
                       "vCAF",
                       "iCAF",
                       "ROS5/1/4",
                       "apCAF",
                       "Mast cells",
                       "NK cells",
                       "TAMs",
                       "T cells",
                       "Macrophages(Stress-like)",
                       "Macrophages(Proliferating)",
                       "Macrophages(IFNs response)",
                       "Monocytes",
                       "DCs")
# lapply is used to obtain the expressed genes of each sending cell type respectively
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, combined, 0.05)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()
potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 

length(expressed_genes_sender)

length(potential_ligands)

length(potential_ligands_focused)

####2.lapply is used to obtain the expressed genes of each sending cell type respectively####
condition_oi <-  "High"
condition_reference <- "Low"

seurat_obj_receiver <- subset(combined, idents = receiver)

DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
                                  ident.1 = condition_oi, ident.2 = condition_reference,
                                  group.by = "risk_group",
                                  min.pct = 0.05) %>% rownames_to_column("gene")

geneset_oi <- DE_table_receiver %>% filter(p_val<= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

####3. Define background genes####
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

length(background_expressed_genes)

length(geneset_oi)


####4. Carry out the ligand activity analysis of NicheNet####
ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
ligand_activities

p_hist_lig_activity <- ggplot(ligand_activities, aes(x=aupr_corrected)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(30, aupr_corrected) %>% pull(aupr_corrected))),
             color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()

p_hist_lig_activity

str(ligand_activities)

best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)

vis_ligand_aupr <- ligand_activities %>%
  filter(test_ligand %in% best_upstream_ligands) %>%
  arrange(aupr_corrected) %>%
  column_to_rownames(var = "test_ligand") %>%
  dplyr::select(aupr_corrected) %>% 
  as.matrix(ncol=1)

(make_heatmap_ggplot(vis_ligand_aupr,
                     "Prioritized ligands", "Ligand activity", 
                     legend_title = "AUPR", color = "darkorange") + 
    theme(axis.text.x.top = element_blank()))  

####5. Carry out the ligand activity analysis of NicheNet####
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()

nrow(active_ligand_target_links_df)

head(active_ligand_target_links_df)

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33) 
nrow(active_ligand_target_links)

head(active_ligand_target_links)
order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

p_ligand_target <- make_heatmap_ggplot(vis_ligand_target[,1:50], "Prioritized ligands", "Predicted target genes",
                                       color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")
p_ligand_target

ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

(make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                     y_name = "Ligands", x_name = "Receptors",  
                     color = "mediumvioletred", legend_title = "Prior interaction potential"))
####6. Implement a send-centered approach####

ligand_activities_all <- ligand_activities 
best_upstream_ligands_all <- best_upstream_ligands

ligand_activities <- ligand_activities %>% filter(test_ligand %in% potential_ligands_focused)
best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>%
  pull(test_ligand) %>% unique()

ligand_aupr_matrix <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>%  dplyr::select(aupr_corrected)  %>% arrange(aupr_corrected)
vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1) 

p_ligand_aupr <- make_heatmap_ggplot(vis_ligand_aupr,
                                     "Prioritized ligands", "Ligand activity", 
                                     legend_title = "AUPR", color = "darkorange") + 
  theme(axis.text.x.top = element_blank())

p_ligand_aupr

# Target gene plot
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33) 

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

p_ligand_target <- make_heatmap_ggplot(vis_ligand_target[,1:50], "Prioritized ligands", "Predicted target genes",
                                       color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")

p_ligand_target

# Receptor plot
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

p_ligand_receptor <- make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                                         y_name = "Ligands", x_name = "Receptors",  
                                         color = "mediumvioletred", legend_title = "Prior interaction potential")

p_ligand_receptor
best_upstream_ligands_all %in% rownames(combined) %>% table()

#order
desired_order <- c( "ROS0/2/3",
                    "myCAF",
                    "vCAF",
                    "iCAF",
                    "ROS5/1/4",
                    "apCAF",
                    "Mast cells",
                    "NK cells",
                    "TAMs",
                    "T cells",
                    "Macrophages(Stress-like)",
                    "Macrophages(Proliferating)",
                    "Macrophages(IFNs response)",
                    "Monocytes",
                    "DCs")
combined$celltype <- factor(combined$celltype, levels = desired_order)
levels(combined$celltype)
Idents(combined) <- combined$celltype

p_dotplot <- DotPlot(subset(combined, celltype %in% sender_celltypes),
                     features = rev(best_upstream_ligands),dot.scale = 4)+
  scale_color_gradientn(colors =c(
    "#FDFEFF",
    "#B2E8EE",
    "#70D1DF",
    "#4CBCD0",
    "#6CCBBB",
    "#ADDAA7",
    "#FBEB96", # 黄
    "#F6C065", # 黄橙
    "#F39A3C", # 橙
    "#E31A1C" )) + 
  coord_flip() +
  scale_y_discrete(position = "right")
p_dotplot

(make_line_plot(ligand_activities = ligand_activities_all,
                potential_ligands = potential_ligands_focused) +
    theme(plot.title = element_text(size=11, hjust=0.1, margin=margin(0, 0, -5, 0))))

####7.Visual summary of the NicheNet analysis####

figures_without_legend <- cowplot::plot_grid(
  p_ligand_aupr + theme(legend.position = "none"),
  p_dotplot + theme(legend.position = "none",
                    axis.ticks = element_blank(),
                    axis.title.y = element_blank(),
                    axis.title.x = element_text(size = 12),
                    axis.text.y = element_text(size = 9),
                    axis.text.x = element_text(size = 9,  angle = 90, hjust = 0)) +
    ylab("Expression in Sender"),
  p_ligand_target + theme(legend.position = "none",
                          axis.title.y = element_blank()),
  align = "hv",
  nrow = 1,
  rel_widths = c(1.5,1.5,4))

legends <- cowplot::plot_grid(
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_dotplot)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target)),
  nrow = 1,
  align = "h", rel_widths = c(3,1,1.5))

combined_plot <-  cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
combined_plot













