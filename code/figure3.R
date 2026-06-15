library(Seurat)
library(monocle3)
library(harmony)
library(CytoTRACE2)
library(patchwork) 
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(DoubletFinder)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSEABase)
library(enrichplot)
library(tidyverse)
library(ComplexHeatmap)
library(SeuratWrappers)
set.seed(1)

seurat_standard_normalize_and_scale <- function(colon, cluster, cluster_resolution){
  # colon is seurat object, 
  colon <- NormalizeData(colon, normalization.method = "LogNormalize", scale.factor = 10000)
  colon <- FindVariableFeatures(colon, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(colon)
  colon <- ScaleData(colon, features = all.genes)
  colon <- RunPCA(colon, features = VariableFeatures(object = colon))
  if (cluster){
    colon <- FindNeighbors(colon, dims = 1:20)
    colon <- FindClusters(colon, resolution = cluster_resolution)
  }
  colon <- RunUMAP(colon, dims = 1:20)
  return(colon)
}
make_seurat_object_and_doublet_removal <- function(data_directory, project_name){
  # function for basic seurat based qc and doubletfinder based doublet removal
  
  setwd("F:\\scRNA\\GSE162454_RAW")
  colon.data <- Read10X(data.dir = data_directory)
  currentSample <- CreateSeuratObject(counts = colon.data, project = project_name, min.cells = 3, min.features = 40)
  currentSample[["mt_percent"]] <- PercentageFeatureSet(currentSample, pattern = "^MT-")
  
  # qc plot-pre filtering
  setwd("F:\\scRNA\\GSE162454_RAW")
  pdf(paste0("./qc_plots_", project_name, "_prefiltered.pdf"))
  print(VlnPlot(currentSample, features = c("nFeature_RNA", "nCount_RNA", "mt_percent"), ncol = 3, pt.size = 0.05))
  dev.off()
  pdf(paste0("./qc_plots_", project_name, "_prefiltered_no_points.pdf"))
  print(VlnPlot(currentSample, features = c("nFeature_RNA", "nCount_RNA", "mt_percent"), ncol = 3, pt.size = 0))
  dev.off()
  
  # filter everything to 400 unique genes/cell
  currentSample <- subset(currentSample, subset =  nFeature_RNA > 200 & 
                            nFeature_RNA < 6000 & 
                            mt_percent < 10 & 
                            nCount_RNA <30000& 
                            nCount_RNA > 1000)
  
  # Normalize and make UMAP
  currentSample <- seurat_standard_normalize_and_scale(currentSample, FALSE)
  
  # Run doublet finder
  nExp_poi <- round(0.075 * ncol(currentSample))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  seu_colon <- doubletFinder(currentSample, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  print(head(seu_colon@meta.data))
  
  # rename columns
  seu_colon$doublet.class <- seu_colon[[paste0("DF.classifications_0.25_0.09_",nExp_poi)]]
  seu_colon[[paste0("DF.classifications_0.25_0.09_",nExp_poi)]] <- NULL
  pann <- grep(pattern="^pANN", x=names(seu_colon@meta.data), value=TRUE)
  seu_colon$pANN <- seu_colon[[pann]]
  seu_colon[[pann]] <- NULL
  
  # plot pre and post doublet finder results
  pdf(paste0("./UMAP_pre_double_removal", project_name, ".pdf"))
  print(DimPlot(seu_colon, reduction = "umap", group.by = "doublet.class", cols = c("#D51F26", "#272E6A")))
  dev.off()
  seu_colon <- subset(seu_colon, subset = doublet.class != "Doublet")
  pdf(paste0("./UMAP_post_double_removal", project_name, ".pdf"))
  print(DimPlot(seu_colon, reduction = "umap", cols = c("#D51F26")))
  dev.off()
  
  # Remove extra stuff and return filtered Seurat object
  seu_colon <- DietSeurat(seu_colon, counts=TRUE, data=TRUE, scale.data=FALSE, assays="RNA")
  return(seu_colon)
}

seurat_qc_plots <- function(colon, sample_name){
  # Make some basic qc plots
  pdf(paste0("./seurat_nFeature_plots_", sample_name, ".pdf"), width = 40, height = 15)
  print(VlnPlot(colon, features = c("nFeature_RNA"), ncol = 1, pt.size = 0.2))
  dev.off()
  
  pdf(paste0("./seurat_nCount_plots_", sample_name, ".pdf"), width = 40, height = 15)
  print(VlnPlot(colon, features = c("nCount_RNA"), ncol = 1, pt.size = 0.2))
  dev.off()
  
  pdf(paste0("./seurat_pMT_plots_", sample_name, ".pdf"), width = 40, height = 15)
  print(VlnPlot(colon, features = c("mt_percent"), ncol = 1, pt.size = 0.2))
  dev.off()
}


data_directory=c("GSE162454_OS_1_/","GSE162454_OS_2_/","GSE162454_OS_3_/","GSE162454_OS_4_/","GSE162454_OS_5_/","GSE162454_OS_6_/")
project_name=c("sample1","sample2","sample3","sample4","sample5","sample6")


samples <- project_name

sample1 <- make_seurat_object_and_doublet_removal(data_directory[1], samples[1])


### Merge multiple samples
seu_list <- sample1
for (i in 2:length(samples)){
  
  
  
  sc.i = make_seurat_object_and_doublet_removal(data_directory[i], samples[i])
  seu_list=merge(seu_list,sc.i)
  
}

table(seu_list$orig.ident)


scRNA_harmony=seu_list
scRNA_harmony  <- NormalizeData(scRNA_harmony ) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)

scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "orig.ident")


ElbowPlot(scRNA_harmony,ndims = 50)

scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution =seq(from=0.1,to=1,by=0.1))

Idents(scRNA_harmony) <- "RNA_snn_res.0.2"

scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:20)
scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = 1:20)
DimPlot(scRNA_harmony , reduction = "umap",label = T) 
DimPlot(scRNA_harmony , reduction = "tsne",label = T) 
DimPlot(scRNA_harmony, reduction = "umap", split.by ='orig.ident')

DimPlot(scRNA_harmony, reduction = "umap", group.by='orig.ident')
table(scRNA_harmony$orig.ident)  

scRNA_harmony=JoinLayers(scRNA_harmony)


#Create a collection of markers based on the literature
marker <- c('ALPL','RUNX2', 'IBSP',#Osteoblastic OS cells 
            'LYZ','CD68',"CD163",#Myeloid cells
            #"ACAN","COL2A1","SOX9",#chondroblasts
            'ACP5', 'CTSK',# Osteoclastic cells 
            'FBLN1','ACTA2','TAGLN','COL3A1', 'COL6A1',#Fibroblasts (CAFs)
            'CD2', 'CD3D','CD3E', 'CD3G','GNLY', 'NKG7', 'KLRD1', 'KLRB1',#NK/T cells
            'EGFL7', 'PLVAP',#endothelial cells
            'MS4A1', 'CD79A',#B cells 
            'IGHG1', 'MZB1'#plasma cells
)

par(mar = c(5, 4, 8, 3))
DotPlot(scRNA_harmony,features = marker)+RotatedAxis()+  scale_color_gradientn(colors = c(
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

DimPlot(scRNA_harmony , reduction = "umap",label = T) 


scRNA_harmony$celltype <- recode(scRNA_harmony@meta.data$seurat_clusters,
                                 "0"="Myeloid cells","1"="NK/T cells","2"="Osteoblastic OS cells",
                                 "3"="Myeloid cells","4"="Osteoblastic OS cells","5"="Myeloid cells",
                                 "6"="Fibroblasts (CAFs)","7"="Plasma cells","8"="Osteoclastic cells",
                                 "9"="B cells","10"="Endothelial cells")

DimPlot(scRNA_harmony , reduction = "umap",group.by = "celltype",label = T) 
table(scRNA_harmony$celltype)
custom_colors <- c(
  "Osteoblastic OS cells" = "#54A849",
  "Myeloid cells" = "#BC86C4",
  "Osteoclastic cells" = "#EB6D65",
  "Fibroblasts (CAFs)" = "#E572AF",
  "NK/T cells" = "#7084CB",
  "Endothelial cells" = "#69B5E8",
  "B cells" = "#BF811D",
  "Plasma cells" = "#57B094"
)

#Extract the expressed data
gene_expr <- FetchData(scRNA_harmony, vars = c("BNIP3", "MYC", "SAR1A", "PEA15"))

#Calculate the weighted risk score
scRNA_harmony$RiskScore <- with(gene_expr,
                                0.4775 * BNIP3 +
                                  0.6158 * MYC +
                                  0.4827 * SAR1A -
                                  0.6128 * PEA15)
risk_df <- data.frame(
  cluster = scRNA_harmony$celltype,  
  RiskScore = scRNA_harmony$RiskScore
)

#Sort the cluster order (by mean risk from low to high)
cluster_order <- risk_df %>%
  group_by(cluster) %>%
  summarise(mean_risk = mean(RiskScore)) %>%
  arrange(mean_risk) %>%
  pull(cluster)
risk_df$cluster <- factor(risk_df$cluster, levels = cluster_order)

write.csv(risk_cluster_df,"risk_cluster_df.csv",row.names = F,quote = F)

#Visualize the risk score
VlnPlot(scRNA_harmony, features = "RiskScore", group.by = "celltype", pt.size = 0.1) +
  ggtitle("Cox Model-based Prognostic Risk Score per Cluster")


#Calculate the average RiskScore of each cluster
avg_risk_by_cluster <- scRNA_harmony@meta.data %>%
  group_by(cluster = celltype) %>%
  summarise(mean_risk = mean(RiskScore, na.rm = TRUE),
            sd_risk = sd(RiskScore, na.rm = TRUE),
            n = n(),
            se = sd_risk / sqrt(n),
            lower = mean_risk - 1.96 * se,
            upper = mean_risk + 1.96 * se)

head(avg_risk_by_cluster)

ggplot(risk_df, aes(x = cluster, y = RiskScore, fill = cluster)) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "red") +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2, fill = "red") +
  scale_fill_manual(values = cluster_colors) +
  labs(x = NULL, y = "Risk Score") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )


#Specified color
custom_colors <- c(
  "Osteoblastic OS cells" = "#54A849",
  "Myeloid cells"         = "#BC86C4",
  "Osteoclastic cells"    = "#EB6D65",
  "Fibroblasts (CAFs)"    = "#57B094",
  "NK/T cells"            = "#7084CB",
  "Endothelial cells"     = "#69B5E8",
  "B cells"               = "#BF811D",
  "Plasma cells"          = "#E572AF"
)

DimPlot(scRNA_harmony,reduction = "umap",label = T,group.by = "celltype",cols = custom_colors,pt.size = 0.6,alpha = 1,stroke= 0)
ggsave("celltype.pdf",width = 12,height = 6)

DotPlot(scRNA_harmony,features = marker,group.by = "celltype")+
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
    "#E31A1C"))+  
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.8))
ggsave("Total marker gene.pdf",width = 10,height = 5.6)

save(scRNA_harmony,file = "scRNA_harmony_celltype.Rdata")
load("scRNA_harmony_celltype.Rdata")

#Cell stack bar chart####
#Calculate the quantity of each cell type in different groups
cellnum <- table(scRNA_harmony$celltype, scRNA_harmony$orig.ident)
cell.prop <- as.data.frame(prop.table(cellnum))

#Set the column name
colnames(cell.prop) <- c("Celltype", "Group", "Proportion")

#Draw a stacked bar chart
p.bar <- ggplot(cell.prop, aes(x = Group, y = Proportion, fill = Celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = custom_colors) + 
  ggtitle("Cell Type Proportion by Group") +
  theme_bw() +
  theme(axis.ticks.length = unit(0.5, 'cm')) +
  guides(fill = guide_legend(title = NULL))

print(p.bar)
ggsave("Cell proportion bar chart.pdf",width = 4.5,height = 5.6)

####sub Osteoblastic OS cells####
os_sc <- subset(scRNA_harmony,ident='Osteoblastic OS cells')
os_sc <- subset(os_sc, subset = CD3D < 0.1 & CD3E < 0.1 & CD3G < 0.1 & CD2 < 0.1)
#creat seurat Object
os_sc <- CreateSeuratObject(counts = GetAssayData(os_sc,assay = 'RNA',layer = 'counts'),meta.data = os_sc@meta.data)

#Complete the preprocessing：
os_sc <- NormalizeData(os_sc)
os_sc <- FindVariableFeatures(os_sc)
os_sc <- ScaleData(os_sc)
os_sc <- RunPCA(os_sc)

ElbowPlot(os_sc, ndims = 50)
DimPlot(os_sc, group.by = "orig.ident", reduction = "pca")  # Before going to batch
os_sc <- RunHarmony(os_sc, group.by.vars = "orig.ident")
# Replace PCA with Harmony for result display
DimPlot(os_sc, reduction = "harmony", group.by = "orig.ident") #After batch removal

os_sc <- RunUMAP(os_sc, reduction = "harmony", dims = 1:20)
os_sc <- RunTSNE(os_sc, reduction = "harmony", dims = 1:20)
os_sc <- FindNeighbors(os_sc, reduction = "harmony", dims = 1:20)

os_sc <- FindClusters(os_sc, resolution =seq(from=0.01,to=0.1,by=0.01))


Idents(os_sc) <- "RNA_snn_res.0.1"
colors <- c(
  "0"="#8AC760",
  "1"= "#899FD1",
  "2"= "#12737C",
  "3"= "#BF6DAB",
  "4"= "#F47D2B",
  "5"= "#3DBCA9"
)
os_sc$seurat_clusters <-  Idents(os_sc)
os_sc$celltype <- recode(os_sc@meta.data$seurat_clusters,
                         "HighRiskCluster"="ROS0/2/3",
                         "LowRiskCluster"="ROS5/1/4")
Idents(os_sc) <- "celltype"
os_sc@meta.data[["SingleR_Main"]] <- NULL

DimPlot(os_sc, reduction = "umap", group.by = "seurat_clusters",label = T,pt.size = 0.3,cols = colors,alpha = 1,stroke=0)
ggsave("os subcluster umap.pdf",width = 12,height = 8)

DimPlot(os_sc, reduction = "tsne", group.by = "seurat_clusters",label = T,pt.size = 1.5)
VlnPlot(os_sc,features = c( "MYC", "BNIP3", "SAR1A", "PEA15","MKI67","FBLN2"),group.by = "cluster",pt.size = 0.06,ncol = 3,alpha = 0.1,cols =colors)

ggsave("6_markers vinplot.pdf",width = 22,height = 10)
#Search for each cluster markergene
markers <- FindAllMarkers(object = os_sc, test.use="wilcox",
                          only.pos =TRUE,
                          logfc.threshold=0.25)

write.csv(markers,"marker.csv")

#top10 <-  markers  %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top5 <- markers %>%
  filter((pct.1 - pct.2) > 0.25) %>%  
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)      


sc.s <- subset(os_sc, downsample = min(table(Idents(os_sc))))
table(os_sc$RNA_snn_res.0.1)
DoHeatmap(sc.s,features = top5$gene,label =T,slot = "scale.data",group.colors =colors )+scale_fill_gradientn(
  colors = colorRampPalette(c("#410a5e", "#589088", "#eaf24e"))(100)
)
ggsave("top 10 of Os marker gene pheatmap.pdf",width = 8,height = 8)

####Construct a linear weighted prognostic risk score based on the risk scoring formula constructed by the Cox model####
gene_expr <- FetchData(os_sc, vars = c("BNIP3", "MYC", "SAR1A", "PEA15"))

#Calculate the weighted risk score
os_sc$RiskScore <- with(gene_expr,
                        0.4775 * BNIP3 +
                          0.6158 * MYC +
                          0.4827 * SAR1A -
                          0.6128 * PEA15)

risk_cluster_os <- data.frame(
  Celltype =os_sc$cluster,
  RiskScore = os_sc$RiskScore
)
risk_cluster_os$Celltype <- factor(risk_cluster_os$Celltype, 
                                   levels = c("5", "1", "4", "0", "2", "3"))
write.csv(risk_cluster_os,"risk_cluster_os.csv",row.names = F,quote = F)


VlnPlot(os_sc, features = "RiskScore", group.by = "seurat_clusters", pt.size = 0.1) +
  ggtitle("Cox Model-based Prognostic Risk Score per Cluster")
#Classify high and low risks based on the median
os_sc$risk_group <- ifelse(os_sc$RiskScore > median(os_sc$RiskScore), "High", "Low")
#Define high-risk clusters
high_risk_clusters <- c("0", "2","3")

#Add a risk group label to each cell
os_sc$cluster_risk_group <- ifelse(os_sc$seurat_clusters %in% high_risk_clusters, 
                                   "ROS0/2/3", 
                                   "ROS5/1/4")
# Visualize the distribution on UMAP
DimPlot(os_sc, group.by = "cluster_risk_group", label = TRUE, reduction = "umap",pt.size = 0.1) +
  ggtitle("ROS0/2/3 vs ROS5/1/4")


#Make sure the cluster information is in the meta.data
os_sc$cluster <- Idents(os_sc)

#Calculate the average RiskScore of each cluster
avg_risk_by_cluster <- os_sc@meta.data %>%
  group_by(cluster = seurat_clusters) %>%
  summarise(mean_risk = mean(RiskScore, na.rm = TRUE),
            sd_risk = sd(RiskScore, na.rm = TRUE),
            n = n(),
            se = sd_risk / sqrt(n),
            lower = mean_risk - 1.96 * se,
            upper = mean_risk + 1.96 * se)

print(avg_risk_by_cluster)
#Violin painting display
ggplot(risk_cluster_os, aes(x = Celltype, y = RiskScore, fill = Celltype)) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "red") +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2, fill = "red") +
  scale_fill_manual(values = colors) +
  labs(x = NULL, y = "Risk Score") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

save(os_sc,file = "os_sc.Rdata")

Idents(os_sc) <- os_sc$cluster_risk_group



####Analysis of differences in OS subgroups####
Idents(os_sc) <- os_sc$cluster_risk_group

#Search for marker genes of high-risk clusters
markers <- FindMarkers(os_sc, ident.1 = "ROS0/2/3", ident.2 = "ROS5/1/4",
                       logfc.threshold = 0.25, min.pct = 0.25)

highrisk_markers <- markers[markers$avg_log2FC > 0 & markers$p_val_adj < 0.05, ]
lowrisk_markers <- markers[markers$avg_log2FC < 0 & markers$p_val_adj < 0.05, ]

head(markers[markers$avg_log2FC > 0, ])

#Extract significantly upregulated genes
up_genes <- rownames(markers)[markers$p_val_adj < 0.05 & markers$avg_log2FC > 0.25]
down_genes <- rownames(markers)[markers$p_val_adj < 0.05 & markers$avg_log2FC < -0.25]
# SYMBOL → ENTREZ ID
gene_entrez <- bitr(down_genes, fromType = "SYMBOL", #down_genes
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)

go_result <- enrichGO(gene         = gene_entrez$ENTREZID,
                      OrgDb        = org.Hs.eg.db,
                      keyType      = "ENTREZID",
                      ont          = "ALL",
                      pAdjustMethod= "BH",
                      qvalueCutoff = 0.05)
#Split the geneID field in the enrichGO result and convert it to Symbol
go_df <- as.data.frame(go_result)

#Define a function to convert the Entrez ID string to the Gene Symbol string
convert_entrez_to_symbol <- function(entrez_string) {
  entrez_ids <- unlist(strsplit(entrez_string, "/"))
  symbols <- bitr(entrez_ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)$SYMBOL
  paste(unique(symbols), collapse = "/")
}

#Application function transformation
go_df$GeneSymbols <- sapply(go_df$geneID, convert_entrez_to_symbol)

head(go_df[, c("ID", "Description", "geneID", "GeneSymbols")])

pdf("GO BP Enrichment in ROS0/2/3.pdf",width =6,height =6.5)
pdf("GO BP Enrichment in ROS5/1/4.pdf",width =6,height =6.5)
dotplot(go_result, showCategory = 20, title = "GO BP Enrichment in ROS0/2/3",color="qvalue",font.size = 10)
dotplot(go_result, showCategory = 20, title = "GO BP Enrichment in ROS5/1/4",color="qvalue",font.size = 10)
dev.off()

kegg_result <- enrichKEGG(gene = gene_entrez$ENTREZID,
                          organism = 'hsa',
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05)

#Conversion function: Convert the Entrez ID string to the Gene Symbol string
convert_entrez_to_symbol <- function(entrez_string) {
  entrez_ids <- unlist(strsplit(entrez_string, "/"))
  symbols <- bitr(entrez_ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)$SYMBOL
  paste(unique(symbols), collapse = "/")
}


kegg_df <- as.data.frame(kegg_result)

# Add a column of gene symbols
kegg_df$GeneSymbols <- sapply(kegg_df$geneID, convert_entrez_to_symbol)


head(kegg_df[, c("ID", "Description", "geneID", "GeneSymbols")])

pdf("KEGG Enrichment in ROS0/2/3.pdf",width =6,height =6.5)
pdf("KEGG Enrichment in ROS5/1/4.pdf",width =6,height =6.5)
dotplot(kegg_result, showCategory = 20,color="qvalue", title = "KEGG Enrichment in ROS0/2/3",font.size = 10)
dotplot(kegg_result, showCategory = 20,color="qvalue", title = "KEGG Enrichment in ROS5/1/4",font.size = 10)
dev.off()

####GSEA####
Idents(os_sc) <- os_sc$cluster_risk_group
#find markergene#
markers <- FindMarkers(os_sc, ident.1 = "HighRiskCluster", ident.2 = "LowRiskCluster",
                       logfc.threshold = 0.25, min.pct = 0.25)

write.csv(markers,"risk_cluster_marker.csv")
markers$gene <- rownames(markers)

gene_list <- markers %>%
  arrange(desc(avg_log2FC)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  pull(avg_log2FC)
names(gene_list) <- markers %>%
  arrange(desc(avg_log2FC)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  pull(gene)


head(gene_list)
gmt_file <- "h.all.v2025.1.Hs.symbols.gmt" #Hallmark
gene_sets <- getGmt(gmt_file)

# Extract the terms and genes of gene sets
geneset <- data.frame(
  term = rep(sapply(gene_sets, function(x) x@setName), 
             times = sapply(gene_sets, function(x) length(geneIds(x)))),
  gene = unlist(sapply(gene_sets, geneIds))
)

# Remove the "HALLMARK_"
geneset$term <- gsub(pattern ="HALLMARK_","", geneset$term)
geneset$term <- gsub(pattern ="KEGG_","", geneset$term)

egmt <- GSEA(gene_list,
             TERM2GENE = geneset,
             pvalueCutoff =1, # Set it to 1 to obtain all results
             minGSSize =1,
             maxGSSize =500000)

egmt <- egmt@result

data<-egmt[,c("ID","NES","setSize","pvalue")]
data$Gene_Number<-data$setSize
head(data)

#Set the path names on the Y-axis to factors and sort them
data<-data[order(data$NES,decreasing=T),]
data$ID<-factor(data$ID,levels=data$ID)
data$xlab<-1:50
head(data)
summary(data$NES)
summary(data$Gene_Number)

#The names of the channels marked in the figure
label<-c( 
  "HYPOXIA",
  "EPITHELIAL_MESENCHYMAL_TRANSITION",
  "ANGIOGENESIS",
  
  "MYC_TARGETS_V1",
  "MTORC1_SIGNALING",
  "PEROXISOME",
  "OXIDATIVE_PHOSPHORYLATION",
  "TNFA_SIGNALING_VIA_NFKB",
  "IL2_STAT5_SIGNALING",
  "KRAS_SIGNALING_UP",
  "IL6_JAK_STAT3_SIGNALING",
  
  "E2F_TARGETS",
  "G2M_CHECKPOINT",
  "MITOTIC_SPINDLE"
)

#Extract the corresponding path
data_label<-data[data$ID%in%label,]
data_label
#The color of the passage in the picture
data_label$col<-c("purple","purple","purple","grey","grey","grey","grey","grey","grey","grey","grey","#EBBA37","#EBBA37","#EBBA37")

p<-ggplot(data=data,aes(x=xlab,y=NES))+
  geom_point(aes(size=Gene_Number,alpha=-log10(pvalue)),shape=21,stroke=0,fill="#0000ff",colour="black")+
  scale_size_continuous(range=c(4,12))+
  xlab(label="Hallmark gene sets")+
  ylab(label="Normalized enrichment score (NES)")+
  theme_classic(base_size=15)+
  scale_x_continuous(breaks=seq(0,50,by=10),labels=seq(0,50,by=10))+
  scale_y_continuous(breaks=seq(-3,3,by=1),labels=seq(-3,3,by=1))+
  guides(size=guide_legend(title="Gene Number"),
         alpha=guide_legend(title="-log10(pvalue)"))+
  theme(
    axis.line=element_line(color="black",size=0.6),
    axis.text=element_text(face="bold"),
    axis.title=element_text(size=13)
  )

p

p1<-p+
  geom_text_repel(data=data_label,aes(x=xlab,y=NES,label=ID),size=3,color=data_label$col,
                  force=20,
                  point.padding=0.5,
                  min.segment.length=0,
                  hjust=1.2,
                  segment.color="grey20",
                  segment.size=0.3,
                  segment.alpha=1,
                  nudge_y=-0.2,
                  nudge_x = 0.1
  )
p1

####cytotrace2####
#Get the original counts
expr_matrix <- as.matrix(GetAssayData(os_sc,assay = 'RNA',slot = 'counts'))
# Run CytoTRACE
results <- cytotrace2(expr_matrix, species = "human")

# Add the score to meta.data
os_sc$CytoTRACE_score <- cyto_out$CytoTRACE[colnames(os_sc)]

os_sc <- AddMetaData(os_sc, metadata = results)

# Convert to a data frame
head(os_sc$seurat_clusters)
os_sc$seurat_clusters <- Idents(os_sc)
annotation_df <- data.frame(cluster = os_sc$seurat_clusters)
rownames(annotation_df) <- names(os_sc$seurat_clusters)

plots <- plotData(
  cytotrace2_result = results,
  expression_data = expr_matrix,
  annotation = annotation_df,   
  is_seurat = FALSE
)

# Add CytoTRACE2 results
os_sc$CytoTRACE2_Score <- results$CytoTRACE2_Score
os_sc$CytoTRACE2_Potency <- results$CytoTRACE2_Potency
os_sc$CytoTRACE2_Relative <- results$CytoTRACE2_Relative

head(os_sc@meta.data[, c("CytoTRACE2_Score", "CytoTRACE2_Potency", "CytoTRACE2_Relative")])

FeaturePlot(os_sc, features = "CytoTRACE2_Score", reduction = "umap") +
  scale_color_viridis_c(option = "magma") +
  ggtitle("CytoTRACE2 Score on UMAP")

DimPlot(os_sc, group.by = "CytoTRACE2_Potency", reduction = "umap", label = TRUE) +
  ggtitle("CytoTRACE2 Potency Category on UMAP")

# Extract UMAP and CytoTRACE2 information
umap_df <- Embeddings(os_sc, "umap") %>% as.data.frame()
umap_df$CytoTRACE2 <- os_sc$CytoTRACE2_Score  
umap_df$pseudotime <- umap_df$CytoTRACE2

ggplot(umap_df, aes(x = umap_1, y = umap_2, color = pseudotime)) +
  geom_point(size =3, alpha = 0.2,stroke=0) +
  scale_color_gradientn(colors = c("#03A750","#F7B90F", "#B91A2E")) +
  ggtitle("CytoTRACE2 Score on UMAP") +
  theme_void()
ggsave("cytoTRACE2.pdf",width = 10,height = 5)


print(plots$CytoTRACE2_UMAP)           
print(plots$CytoTRACE2_Potency_UMAP)   
print(plots$CytoTRACE2_Relative_UMAP) 


cds$cyto_score <- os_sc$CytoTRACE2_Score[colnames(cds)]
plot_cells(cds, color_cells_by = "cyto_score")


####monocle3####
cds <- readRDS("OS_cell_monocle3_result.rds")
p1 <- plot_cells(cds,
                 color_cells_by = "pseudotime",
                 label_branch_points = FALSE,
                 label_leaves = FALSE,
                 cell_stroke = 0,
                 cell_size = 2.5,alpha = 0.3,
                 trajectory_graph_color = "black",
                 label_roots = T) + 
  scale_color_gradientn(colors = c("#03A750","#F7B90F", "#B91A2E")) +
  theme_void() 
p1
ggsave("monocle3 trajectory.pdf",width = 10,height = 5)

###CNV####
###Endothelial cells are used as reference cells####
table(scRNA_harmony$celltype)

#Extract Endothelial cells
endo_subset <- subset(scRNA_harmony, ident = "Endothelial cells")

# Make sure Endothelial cells are not repeated in os_sc
endo_subset <- subset(endo_subset, cells = setdiff(colnames(endo_subset), colnames(os_sc)))

# Merge into os_sc
os_sc <- merge(os_sc, endo_subset)
os_sc=JoinLayers(os_sc)
# Create CNV_group (original CNV grouping)
os_sc$CNV_group <- ifelse(os_sc$seurat_clusters %in% c("1", "4", "5"), "ROS5/1/4",
                          ifelse(os_sc$seurat_clusters %in% c("0", "2", "3"), "ROS0/2/3", NA))

#Label Endothelial cells with reference_OB
endo_in_os <- colnames(endo_subset)
os_sc$CNV_group[endo_in_os] <- "reference_Endothelial"


table(os_sc$CNV_group)

# Construct the annotation file
annotation_df <- data.frame(
  cell_id = colnames(os_sc),
  group   = os_sc$CNV_group
)

# Save
write.table(annotation_df, "annotation.txt", sep = "\t",  row.names = FALSE, col.names = FALSE, quote = FALSE)
#Expression matrix (original counts are recommended)
expr_mat <- GetAssayData(os_sc, assay = "RNA", slot = "counts")
write.table(as.matrix(expr_mat), "os_expr_counts.txt", sep = "\t", quote = FALSE)

#Gene location information file (gene annotation)
library(biomaRt)
mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
gene_info <- getBM(attributes = c("external_gene_name", "chromosome_name", "start_position", "end_position"),
                   filters = "external_gene_name",
                   values = rownames(expr_mat),
                   mart = mart)


#Permitted chromosome names
valid_chr <- as.character(c(1:22))
#rename：Gene, Chr, Start, End
colnames(gene_info) <- c("Gene","Chr","Start","End")


# Filter out non-standard chromosome rows
gene_info <- gene_info[gene_info$Chr %in% valid_chr, ]
# Convert Chr to factor and set levels to the character form from 1 to 22
gene_info$Chr <- factor(gene_info$Chr, levels = as.character(1:22))
#Sort the gene_info data frames by Chr
gene_info <- gene_info[order(gene_info$Chr), ]
table(gene_info$Chr)

write.table(gene_info, file = "gene_pos.txt", sep = "\t", row.names =F, col.names = F, quote = FALSE)

#Run the main function of inferCNV
library(infercnv)

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = "os_expr_counts.txt",
  annotations_file =  "annotation.txt",
  delim = "\t",
  gene_order_file = "gene_pos.txt",   # The location information generated in the previous step
  ref_group_names = c("reference_Endothelial") # The name of the set reference group
)

# Start running the CNV inference analysis
infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff =0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics 
                              out_dir = "infercnv_output_endo2",
                              cluster_by_groups = TRUE,
                              denoise = TRUE,
                              output_format = "pdf",
                              HMM = F,num_threads = 1)  

