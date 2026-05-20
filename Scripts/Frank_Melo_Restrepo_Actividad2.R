


library(tximport)
library(readr)
library(DESeq2)
library(tidyverse)

#setwd("//wsl.localhost/Ubuntu-24.04/home/frank/Actividad5/Fastqs")
getwd()

dir()
samples <- c("AbrahamSimpson", "HomerSimpson", "MargeSimpson",
             "PattyBouvier", "SelmaBouvier")

files <- file.path("salida", samples, "quant.sf")
names(files) <- samples

tx2gene <- read.delim("Transcrito_a_Gen.tsv", header = F)

colnames(tx2gene) <- c("TXNAME", "GENEID")

tx2gene

####import salmon
txi <- tximport(files, type="salmon", tx2gene=tx2gene, dropInfReps = TRUE)
all(file.exists(files))
dir()


Design <- read.csv("Design.csv")

coldata <- Design %>% filter(Condition!="Normopeso")
coldata


coldata$Condition[coldata$Condition == "Sobrepeso/Obeso1"] <- "Obeso1"
coldata$Condition[coldata$Condition == "Sobrepeso/Obeso2"] <- "Obeso2"

coldata$Condition <- factor(coldata$Condition)
coldata$Sexo <- factor(coldata$Sexo)




str(coldata)
####dds
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~Condition)


dds <- DESeq(dds)
res <- results(dds)
View(dds)


nrow(dds)# numero de genes
names(dds)# nobre de los genes



res = results(dds, contrast=c( "Condition", "Obeso1", "Obeso2"), alpha=1e-5) 
res

####res


res_df= as.data.frame(res)


res_df <- data.frame(Gene.Name = rownames(res_df), res_df)

# res_df=merge(res_df, genes, by="row.names")
# res_df
# head(res_df)
# 
# genes_to_check = c("THY1", "SFMBT2", "PASD1", "SNAI1")
# res_df[res_df$Gene.Name %in% genes_to_check, ]
#6. Visualización de los resultados
# ------------------------------------------------------------------------------


# MA plot
plotMA(res) # Los valores azules son los genes que pasan el umbral que hemos especificado del alpha

# Volcano plot
BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
EnhancedVolcano(res_df, lab=res_df$Gene.Name, x='log2FoldChange', y='pvalue', labSize = 3, axisLabSize = 10)

### O siendo más estrictos
EnhancedVolcano(res_df, lab=res_df$Gene.Name, x='log2FoldChange', y='pvalue', pCutoff=1e-20, FCcutoff=5, labSize = 3, axisLabSize = 10)

### O prefiltrando genes no significativos para aligerar la nube de puntos:
res_filt <- res_df[
  res_df$pvalue < 0.001 & abs(res_df$log2FoldChange) > 2, ]
EnhancedVolcano(res_filt, lab=res_filt$Gene.Name, x='log2FoldChange', y='pvalue', pCutoff=1e-20, FCcutoff=5, labSize = 3, axisLabSize = 10)

# HeatMap-------------------



####-------

if(!require(pheatmap)){
  install.packages("pheatmap")
  library(pheatmap)
}
vsd <- varianceStabilizingTransformation(dds, blind = FALSE) # robusto por los pocos genes
# vsd <- vst(dds, blind = FALSE) ## Aplica la transformación vst() (variance-stabilizing transformation) de DESeq2 al objeto dds
mat <- assay(vsd)[(rownames(res)), ] ## Extrae la matriz de expresión transformada (filas = genes, columnas = muestras)
rownames(mat) <- res_df$Gene.Name ## Reemplaza los identificadores de fila por nombres de genes
mat_scaled <- t(scale(t(mat))) # Escala y transpone la matriz

pheatmap(mat_scaled,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = as.data.frame(colData(dds)),
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 3,
         fontsize_col = 8,
         color = colorRampPalette(c("blue", "white", "red"))(50))

# Extracción de lista top genes diferencialmente expresados

# Ordenar por significancia ajustada
res_ordered <- res_df[order(res_df$padj), ]
# Seleccionar los 10 genes más significativos
top10_genes <- head(res_ordered, 10)
top10_genes

mat_subset <- mat_scaled[(head(order(res$padj), 10)), ]

pheatmap(mat_subset,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = as.data.frame(colData(dds)),
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(50))
###plotpca

plotPCA(vsd, intgroup="Condition")



# Extraer datos del PCA
pca <- plotPCA(vsd, intgroup = "Condition", returnData = TRUE)
pca

# Calcular porcentaje de varianza explicada
percentVar <- round(100 * attr(pca, "percentVar"))



library(ggrepel)

ggplot(pca, aes(x = PC1, y = PC2, color = Condition, label = Sample)) +
  geom_point(size = 4) +
  geom_text_repel(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% varianza")) +
  ylab(paste0("PC2: ", percentVar[2], "% varianza")) +
  ggtitle("Análisis de Componentes Principales (PCA)") +
  theme_gray()



#7. Análisis de enriquecimiento----
# ------------------------------------------------------------------------------

pkgs_bioc <- c("clusterProfiler","enrichplot", "org.Hs.eg.db", "ReactomePA", "enrichplot")
for (p in pkgs_bioc) {
  if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p)
  library(p, character.only = TRUE)
}
BiocManager::install("enrichplot")
library("enrichplot")

## Creamos un conjunto de datos con símbolos únicos y los resultados de significancia
genes_df <- unique(data.frame(symbol = res_df$Gene.Name,
                              log2FC = res_df$log2FoldChange,
                              padj = res_df$padj,
                              stringsAsFactors = FALSE))
genes_df

## Eliminamos aquellos genes que no tienen nomenclatura, podríamos eliminar incluso con algún valor faltante
genes_df <- genes_df[!is.na(genes_df$symbol) & genes_df$symbol != "", ]

## Asociamos cada nombre de gen con el identificador de la base de datos ENTREZID
map <- bitr(genes_df$symbol,
            fromType = "SYMBOL",
            toType   = c("ENTREZID"),
            OrgDb    = org.Hs.eg.db)

map

## Unimos ambos datos
genes_mapped <- merge(genes_df, map, by.x = "symbol", by.y = "SYMBOL")
nrow(genes_df); nrow(genes_mapped) # Comprobamos duplicados: un nombre puede mapear a >1 ENTREZ (rara vez)
universe_entrez <- unique(map$ENTREZID) # Creamos el objeto de elementos únicos para los análisis

## Over-Representation Analysis (ORA)

sig_genes <- genes_mapped[genes_mapped$padj < 0.05 & abs(genes_mapped$log2FC) >= 1, ]
length(sig_genes$ENTREZID)
sig_genes





#######-----sirve GO------
res_df <- as.data.frame(res)
res_df
res_df <- na.omit(res_df)
genes_sig <- subset(res_df)

genes <- rownames(genes_sig)
genes

map <- bitr(genes, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)

map



ego <- enrichGO(
  gene          = map$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",   # BP = procesos biológicos
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.1,
  readable      = TRUE
)


head(ego)
dotplot(ego_bp)
barplot(ego, showCategory = 10)
ego


dotplot(ego, showCategory = 10) +
  ggtitle("GO Enrichment - Procesos Biológicos")


#enrichPathway-----


reactome <- enrichPathway(
  gene = map$ENTREZID,
  organism = "human",
  pvalueCutoff = 0.1,
  readable = TRUE
)


head(reactome)


dotplot(reactome, showCategory = 10) +
  ggtitle("Reactome Pathway Enrichment")

barplot(reactome, showCategory = 10)

#GSEA-----

gene_list <- res_df$log2FoldChange
names(gene_list) <- rownames(res_df)
gene_list

map <- bitr(names(gene_list),
            fromType = "SYMBOL",
            toType = "ENTREZID", 
            OrgDb = org.Hs.eg.db)



gene_df <- merge(map,
                 data.frame(SYMBOL= names(gene_list),
                            logFC= gene_list),
                 by.x="SYMBOL", by.Y="SYMBOL")
gene_df

gene_list_entrez <- gene_df$logFC

names(gene_list_entrez) <- gene_df$ENTREZID
gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)


gsea_go <- gseGO(
  geneList = gene_list_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pvalueCutoff = 0.1,
  verbose = FALSE
)

head(gsea_go)

dotplot(gsea_go, showCategory=10)+
  ggtitle("GSEA - GO Biological Processes")



dotplot(gsea_go, showCategory = 10) +
  ggtitle("GSEA GO Enrichment") +
  theme_minimal()



gseaplot2(gsea_go, geneSetID = 1)

ridgeplot(gsea_go, showCategory = 10)
gsea_go@result[1:5, ]



#Biomart----

library(biomaRt)

ensembl <- useEnsembl("genes", dataset = "hsapiens_gene_ensembl")

genes <- rownames(res_df)

annot <- getBM(
  attributes = c(
    "external_gene_name",
    "ensembl_gene_id",
    "chromosome_name",
    "start_position",
    "end_position",
    "gene_biotype"
  ),
  filters = "external_gene_name",
  values = genes,
  mart = ensembl
)

head(annot)
nrow(annot)

res_df$external_gene_name <- rownames(res_df)

merged_data <- merge(res_df, annot, by = "external_gene_name")
``
head(merged_data)

###hasta Aca---
