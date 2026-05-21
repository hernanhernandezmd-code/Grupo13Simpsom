############################################################
### Proyecto: Grupo 13 Simpson RNA-seq obesidad
### PARTE 5. ENRIQUECIMIENTO FUNCIONAL Y DISCUSIÓN BIOLÓGICA
### Objetivo: interpretar genes diferenciales mediante GO, KEGG,
### Reactome y revisión funcional gen a gen si el número de genes es bajo
############################################################

### 0. INSTALACIÓN, CARGA DE PAQUETES Y REPRODUCIBILIDAD ####

InsLoad.pks <- function(pqks, force = FALSE) {
     if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")    # Instala gestor BioC si falta
     }
     new.pkg <- pqks[!(pqks %in% installed.packages()[, "Package"])]    # Detecta faltantes
     if (length(new.pkg) > 0) {
          BiocManager::install(new.pkg, update = FALSE, ask = FALSE, force = force)    # Instala faltantes
     }
     loaded <- sapply(pqks, require, character.only = TRUE)    # Carga paquetes
     print(loaded)
     print("Arriba, el listado de paquetes requeridos y con confirmacion TRUE, si fueron cargados:")
     if (length(new.pkg) > 0) {
          return(print(paste("Se hizo tramite de instalacion de", length(new.pkg),
                             "paquetes nuevos:", paste(new.pkg, collapse = ", "))))
     } else {
          return(print("Igual, los paquetes requeridos ya estaban instalados"))
     }
}

pqks <- c("clusterProfiler", "ReactomePA", "org.Hs.eg.db", "AnnotationDbi",
          "enrichplot", "ggplot2", "readr", "dplyr", "tibble")    # Paquetes para enriquecimiento
InsLoad.pks(pqks)

sealseed <- function(seed = 7777777) {
     RNGkind(kind = "Mersenne-Twister",
             normal.kind = "Inversion",
             sample.kind = "Rejection")    # Configuración explícita del RNG
     set.seed(seed)
     message("Semilla fijada en: ", seed)
     invisible(seed)
}

sealseed(123)

dir.create("tables", showWarnings = FALSE)    # Tablas finales y auxiliares
dir.create("graphs", showWarnings = FALSE)    # Figuras del enriquecimiento
dir.create("objects", showWarnings = FALSE)    # Objetos de R reutilizables

### 1. CARGA DE RESULTADOS DIFERENCIALES ####

res_deseq2 <- read.csv("tables/resultados_DESeq2_Obeso1_vs_Obeso2.csv")    # Resultado principal

res_deseq2_edad <- read.csv("tables/resultados_DESeq2_edad_Obeso1_vs_Obeso2.csv")    # Modelo ajustado por edad

res_edger <- read.csv("tables/resultados_edgeR_Obeso1_vs_Obeso2.csv")    # Comparación secundaria

cat("\nGenes evaluados por DESeq2 principal:\n")
print(nrow(res_deseq2))

cat("\nGenes evaluados por edgeR:\n")
print(nrow(res_edger))

### 2. GENES CANDIDATOS PARA ENRIQUECIMIENTO Y DISCUSIÓN ####

genes_deseq2_sig <- res_deseq2[!is.na(res_deseq2$padj) &
                                    res_deseq2$padj < 0.05, ]    # Genes con FDR < 0.05

genes_deseq2_sig <- genes_deseq2_sig[order(genes_deseq2_sig$padj), ]    # Más significativos arriba

genes_candidatos <- unique(genes_deseq2_sig$gene_id)    # Genes DESeq2 que entran al enriquecimiento

universo_genes <- unique(res_deseq2$gene_id)    # Fondo: genes realmente evaluados por DESeq2

genes_edger_sig <- res_edger[!is.na(res_edger$FDR) &
                                  res_edger$FDR < 0.05, ]    # Genes significativos por edgeR

genes_nucleo <- intersect(genes_candidatos,
                          genes_edger_sig$gene_id)    # Genes reproducidos por DESeq2 y edgeR

genes_exploratorios <- union(genes_candidatos,
                             genes_edger_sig$gene_id)    # Lista amplia para lectura funcional exploratoria

if (file.exists("tables/genes_nucleo_DESeq2_edgeR.csv")) {
     genes_nucleo <- read.csv("tables/genes_nucleo_DESeq2_edgeR.csv")$gene_id
}    # Usa la lista integrada de parte 4 si ya fue generada

if (file.exists("tables/genes_exploratorios_DESeq2_edgeR.csv")) {
     genes_exploratorios <- read.csv("tables/genes_exploratorios_DESeq2_edgeR.csv")$gene_id
}    # Usa la unión integrada de parte 4 si ya fue generada

genes_nucleo <- unique(genes_nucleo[!is.na(genes_nucleo) & genes_nucleo != ""])
genes_exploratorios <- unique(genes_exploratorios[!is.na(genes_exploratorios) & genes_exploratorios != ""])

resumen_candidatos <- data.frame(
     metodo = c("DESeq2 principal", "DESeq2 ajustado por edad", "edgeR secundario"),
     genes_evaluados = c(nrow(res_deseq2), nrow(res_deseq2_edad), nrow(res_edger)),
     genes_significativos_FDR_0.05 = c(nrow(genes_deseq2_sig),
                                       sum(!is.na(res_deseq2_edad$padj) &
                                                res_deseq2_edad$padj < 0.05),
                                       nrow(genes_edger_sig)))    # Resumen para discusión

cat("\nGenes significativos por DESeq2 principal:\n")
print(genes_candidatos)

cat("\nGenes núcleo DESeq2-edgeR para enriquecimiento exploratorio:\n")
print(genes_nucleo)

cat("\nGenes exploratorios DESeq2 o edgeR:\n")
print(genes_exploratorios)

write.csv(genes_deseq2_sig,
          file = "tables/genes_candidatos_DESeq2.csv",
          row.names = FALSE)    # Genes usados como entrada principal

write.csv(genes_edger_sig,
          file = "tables/genes_candidatos_edgeR.csv",
          row.names = FALSE)    # Genes significativos por edgeR, si los hay

write.csv(resumen_candidatos,
          file = "tables/resumen_candidatos_enriquecimiento.csv",
          row.names = FALSE)    # Resumen de métodos

write.csv(data.frame(gene_id = genes_nucleo),
          file = "tables/genes_nucleo_DESeq2_edgeR.csv",
          row.names = FALSE)    # Núcleo reproducido por métodos

write.csv(data.frame(gene_id = genes_exploratorios),
          file = "tables/genes_exploratorios_DESeq2_edgeR.csv",
          row.names = FALSE)    # Lista amplia para análisis funcional exploratorio

### 3. ENRIQUECIMIENTO GO CON SÍMBOLOS GÉNICOS ####

enrich_go_bp <- clusterProfiler::enrichGO(gene = genes_candidatos,    # Genes candidatos
                                          universe = universo_genes,    # Fondo: genes evaluados
                                          OrgDb = org.Hs.eg.db::org.Hs.eg.db,    # Anotación humana
                                          keyType = "SYMBOL",    # Identificadores como símbolos génicos
                                          ont = "BP",    # Procesos biológicos
                                          pAdjustMethod = "BH",    # Corrección FDR
                                          pvalueCutoff = 0.05,    # Umbral nominal
                                          qvalueCutoff = 0.20,    # Umbral FDR flexible
                                          minGSSize = 2,    # Ajustado por universo pequeño
                                          maxGSSize = 500,    # Tamaño máximo del término
                                          readable = TRUE)    # Mantiene genes legibles

enrich_go_mf <- clusterProfiler::enrichGO(gene = genes_candidatos,    # Misma lista candidata
                                          universe = universo_genes,    # Mismo fondo génico
                                          OrgDb = org.Hs.eg.db::org.Hs.eg.db,    # Anotación humana
                                          keyType = "SYMBOL",    # Identificadores como símbolos génicos
                                          ont = "MF",    # Función molecular
                                          pAdjustMethod = "BH",    # Corrección FDR
                                          pvalueCutoff = 0.05,    # Umbral nominal
                                          qvalueCutoff = 0.20,    # Umbral FDR flexible
                                          minGSSize = 2,    # Ajustado por universo pequeño
                                          maxGSSize = 500,    # Tamaño máximo del término
                                          readable = TRUE)    # Mantiene genes legibles

enrich_go_cc <- clusterProfiler::enrichGO(gene = genes_candidatos,    # Misma lista candidata
                                          universe = universo_genes,    # Mismo fondo génico
                                          OrgDb = org.Hs.eg.db::org.Hs.eg.db,    # Anotación humana
                                          keyType = "SYMBOL",    # Identificadores como símbolos génicos
                                          ont = "CC",    # Componente celular
                                          pAdjustMethod = "BH",    # Corrección FDR
                                          pvalueCutoff = 0.05,    # Umbral nominal
                                          qvalueCutoff = 0.20,    # Umbral FDR flexible
                                          minGSSize = 2,    # Ajustado por universo pequeño
                                          maxGSSize = 500,    # Tamaño máximo del término
                                          readable = TRUE)    # Mantiene genes legibles

res_go_bp <- as.data.frame(enrich_go_bp)    # Resultado GO BP como tabla
res_go_mf <- as.data.frame(enrich_go_mf)    # Resultado GO MF como tabla
res_go_cc <- as.data.frame(enrich_go_cc)    # Resultado GO CC como tabla

write.csv(res_go_bp,
          file = "tables/enriquecimiento_GO_BP_DESeq2.csv",
          row.names = FALSE)    # Guarda GO Biological Process

write.csv(res_go_mf,
          file = "tables/enriquecimiento_GO_MF_DESeq2.csv",
          row.names = FALSE)    # Guarda GO Molecular Function

write.csv(res_go_cc,
          file = "tables/enriquecimiento_GO_CC_DESeq2.csv",
          row.names = FALSE)    # Guarda GO Cellular Component

### 4. CONVERSIÓN A ENTREZ PARA REACTOME Y KEGG ####

genes_entrez <- clusterProfiler::bitr(geneID = genes_candidatos,    # Genes candidatos en símbolo
                                      fromType = "SYMBOL",    # Formato inicial
                                      toType = "ENTREZID",    # Formato requerido por rutas
                                      OrgDb = org.Hs.eg.db::org.Hs.eg.db)    # Base humana

universo_entrez <- clusterProfiler::bitr(geneID = universo_genes,    # Todos los genes evaluados
                                         fromType = "SYMBOL",    # Formato inicial
                                         toType = "ENTREZID",    # Formato requerido por rutas
                                         OrgDb = org.Hs.eg.db::org.Hs.eg.db)    # Base humana

write.csv(genes_entrez,
          file = "tables/genes_candidatos_entrez.csv",
          row.names = FALSE)    # Equivalencia de genes candidatos

write.csv(universo_entrez,
          file = "tables/universo_genes_entrez.csv",
          row.names = FALSE)    # Equivalencia del universo usado

### 5. ENRIQUECIMIENTO REACTOME CONSERVADOR ####

enrich_reactome <- tryCatch(
     ReactomePA::enrichPathway(gene = genes_entrez$ENTREZID,    # Genes candidatos en Entrez
                               universe = universo_entrez$ENTREZID,    # Fondo evaluado en Entrez
                               organism = "human",    # Reactome para humano
                               pvalueCutoff = 0.05,    # Umbral nominal
                               pAdjustMethod = "BH",    # Corrección por múltiples pruebas
                               minGSSize = 2,    # Ajustado por pocos genes
                               maxGSSize = 500,    # Tamaño máximo
                               readable = TRUE),    # Muestra símbolos legibles
     error = function(e) {
          message("Reactome conservador no pudo ejecutarse: ", conditionMessage(e))
          NULL
     })

estado_reactome_conservador <- if (is.null(enrich_reactome)) {
     "no ejecutado o sin resultados por error"
} else {
     "ejecutado"
}

res_reactome <- if (is.null(enrich_reactome)) data.frame() else as.data.frame(enrich_reactome)    # Reactome como tabla común

write.csv(res_reactome,
          file = "tables/enriquecimiento_Reactome_DESeq2.csv",
          row.names = FALSE)    # Guarda Reactome

### 6. ENRIQUECIMIENTO KEGG CONSERVADOR ####

enrich_kegg <- tryCatch(
     clusterProfiler::enrichKEGG(gene = genes_entrez$ENTREZID,    # Genes candidatos en Entrez
                                 universe = universo_entrez$ENTREZID,    # Fondo evaluado en Entrez
                                 organism = "hsa",    # hsa: Homo sapiens
                                 pvalueCutoff = 0.05,    # Umbral nominal
                                 pAdjustMethod = "BH",    # Corrección por múltiples pruebas
                                 qvalueCutoff = 0.20,    # Umbral FDR flexible
                                 minGSSize = 2,    # Ajustado por universo pequeño
                                 maxGSSize = 500),    # Tamaño máximo
     error = function(e) {
          message("KEGG conservador no pudo ejecutarse: ", conditionMessage(e))
          NULL
     })

estado_kegg_conservador <- if (is.null(enrich_kegg)) {
     "no ejecutado por conexion o error"
} else {
     "ejecutado"
}

res_kegg <- if (is.null(enrich_kegg)) data.frame() else as.data.frame(enrich_kegg)    # KEGG como tabla común

write.csv(res_kegg,
          file = "tables/enriquecimiento_KEGG_DESeq2.csv",
          row.names = FALSE)    # Guarda KEGG

### 6B. ENRIQUECIMIENTO EXPLORATORIO CON FONDO AMPLIO ####

genes_exploratorios_entrada <- if (length(genes_exploratorios) >= 2) {
     genes_exploratorios
} else {
     genes_nucleo
}    # Lista exploratoria: genes significativos por DESeq2 o edgeR

genes_entrez_exploratorios <- clusterProfiler::bitr(
     geneID = genes_exploratorios_entrada,
     fromType = "SYMBOL",
     toType = "ENTREZID",
     OrgDb = org.Hs.eg.db::org.Hs.eg.db)    # Conversión para rutas exploratorias

genes_entrez_exploratorios <- genes_entrez_exploratorios[
     !duplicated(genes_entrez_exploratorios$ENTREZID), ]    # Evita duplicados por mapeo

enrich_go_bp_exploratorio <- clusterProfiler::enrichGO(
     gene = genes_exploratorios_entrada,
     OrgDb = org.Hs.eg.db::org.Hs.eg.db,
     keyType = "SYMBOL",
     ont = "BP",
     pAdjustMethod = "BH",
     pvalueCutoff = 0.10,
     qvalueCutoff = 0.20,
     minGSSize = 2,
     maxGSSize = 500,
     readable = TRUE)    # Fondo amplio: lectura exploratoria, no evidencia causal

res_go_bp_exploratorio <- as.data.frame(enrich_go_bp_exploratorio)

write.csv(res_go_bp_exploratorio,
          file = "tables/enriquecimiento_GO_BP_exploratorio_DESeq2_edgeR.csv",
          row.names = FALSE)    # GO BP exploratorio para discusión funcional

enrich_reactome_exploratorio <- tryCatch(
     ReactomePA::enrichPathway(gene = genes_entrez_exploratorios$ENTREZID,
                               organism = "human",
                               pvalueCutoff = 0.10,
                               pAdjustMethod = "BH",
                               minGSSize = 2,
                               maxGSSize = 500,
                               readable = TRUE),
     error = function(e) {
          message("Reactome exploratorio no pudo ejecutarse: ", conditionMessage(e))
          NULL
     })

estado_reactome_exploratorio <- if (is.null(enrich_reactome_exploratorio)) {
     "no ejecutado o sin resultados por error"
} else {
     "ejecutado"
}

res_reactome_exploratorio <- if (is.null(enrich_reactome_exploratorio)) {
     data.frame()
} else {
     as.data.frame(enrich_reactome_exploratorio)
}

write.csv(res_reactome_exploratorio,
          file = "tables/enriquecimiento_Reactome_exploratorio_DESeq2_edgeR.csv",
          row.names = FALSE)    # Reactome exploratorio

enrich_kegg_exploratorio <- tryCatch(
     clusterProfiler::enrichKEGG(gene = genes_entrez_exploratorios$ENTREZID,
                                 organism = "hsa",
                                 pvalueCutoff = 0.10,
                                 pAdjustMethod = "BH",
                                 qvalueCutoff = 0.20,
                                 minGSSize = 2,
                                 maxGSSize = 500),
     error = function(e) {
          message("KEGG exploratorio no pudo ejecutarse: ", conditionMessage(e))
          NULL
     })

estado_kegg_exploratorio <- if (is.null(enrich_kegg_exploratorio)) {
     "no ejecutado por conexion o error"
} else {
     "ejecutado"
}

res_kegg_exploratorio <- if (is.null(enrich_kegg_exploratorio)) {
     data.frame()
} else {
     as.data.frame(enrich_kegg_exploratorio)
}

write.csv(res_kegg_exploratorio,
          file = "tables/enriquecimiento_KEGG_exploratorio_DESeq2_edgeR.csv",
          row.names = FALSE)    # KEGG exploratorio

if (nrow(res_go_bp_exploratorio) > 0) {
     grafico_go_exploratorio <- enrichplot::dotplot(enrich_go_bp_exploratorio,
                                                    showCategory = 12) +
          ggplot2::ggtitle("GO BP exploratorio, DESeq2-edgeR")

     ggplot2::ggsave(filename = "graphs/enriquecimiento_GO_BP_exploratorio_DESeq2_edgeR.png",
                     plot = grafico_go_exploratorio,
                     width = 9,
                     height = 7,
                     dpi = 150)    # Figura para póster

     ggplot2::ggsave(filename = "graphs/enriquecimiento_GO_BP_exploratorio_DESeq2_edgeR.pdf",
                     plot = grafico_go_exploratorio,
                     width = 9,
                     height = 7)    # Versión vectorial
}

### 7. TABLA PARA DISCUSIÓN GEN A GEN ####

discusion_genes <- genes_deseq2_sig[, c("gene_id", "baseMean", "log2FoldChange", "pvalue", "padj")]    # Columnas clave

discusion_genes$direccion <- ifelse(discusion_genes$log2FoldChange > 0,
                                    "Mayor en Obeso1",
                                    "Mayor en Obeso2")    # Sentido del cambio

discusion_genes$interpretacion_biologica <- ""    # Espacio para curaduría manual

discusion_genes$fuente_sugerida <- "GeneCards/PubMed/OMIM/Reactome/KEGG"    # Bases para revisión

write.csv(discusion_genes,
          file = "tables/tabla_discusion_gen_a_gen.csv",
          row.names = FALSE)    # Tabla para redactar discusión

### 8. GRÁFICOS DE ENRIQUECIMIENTO ####

cat("\nNúmero de términos conservadores por categoría:\n")
cat("GO BP:", nrow(res_go_bp), "\n")
cat("GO MF:", nrow(res_go_mf), "\n")
cat("GO CC:", nrow(res_go_cc), "\n")
cat("Reactome:", nrow(res_reactome), "\n")
cat("KEGG:", nrow(res_kegg), "\n")

cat("\nNúmero de términos exploratorios por categoría:\n")
cat("GO BP exploratorio:", nrow(res_go_bp_exploratorio), "\n")
cat("Reactome exploratorio:", nrow(res_reactome_exploratorio), "\n")
cat("KEGG exploratorio:", nrow(res_kegg_exploratorio), "\n")

### 9. RESUMEN Y OBJETOS DE ENRIQUECIMIENTO ####

resumen_enriquecimiento <- data.frame(
     enfoque = c("Conservador", "Conservador", "Conservador", "Conservador", "Conservador",
                 "Exploratorio", "Exploratorio", "Exploratorio"),
     analisis = c("GO_BP", "GO_MF", "GO_CC", "Reactome", "KEGG",
                  "GO_BP", "Reactome", "KEGG"),
     universo = c(rep("37 genes evaluados por DESeq2", 5),
                  rep("Fondo amplio de anotación", 3)),
     estado = c(rep("ejecutado", 3),
                estado_reactome_conservador,
                estado_kegg_conservador,
                "ejecutado",
                estado_reactome_exploratorio,
                estado_kegg_exploratorio),
     terminos_enriquecidos = c(nrow(res_go_bp), nrow(res_go_mf), nrow(res_go_cc),
                               nrow(res_reactome), nrow(res_kegg),
                               nrow(res_go_bp_exploratorio),
                               nrow(res_reactome_exploratorio),
                               nrow(res_kegg_exploratorio)))    # Conteo de resultados por enfoque

write.csv(resumen_enriquecimiento,
          file = "tables/resumen_enriquecimiento_funcional.csv",
          row.names = FALSE)    # Resumen para auditoría

nota_enriquecimiento <- data.frame(
     punto = c("Resultado conservador", "Resultado exploratorio", "Interpretación"),
     texto = c("Con universo restringido a los 37 genes evaluados, el análisis puede no detectar términos enriquecidos.",
               "Con fondo amplio de anotación, pueden aparecer rutas coherentes con obesidad y señalización hormonal.",
               "El enriquecimiento se usa como apoyo funcional descriptivo; no constituye evidencia causal, incluso si resulta positivo. KEGG depende de conexión externa y puede quedar no ejecutado."))    # Salvedad para póster

write.csv(nota_enriquecimiento,
          file = "tables/nota_interpretacion_enriquecimiento.csv",
          row.names = FALSE)    # Texto base para discusión

saveRDS(object = enrich_go_bp,
        file = "objects/enrich_go_bp.rds")    # Objeto GO BP

saveRDS(object = enrich_go_mf,
        file = "objects/enrich_go_mf.rds")    # Objeto GO MF

saveRDS(object = enrich_go_cc,
        file = "objects/enrich_go_cc.rds")    # Objeto GO CC

saveRDS(object = enrich_reactome,
        file = "objects/enrich_reactome.rds")    # Objeto Reactome

saveRDS(object = enrich_kegg,
        file = "objects/enrich_kegg.rds")    # Objeto KEGG

saveRDS(object = enrich_go_bp_exploratorio,
        file = "objects/enrich_go_bp_exploratorio.rds")    # Objeto GO BP exploratorio

saveRDS(object = enrich_reactome_exploratorio,
        file = "objects/enrich_reactome_exploratorio.rds")    # Objeto Reactome exploratorio

saveRDS(object = enrich_kegg_exploratorio,
        file = "objects/enrich_kegg_exploratorio.rds")    # Objeto KEGG exploratorio

### REVISIÓN RÁPIDA DE RESULTADOS DE ENRIQUECIMIENTO ####



cat("\nPARTE 5 FINALIZADA\n")
cat("Se ejecutó enriquecimiento conservador y exploratorio.\n")
cat("El resultado exploratorio se interpreta como apoyo funcional descriptivo, no como evidencia causal.\n")
