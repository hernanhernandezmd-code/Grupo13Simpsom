# Grupo 13 Simpson - RNA-seq (Obeso1 vs Obeso2)

Repositorio del analisis grupal de expresion diferencial de genes relacionados con obesidad usando datos simulados de RNA-seq.

## Integrantes

- Hernán Guillermo Hernández
- Frank Melo Restrepo
- Carlos Hernán Sierra Torres
- Andrés Ricardo López Preciado
- Julián Andrés Silva Rojas

## Objetivo

Ejecutar un flujo reproducible para comparar los perfiles de expresion entre los grupos Sobrepeso/Obeso1 y Sobrepeso/Obeso2. El pipeline principal realiza control de calidad, cuantificacion directa con Salmon, agrupacion transcrito-gen, analisis diferencial, visualizacion e interpretacion funcional.

## Flujo principal

El analisis principal sigue esta ruta:

1. FastQC y MultiQC sobre FASTQ crudos.
2. Cuantificacion directa con Salmon, sin trimming previo.
3. Agrupacion de transcritos a genes con `tximport` y `Transcrito_a_Gen.tsv`.
4. Analisis diferencial con DESeq2.
5. Comparacion secundaria con edgeR.
6. Visualizaciones finales: volcano plot, heatmap, PCA y tablas para poster.
7. Enriquecimiento funcional conservador y exploratorio.

El trimming se conserva solo como analisis complementario o de sensibilidad, porque los FASTQ son simulados y presentan calidad suficiente para el flujo principal.

## Estructura

- `Fastqs/`: lecturas FASTQ simuladas paired-end.
- `Genes/`: archivos de referencia por gen incluidos en la actividad.
- `1_fastqc/`: reportes individuales de FastQC.
- `2_multiqc/`: reporte integrado de MultiQC.
- `4_indice_transcriptoma/`: indice principal de Salmon construido desde `Referencia.fasta`.
- `5_cuantificacion_salmon/`: cuantificaciones principales de Salmon (`quant.sf`) por muestra.
- `3_lecturas_limpias/`: lecturas generadas solo si se ejecuta el analisis complementario con trimming.
- `tables/`: matrices, resultados diferenciales, enriquecimiento y tablas para poster.
- `graphs/`: figuras finales en PNG/PDF.
- `objects/`: objetos `.rds` para reutilizar resultados intermedios.
- `Scripts/paracomplementarloprincipal/`: material complementario que no pertenece al pipeline principal.
- `docs/`: textos breves para metodologia, discusion y conclusiones del poster.
- `Design.csv`: metadatos de muestra, condicion, edad y sexo.
- `Referencia.fasta`: transcriptoma de referencia usado por Salmon.
- `Transcrito_a_Gen.tsv`: correspondencia transcrito-gen para `tximport`.

## Scripts

- `Scripts/parte1yParte2.sh`: entorno, FastQC, MultiQC, indice Salmon y cuantificacion directa sin trimming.
- `Scripts/parte3.r`: importacion de `quant.sf` con `tximport` y matrices por gen.
- `Scripts/parte4.R`: DESeq2, modelo ajustado por edad, edgeR y tablas comparativas.
- `Scripts/parte5.R`: enriquecimiento conservador y exploratorio, mas tabla para discusion gen a gen.
- `Scripts/parte6.R`: volcano plot, heatmap, PCA y tablas finales para poster.
- `Scripts/paracomplementarloprincipal/analisis_complementario_trimming.sh`: flujo opcional con Trimmomatic para evaluar sensibilidad al preprocesamiento.
- `Scripts/paracomplementarloprincipal/Frank_Melo_Restrepo_Actividad2.R`: aporte complementario/exploratorio de Frank.
- `docs/guia_poster.md`: propuesta breve de texto para poster.

## Ejecucion

Desde la raiz del repositorio:

```bash
bash Scripts/parte1yParte2.sh
Rscript Scripts/parte3.r
Rscript Scripts/parte4.R
Rscript Scripts/parte5.R
Rscript Scripts/parte6.R
```

El analisis complementario con trimming se ejecuta aparte:

```bash
bash Scripts/paracomplementarloprincipal/analisis_complementario_trimming.sh
```

## Notas metodologicas

### Trimming

El flujo principal no aplica trimming. Esta decision esta alineada con el caracter simulado de los datos y con la buena calidad global observada en el control inicial. El trimming se mantiene como analisis complementario, no como base del resultado presentado en el poster.

### Normalizacion

DESeq2 y edgeR modelan conteos crudos mediante distribucion binomial negativa e incorporan normalizacion interna por profundidad de libreria. Por eso, las pruebas diferenciales se hacen sobre conteos no transformados. Las matrices normalizadas y TPM se usan para descripcion, tablas y visualizacion.

### Enriquecimiento funcional

Se reportan dos enfoques:

- Conservador: usa como universo los 37 genes evaluados. Este enfoque no detecto enriquecimiento funcional significativo.
- Exploratorio: usa genes significativos por DESeq2 o edgeR con fondo amplio de anotacion. Este enfoque recupera terminos biologicamente coherentes con obesidad, apetito y senalizacion hormonal.

El enriquecimiento se interpreta como apoyo funcional descriptivo, no como evidencia causal. Incluso cuando aparecen terminos positivos, el resultado depende del universo de referencia, del tamano muestral y del caracter simulado del conjunto de datos.

Nota: KEGG consulta recursos online desde R. Si no hay conexion o DNS disponible, el script deja KEGG en cero y continua con GO/Reactome.

## Resultados clave

El contraste principal es `Obeso1 vs Obeso2`. Un log2FoldChange positivo indica mayor expresion en Obeso1; un valor negativo indica mayor expresion en Obeso2.

Tablas principales:

- `tables/resultados_DESeq2_Obeso1_vs_Obeso2.csv`
- `tables/resultados_edgeR_Obeso1_vs_Obeso2.csv`
- `tables/comparacion_DESeq2_edgeR_Obeso1_vs_Obeso2.csv`
- `tables/tabla_poster_genes_DESeq2.csv`
- `tables/resumen_enriquecimiento_funcional.csv` (incluye estado de ejecucion de KEGG/Reactome)
- `tables/nota_interpretacion_enriquecimiento.csv`

Figuras principales:

- `graphs/volcano_DESeq2_Obeso1_vs_Obeso2.png`
- `graphs/heatmap_genes_significativos_DESeq2.png`
- `graphs/PCA_DESeq2_Obeso1_vs_Obeso2.png`
- `graphs/enriquecimiento_GO_BP_exploratorio_DESeq2_edgeR.png`
