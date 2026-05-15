# Grupo 13 Simpson RNA-seq obesidad

Repositorio del analisis de expresion diferencial de genes relacionados con obesidad usando datos simulados de RNA-seq.

## Objetivo

Procesar archivos FASTQ simulados, realizar control de calidad, cuantificacion con Salmon, agrupacion de transcritos a genes, analisis diferencial y visualizacion de resultados para la comparativa Obeso1 vs Obeso2.

## Scripts

- `parte1yParte2.sh`: entorno Conda, control de calidad, limpieza de lecturas, construcción del índice de Salmon y cuantificación por muestra.
- `parte3.r`: agrupación de transcritos a genes usando `Transcrito_a_Gen.tsv` y los archivos `quant.sf` de Salmon.
- `Parte4.R`: análisis diferencial Obeso1 vs Obeso2 con DESeq2 y comparación secundaria con edgeR.

## Nota
Los archivos FASTQ incluidos corresponden a datos simulados con fines docentes. Las salidas intermedias regenerables, como reportes FastQC/MultiQC, lecturas limpias, índices de Salmon y carpetas completas de cuantificación, pueden reconstruirse ejecutando los scripts del pipeline.