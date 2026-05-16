# Grupo 13 Simpson — RNA‑seq (obesidad)

Repositorio del análisis de expresión diferencial de genes relacionados con obesidad usando datos simulados de RNA‑seq.

## Objetivo

Procesar los FASTQ simulados, ejecutar control de calidad, cuantificación con Salmon, agregar transcritos a genes, realizar análisis diferencial y producir visualizaciones para la comparación Obeso1 vs Obeso2.

## Estructura del repositorio
- `Fastqs/`: lecturas FASTQ simuladas (entradas del pipeline).
- `Genes/`: archivos FASTA por gen usados como referencia en el ejercicio.
- `1_fastqc/`: reportes FastQC por muestra.
- `2_multiqc/`: reporte MultiQC y datos agregados.
- `3_lecturas_limpias/`: lecturas tras el filtrado/recorte.
- `4_indice_transcriptoma/`: índice de Salmon generado desde `Referencia.fasta`.
- `5_cuantificacion_salmon/`: resultados de Salmon (`quant.sf`) por muestra.
- `tables/`: tablas de salida y matrices usadas para análisis y visualización. Contiene, entre otros:
	- `matriz_conteos_genes_todas_muestras.csv`: matriz de conteos brutos (por gen y por muestra).
	- `matriz_TPM_genes_todas_muestras.csv`: matriz de TPM calculada a partir de las cuantificaciones.
	- `matriz_conteos_normalizados_DESeq2.csv`: matriz de conteos normalizados (salida generada con DESeq2, para visualización).
- `objects/`: objetos R (`.rds`) para reutilizar resultados intermedios (p. ej. `txi_genes.rds`, `ajuste_edger.rds`).
- `results/`: salidas finales del análisis y figuras.
- `tables/`: (repetido) tablas intermedias y finales.
- `Design.csv`: diseño experimental con columnas como `sample`, `condition`, `age`, `sex`.
- `Referencia.fasta`: secuencias de referencia usadas para construir el índice.
- `Transcrito_a_Gen.tsv`: mapeo transcrito → gen usado para agregar abundancias.

## Scripts

- `parte1yParte2.sh`: prepara entorno Conda y ejecuta pasos iniciales (FastQC, limpieza, índice y cuantificación).
- `parte3.r`: importa `quant.sf`, usa `Transcrito_a_Gen.tsv` y genera matrices por gen.
- `Parte4.R`: análisis diferencial (DESeq2 principal, análisis ajustado por edad y comparación con edgeR).

## Conteos normalizados vs modelado estadístico

- Se incluye en `tables/matriz_conteos_normalizados_DESeq2.csv` una matriz de conteos normalizados generada con las funciones de DESeq2 (`counts(dds, normalized=TRUE)`). Este archivo se proporciona para facilitar visualizaciones y exportes (p. ej. generación de heatmaps o tablas por gen).
- Importante: los modelos estadísticos de DESeq2 y edgeR se construyen a partir de los conteos brutos y de los factores/size factors apropiados; la matriz normalizada aquí incluida NO sustituye al uso correcto de los conteos crudos dentro del flujo de modelado estadístico (es decir, no se recomienda usar la matriz normalizada como entrada directa para recalcular pruebas de DE).

## Nota sobre reproducibilidad

Las salidas regenerables (FastQC/MultiQC, lecturas limpias, índices de Salmon, cuantificación por muestra) pueden reconstruirse ejecutando los scripts incluidos. Los objetos `.rds` almacenados en `objects/` permiten retomar análisis intermedios sin volver a ejecutar pasos costosos.

Si necesitas que regenere la matriz normalizada a partir de los objetos R originales o que ejecute el pipeline para producir nuevas salidas, indícalo y lo preparo.