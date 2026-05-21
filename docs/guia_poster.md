# Guia breve para el poster

## Metodologia

Se analizaron FASTQ simulados de RNA-seq para comparar los grupos Sobrepeso/Obeso1 y Sobrepeso/Obeso2. El flujo principal incluyo control de calidad con FastQC/MultiQC, cuantificacion directa con Salmon sobre lecturas crudas, agrupacion transcrito-gen con `tximport`, analisis diferencial con DESeq2 y contraste complementario con edgeR. No se aplico trimming en el analisis principal porque los datos son simulados y presentaron calidad suficiente; el trimming se conserva como analisis de sensibilidad.

## Resultados principales

DESeq2 identifico cinco genes diferencialmente expresados en el analisis principal sin trimming: `FTO` y `MC4R` con mayor expresion relativa en Obeso2, y `LEPR`, `SH2B1` y `POMC` con mayor expresion relativa en Obeso1. edgeR se uso como contraste secundario para evaluar consistencia metodologica y amplio la lista de genes candidatos para la interpretacion exploratoria.

## Enriquecimiento

El enriquecimiento conservador, usando como universo los 37 genes evaluados, no detecto terminos funcionales significativos. En cambio, el enriquecimiento exploratorio con fondo amplio recupero procesos GO y rutas Reactome coherentes con apetito, senalizacion hormonal, leptina, melanocortinas y metabolismo energetico. Este resultado debe presentarse como apoyo funcional descriptivo y no como evidencia causal.

## Discusion

Los datos son simulados y el panel de genes ya esta dirigido a obesidad, por lo que la interpretacion debe ser prudente. El valor del analisis no esta en demostrar causalidad biologica, sino en mostrar que el flujo bioinformatico permite cuantificar expresion, resumir transcritos a genes, detectar diferencias entre grupos y conectar los genes principales con funciones biologicas plausibles.

## Conclusion sugerida

El pipeline principal permitio diferenciar perfiles transcriptomicos simulados entre Obeso1 y Obeso2. Los genes identificados se relacionan con regulacion del apetito, senalizacion de leptina y control neuroendocrino del balance energetico. Las rutas enriquecidas apoyan la interpretacion funcional de forma exploratoria, con resultados sensibles al universo de referencia y al caracter simulado del conjunto de datos.
