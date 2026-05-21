#!/bin/bash
set -Eeuo pipefail

############################################################
### Proyecto: Grupo 13 Simpson RNA-seq obesidad
### PARTE 1 Y PARTE 2. ENTORNO, CONTROL DE CALIDAD Y SALMON
### Parte 1: entorno Conda y herramientas bioinformáticas
### Parte 2: FastQC, MultiQC y cuantificación directa con Salmon
############################################################

### 0. CONFIGURACIÓN GENERAL ####

unset LANGUAGE
export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8
export LC_CTYPE=en_US.UTF-8

ENV_NAME="analisis_obesidad"
THREADS=4
IDX="4_indice_transcriptoma/transcriptoma_index"

### 1. ENTORNO CONDA Y PROGRAMAS ####

command -v conda >/dev/null 2>&1 || { echo "ERROR: Conda no esta disponible en esta terminal."; exit 1; }
eval "$(conda shell.bash hook)"    # Permite activar Conda dentro del script

if ! conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
     conda create -n "$ENV_NAME" -c conda-forge -c bioconda \
          fastqc multiqc salmon python=3.10 -y    # Crea entorno si no existe
fi

conda activate "$ENV_NAME"    # Activa el entorno del análisis

unset LANGUAGE
export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8
export LC_CTYPE=en_US.UTF-8

### 2. INSUMOS Y CARPETAS DE TRABAJO ####

[[ -d "Fastqs" ]] || { echo "ERROR: No existe la carpeta Fastqs en $(pwd)."; exit 1; }
[[ -f "Referencia.fasta" ]] || { echo "ERROR: No existe Referencia.fasta en $(pwd)."; exit 1; }
[[ -f "Transcrito_a_Gen.tsv" ]] || { echo "ERROR: No existe Transcrito_a_Gen.tsv en $(pwd)."; exit 1; }

shopt -s nullglob
r1_files=(Fastqs/*_R1.fastq.gz)
[[ ${#r1_files[@]} -gt 0 ]] || { echo "ERROR: No se encontraron archivos *_R1.fastq.gz."; exit 1; }

rm -rf 1_fastqc 2_multiqc 4_indice_transcriptoma 5_cuantificacion_salmon    # Limpia salidas principales previas
mkdir -p 1_fastqc 2_multiqc 4_indice_transcriptoma 5_cuantificacion_salmon    # Crea salidas principales nuevas

echo "Inicio del análisis en: $(pwd)"
echo "Entorno activo: $CONDA_DEFAULT_ENV"
echo "Salmon: $(salmon --version)"


### FIN DE LA PARTE 1  ####
### INICIO DE LA PARTE 2  ####
### 3. CONTROL DE CALIDAD INICIAL ####

fastqc Fastqs/*.fastq.gz -o 1_fastqc/    # Revisa calidad de lecturas crudas
multiqc 1_fastqc/ -o 2_multiqc/ --force    # Resume reportes FastQC

### 4. DECISIÓN DE PREPROCESAMIENTO ####

echo "Analisis principal sin trimming: FASTQ simulados de buena calidad, segun control FastQC/MultiQC."
echo "El trimming conservador queda como analisis complementario en Scripts/paracomplementarloprincipal/."

### 5. ÍNDICE DEL TRANSCRIPTOMA CON SALMON ####

rm -rf "$IDX"    # Evita reutilizar un índice incompleto

env -u LANGUAGE LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 LC_CTYPE=en_US.UTF-8 salmon index \
     -t Referencia.fasta \
     -i "$IDX" \
     -k 31    # Construye índice desde la referencia

[[ -f "$IDX/versionInfo.json" ]] || {
     echo "ERROR: Salmon no genero versionInfo.json."
     echo "Revise si Referencia.fasta es valido y si salmon index termino sin errores."
     exit 1
}

### 6. CUANTIFICACIÓN POR MUESTRA CON SALMON ####

for r1 in "${r1_files[@]}"; do
     base=$(basename "$r1" _R1.fastq.gz)    # Nombre base de la muestra
     r2="Fastqs/${base}_R2.fastq.gz"    # Pareja reversa
     [[ -f "$r2" ]] || { echo "ERROR: Falta pareja R2 para $base."; exit 1; }

     env -u LANGUAGE LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 LC_CTYPE=en_US.UTF-8 salmon quant \
          -i "$IDX" \
          -l A \
          -1 "$r1" \
          -2 "$r2" \
          -p "$THREADS" \
          --validateMappings \
          -o "5_cuantificacion_salmon/${base}"    # Genera quant.sf por muestra desde lecturas crudas

     [[ -f "5_cuantificacion_salmon/${base}/quant.sf" ]] || {
          echo "ERROR: No se genero quant.sf para $base."
          exit 1
     }
done

### 7. REPORTE FINAL Y ARCHIVOS GENERADOS ####

multiqc 1_fastqc 5_cuantificacion_salmon -o 2_multiqc --force    # Integra QC y Salmon

echo "Proceso terminado sin errores."
echo "Archivos quant.sf generados:"
find 5_cuantificacion_salmon -name quant.sf | sort
echo "Siguiente paso: importar quant.sf en R y agrupar transcritos a genes con Transcrito_a_Gen.tsv."
### FIN DE LA PARTE 2: FastQC, MultiQC y cuantificación directa con Salmon finalizados exitosamente ###
