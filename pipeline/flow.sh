#!/bin/bash
set -e

DATA_DIR="/home/rachel/cattlealgua/X204SC22114535-Z01-F001_02"
BEFORE_DATA="$DATA_DIR/before_data/cleaned_files"
AFTER_DATA="$DATA_DIR/after_data/cleaned_files"
OUT_SAM="/home/rachel/out_sam"
BAM_FASTQ="/home/rachel/bamtofastq"
KRAKEN_DIR="/home/rachel/kraken/kraken2"
KRAKEN_DB="$KRAKEN_DIR/kraken2db"
KREPORT_DIR="/home/rachel/kreport"
KRONA_DIR="/home/rachel/kronafile"
KRONA_HTML="/home/rachel/kronahtml"

source /home/rachel/anaconda3/etc/profile.d/conda.sh

for SAMPLE in {13..24}; do
    echo "Processing sample A$SAMPLE..."
    
    if [[ $SAMPLE -ge 13 && $SAMPLE -le 18 ]]; then
        CURRENT_BEFORE_DATA="$BEFORE_DATA"
    else
        CURRENT_BEFORE_DATA="$AFTER_DATA"
    fi

    conda activate rachel-env
    bowtie2 --very-sensitive-local --threads 30 -k 1 -x /home/rachel/human_index \
        -1 "$CURRENT_BEFORE_DATA/trimmed_A${SAMPLE}_1.fastq" \
        -2 "$CURRENT_BEFORE_DATA/trimmed_A${SAMPLE}_2.fastq" \
        -S "$OUT_SAM/A${SAMPLE}.sam"
    conda deactivate

    samtools view -f 4 -bS "$OUT_SAM/A${SAMPLE}.sam" -o "$BAM_FASTQ/unaligned_human_A${SAMPLE}.bam"
    samtools fastq -n -1 "$BAM_FASTQ/unaligned_human_A${SAMPLE}_1.fastq" \
                   -2 "$BAM_FASTQ/unaligned_human_A${SAMPLE}_2.fastq" \
                   -s "$BAM_FASTQ/unaligned_human_A${SAMPLE}_singletons.fastq" \
                   "$BAM_FASTQ/unaligned_human_A${SAMPLE}.bam"

    conda activate rachel-env
    bowtie2 --very-sensitive-local --threads 30 -k 1 -x /home/rachel/bos_taurus_index \
        -1 "$BAM_FASTQ/unaligned_human_A${SAMPLE}_1.fastq" \
        -2 "$BAM_FASTQ/unaligned_human_A${SAMPLE}_2.fastq" \
        -S "$OUT_SAM/A${SAMPLE}_bis.sam"
    conda deactivate

    samtools view -f 4 -bS "$OUT_SAM/A${SAMPLE}_bis.sam" > "$BAM_FASTQ/unaligned_human_cow_A${SAMPLE}.bam"
    samtools fastq -n -1 "$BAM_FASTQ/unaligned_human_cow_A${SAMPLE}_1.fastq" \
                   -2 "$BAM_FASTQ/unaligned_human_cow_A${SAMPLE}_2.fastq" \
                   -s "$BAM_FASTQ/unaligned_human_cow_A${SAMPLE}_singletons.fastq" \
                   "$BAM_FASTQ/unaligned_human_cow_A${SAMPLE}.bam"

    cd "$KRAKEN_DIR"
    ./kraken2 --db "$KRAKEN_DB" --threads 30 --output "A${SAMPLE}_taxonomy" \
            --report "$KREPORT_DIR/A${SAMPLE}_taxonomy_kreport" \
            --paired "$BAM_FASTQ/unaligned_human_cow_A${SAMPLE}_1.fastq" \
                     "$BAM_FASTQ/unaligned_human_cow_A${SAMPLE}_2.fastq"

    conda activate rachel-env
    python "$KRAKEN_DIR/KrakenTools/modifkreport2krona.py" -r "$KREPORT_DIR/A${SAMPLE}_taxonomy_kreport" -o "$KRONA_DIR/A${SAMPLE}_krona"
    ktImportText "$KRONA_DIR/A${SAMPLE}_krona" -o "$KRONA_HTML/A${SAMPLE}.html"
    conda deactivate

done

echo "Pipeline completed for samples A13 to A24."
