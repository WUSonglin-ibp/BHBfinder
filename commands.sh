#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------------------------------------------------
# Edit these values for your settings
# -----------------------------------------------------------------------------

# Paths and files
WORK_DIR=/mnt/e/lab/circular/manual/BHBfinder_v1/demo # Working directory
SPECIES=S_solfataricus # Species name to be used in output files
GENOME_DIR=${WORK_DIR}/genome # directory for genome files
GENOME=${GENOME_DIR}/GCF_000007005.1_ASM700v1_genomic.fna # genome sequencing file
GFF=${GENOME_DIR}/genomic.gff  # genome annotation file in gff format
HISAT2_INDEX=${GENOME_DIR}/hisat2 # output database name for HISAT2
BLAST_DB=${GENOME_DIR}/short_blast # output database name for BLAST

# List of samples for the same species. Sample directory name = dataset name without extension.
SAMPLES=(Ssol_data1 Ssol_data2) 
SPECIFIC=true #  Strand-specific or nonstrand-specific RNA-seq data 

# BHBfinder Path
BHBfinder_DIR=/mnt/e/lab/circular/manual/BHBfinder_v1/ # YOUR PATH TO BHBfinder
FIND_CIRC_SCRIPT=${BHBfinder_DIR}/scripts/find_circrna_v3.py
MERGE_JUNCTION_SCRIPT=${BHBfinder_DIR}/scripts/merge_junction.py
GET_SET2_SCRIPT=${BHBfinder_DIR}/scripts/get_set2.py

# Threads
THREADS=8
BLAST_THREADS=20


cd "${WORK_DIR}"

# -----------------------------------------------------------------------------
# Step 0: Data preparation
# -----------------------------------------------------------------------------
# Prepare RNA-seq data obtained from the NCBI Sequence Read Archive (SRA).
# 1. If necessary, install the SRA Toolkit
#      conda install -c bioconda sra-tools
# 2. Download each SRA Run into a directory named after its accession:
#      prefetch SRR26999894 -O ./Download_dir/SRR26999894
# 3. Convert each downloaded Run into FASTA files
#      fastq-dump --split-3 --fasta --threads "${THREADS}" --outdir ./Download_dir \
#          ./Download_dir/SRR26999894
# 4. For paired-end stranded RNA-seq data, determine which FASTA file contains
#    reads that have the same orientation as the source RNA.
# 5. Copy the selected FASTA file to the sample directory expected by BHBfinder:
#      ${WORK_DIR}/${SAMPLE}/${SAMPLE}.fasta



# -----------------------------------------------------------------------------
# Step 1: normal_alignment
# -----------------------------------------------------------------------------
if ! ls "${HISAT2_INDEX}"*.ht2 >/dev/null 2>&1; then
    hisat2-build "${GENOME}" "${HISAT2_INDEX}"
fi

for SAMPLE in "${SAMPLES[@]}"; do
    SAMPLE_DIR="${WORK_DIR}/${SAMPLE}"
    pushd "${SAMPLE_DIR}" >/dev/null

    hisat2 -p "${THREADS}" -f --no-spliced-alignment -x "${HISAT2_INDEX}" -U "${SAMPLE}.fasta" -S "${SAMPLE}_normal.sam"
    samtools sort "${SAMPLE}_normal.sam" > "${SAMPLE}_normal_sorted.bam"
    samtools flagstat "${SAMPLE}_normal_sorted.bam" > "${SAMPLE}_hisat2_normal_alignment.log"
    samtools view -h -b -F 4 "${SAMPLE}_normal_sorted.bam" > "${SAMPLE}_normal_m.bam"
    samtools index "${SAMPLE}_normal_m.bam"
    bedtools genomecov -dz -ibam "${SAMPLE}_normal_m.bam" > "${SAMPLE}_normal.bedgraph"
    rm "${SAMPLE}_normal.sam"

    popd >/dev/null
done

# -----------------------------------------------------------------------------
# Step 2: spliced_alignment
# -----------------------------------------------------------------------------
if ! ls "${BLAST_DB}.nhr" >/dev/null 2>&1; then
    makeblastdb -in "${GENOME}" -dbtype nucl -out "${BLAST_DB}"
fi

for SAMPLE in "${SAMPLES[@]}"; do
    SAMPLE_DIR="${WORK_DIR}/${SAMPLE}"
    pushd "${SAMPLE_DIR}" >/dev/null

    blastn -db "${BLAST_DB}" -num_threads "${BLAST_THREADS}" -task "blastn-short" -evalue 1e-3 -strand both -query "${SAMPLE}.fasta" -out "${SAMPLE}_R1.tab" -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen"

    popd >/dev/null
done

# -----------------------------------------------------------------------------
# Step 3: junction
# -----------------------------------------------------------------------------
for SAMPLE in "${SAMPLES[@]}"; do
    SAMPLE_DIR="${WORK_DIR}/${SAMPLE}"
    pushd "${SAMPLE_DIR}" >/dev/null

    python "${FIND_CIRC_SCRIPT}" -a circ -basename "${SAMPLE}" -r1 "${SAMPLE}.fasta" -tab1 "${SAMPLE}_R1.tab"
    python "${FIND_CIRC_SCRIPT}" -a junction --anti -basename "${SAMPLE}" -g "${GENOME}" -gff "${GFF}" -c "${SAMPLE}_circ.out" -j1 "${SAMPLE}_circ_R1.fasta" -j2 "${SAMPLE}_circ_R2.fasta" -bedgraph "${SAMPLE}_normal.bedgraph"

    hisat2 -f --no-spliced-alignment -x "${HISAT2_INDEX}" -1 "${SAMPLE}_BHB_R1.fasta" -2 "${SAMPLE}_BHB_R2.fasta" -S "${SAMPLE}_BHB.sam"
    samtools sort "${SAMPLE}_BHB.sam" > "${SAMPLE}_BHB_sorted.bam"
    samtools index "${SAMPLE}_BHB_sorted.bam"
    awk '$1~/^@/ || $1~/_E$/ {print}' "${SAMPLE}_BHB.sam" | samtools sort > "${SAMPLE}_BHB_E_sorted.bam"
    samtools index "${SAMPLE}_BHB_E_sorted.bam"
    awk '$1~/^@/ || $1~/_C$/ {print}' "${SAMPLE}_BHB.sam" | samtools sort > "${SAMPLE}_BHB_C_sorted.bam"
    samtools index "${SAMPLE}_BHB_C_sorted.bam"
    rm "${SAMPLE}_BHB.sam"

    popd >/dev/null
done

# -----------------------------------------------------------------------------
# Step 4: set1
# -----------------------------------------------------------------------------
JUNCTION_BEDS=()
for SAMPLE in "${SAMPLES[@]}"; do
    JUNCTION_BEDS+=("${WORK_DIR}/${SAMPLE}/${SAMPLE}_junction.bed")
done
python "${MERGE_JUNCTION_SCRIPT}" --beds "${JUNCTION_BEDS[@]}" -o "${SPECIES}_Set1.bed"

# -----------------------------------------------------------------------------
# Step 5: set2
# -----------------------------------------------------------------------------
python "${GET_SET2_SCRIPT}" --species "${SPECIES}" --set1 "${SPECIES}_Set1.bed" --genome "${GENOME}" --specific "${SPECIFIC}"
