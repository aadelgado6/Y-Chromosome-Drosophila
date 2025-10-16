#!/bin/bash
#SBATCH --job-name=Y_filtered_extract
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH -o y_inclusive.%j.out
#SBATCH -e y_inclusive.%j.err

# ---------------------------------------------------------
# Competitive + filtered Y-read extraction pipeline (uses combinedISO1.fasta)
# - maps HiFi reads (reads.fastq) to combinedISO1.fasta (whole-genome + new Y)
# - selects reads that have ANY alignment to a Y contig that meets MAPQ and alignment-length filters
# - deduplicates, extracts FASTQ and runs hifiasm
# ---------------------------------------------------------

set -euo pipefail

THREADS=32

# ---------- parameters (tweakable) ----------
MAPQ_MIN=20          # minimum mapping quality to accept an alignment
MIN_REF_ALN=1000     # minimum alignment length on reference (bp)
MIN_PROP=0.80        # OR accept alignments that cover >= MIN_PROP of the read

# ---------- inputs (edit if needed) ----------
# combinedISO1.fasta should be your concatenated whole-genome + new Y file
COMBINED_FA="combinedISO1.fasta"
Y_FA="newiso1Y.fasta"        # used to extract the list of Y contig names internally by awk
ALL_READS="reads.fastq"      # can be .fastq or .fastq.gz

# ---------- outputs (all inside OUTDIR) ----------
OUTDIR="TargetedYAssembly"
mkdir -p "${OUTDIR}"

MINIMAP_INDEX="${COMBINED_FA%.fasta}.mmi"
ALL_TO_COMBINED_BAM="${OUTDIR}/all_to_combined.bam"
Y_FILTERED_NAMES="${OUTDIR}/y_filtered_read_names.txt"
Y_FILTERED_FASTQ="${OUTDIR}/y_reads_filtered.fastq"
HIFIASM_PREFIX="${OUTDIR}/targeted_hifiasmY"
MINIMAP_SUMMARY="${OUTDIR}/minimap_combined.flagstat.txt"
SUMMARY="${OUTDIR}/y_filtered_summary.txt"
HIFIASM_LOG="${OUTDIR}/hifiasm_targeted_hifiasmY.log"

# ---------- modules (edit to your cluster names if needed) ----------
module load GCC/13.2.0 minimap2/2.29 SAMtools/1.21
module load seqtk/1.4
module load hifiasm/0.24.0

echo "Starting filtered Y-read extraction: $(date)"
echo "THREADS=${THREADS} MAPQ_MIN=${MAPQ_MIN} MIN_REF_ALN=${MIN_REF_ALN} MIN_PROP=${MIN_PROP}"
echo "Combined reference: ${COMBINED_FA}"
echo "Y contig source:    ${Y_FA}"

# ---------- 0) sanity checks ----------
if [ ! -f "${COMBINED_FA}" ]; then
  echo "ERROR: combined reference '${COMBINED_FA}' not found. Exiting."
  exit 1
fi
if [ ! -f "${Y_FA}" ]; then
  echo "ERROR: Y fasta '${Y_FA}' not found. Exiting."
  exit 1
fi

# ---------- 1) index combined reference for minimap2 (build if missing) ----------
if [ ! -f "${MINIMAP_INDEX}" ]; then
  echo "Creating minimap2 index: ${MINIMAP_INDEX}"
  minimap2 -d "${MINIMAP_INDEX}" "${COMBINED_FA}"
else
  echo "Found existing minimap2 index: ${MINIMAP_INDEX}"
fi

# ---------- 2) competitive mapping: map all reads to combined ref ----------
minimap2 -t "${THREADS}" -ax map-hifi "${MINIMAP_INDEX}" "${ALL_READS}" \
  | samtools view -b -o "${ALL_TO_COMBINED_BAM}" -

samtools flagstat "${ALL_TO_COMBINED_BAM}" > "${MINIMAP_SUMMARY}" || true
echo "Mapping summary written to ${MINIMAP_SUMMARY}"

# ---------------------------
# 3+4) Select reads with qualifying alignments to Y contigs
#      (reads Y contig names directly from ${Y_FA}; no external list file)
# ---------------------------
echo "Selecting reads with qualifying alignments to Y contigs (reading names from ${Y_FA})..."

samtools view -h "${ALL_TO_COMBINED_BAM}" \
  | awk -v MAPQ_MIN="${MAPQ_MIN}" -v MIN_REF_ALN="${MIN_REF_ALN}" -v MIN_PROP="${MIN_PROP}" -v YFA="${Y_FA}" '
    BEGIN{
      FS = "\t"
      # load Y contig names directly from the Y fasta (first word of each header)
      while ((getline line < YFA) > 0) {
        if (line ~ /^>/) {
          sub(/^>/,"",line)
          split(line, a, /[ \t]/)
          y[a[1]] = 1
        }
      }
      close(YFA)
    }
    /^@/ { next }   # skip SAM header lines
    {
      qname = $1
      rname = $3
      mapq = $5 + 0
      cigar = $6
      seq = $10
      seq_len = (seq=="" ? 0 : length(seq))

      # compute reference-aligned length and read-aligned length from CIGAR
      ref_aln_len = 0
      read_aln_len = 0
      s = cigar
      while (match(s, /([0-9]+)([MIDNSHP=X])/, a)) {
        len = a[1] + 0
        op = a[2]
        if (op == "M" || op == "=" || op == "X" || op == "D" || op == "N") { ref_aln_len += len }
        if (op == "M" || op == "=" || op == "X" || op == "I" ) { read_aln_len += len }
        s = substr(s, RSTART + RLENGTH)
      }

      # only consider alignments to Y contigs loaded from YFA
      if (!(rname in y)) next
      if (mapq < MAPQ_MIN) next

      if ( (ref_aln_len >= MIN_REF_ALN) || ( seq_len > 0 && read_aln_len >= (MIN_PROP * seq_len) ) ) {
        print qname
      }
    }
' \
  | sort -u > "${Y_FILTERED_NAMES}"

echo "Selected $(wc -l < "${Y_FILTERED_NAMES}") unique read names that meet filters."

# ---------- 5) extract reads from original FASTQ by name ----------
if [ ! -s "${Y_FILTERED_NAMES}" ]; then
  echo "ERROR: No reads passed the filtering criteria. Check thresholds or input files."
  exit 1
fi

seqtk subseq "${ALL_READS}" "${Y_FILTERED_NAMES}" > "${Y_FILTERED_FASTQ}"

# ---------- 6) quick checks: read count & total bases ----------
READ_COUNT=$(awk 'END{print NR/4}' "${Y_FILTERED_FASTQ}")
BASES=$(awk 'NR%4==2{b+=length($0)} END{print b+0}' "${Y_FILTERED_FASTQ}")

echo "Filtered reads extracted: ${READ_COUNT}"
echo "Filtered bases total:    ${BASES}"

# ---------- 7) save a small run summary ----------
cat > "${SUMMARY}" <<EOF
Date:               $(date)
Combined reference:  ${COMBINED_FA}
Y assembly:          ${Y_FA}
All reads file:      ${ALL_READS}
All-to-combined BAM: ${ALL_TO_COMBINED_BAM}
Filtered read names: ${Y_FILTERED_NAMES}
Filtered FASTQ:      ${Y_FILTERED_FASTQ}
Read count:          ${READ_COUNT}
Total bases:         ${BASES}
MAPQ threshold:      ${MAPQ_MIN}
Min ref alignment:   ${MIN_REF_ALN}
Min proportion:      ${MIN_PROP}
Hifiasm prefix:      ${HIFIASM_PREFIX}
Threads:             ${THREADS}
Slurm job id:        ${SLURM_JOB_ID:-NA}
EOF

echo "--- summary written to ${SUMMARY} ---"

# ---------- 8) run Hifiasm on filtered reads ----------
if [ ! -s "${Y_FILTERED_FASTQ}" ] ; then
  echo "ERROR: filtered FASTQ is empty. Exiting."
  exit 1
fi

hifiasm -o "${HIFIASM_PREFIX}" -t "${THREADS}" "${Y_FILTERED_FASTQ}" > "${HIFIASM_LOG}" 2>&1

echo "Hifiasm finished. Outputs with prefix: ${HIFIASM_PREFIX}"
ls -lh "${OUTDIR}"/targeted_hifiasmY* || true

echo "Filtered Y-read extraction + assembly completed: $(date)"

