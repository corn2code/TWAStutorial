#!/usr/bin/env bash
set -euo pipefail

# Download 12 FASTQ files (6 paired-end samples) from ENA.
# Run this script from the `scripts/` directory as documented in README.md.

OUT_DIR="../fasta"
mkdir -p "${OUT_DIR}"

URLS=(
  "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR121/ERR12177542/DK3IBZ2_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR121/ERR12177187/38-11_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR121/ERR12177172/2MA22_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR121/ERR12177187/38-11_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR121/ERR12177171/2369_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR121/ERR12177542/DK3IBZ2_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR121/ERR12177184/793_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR121/ERR12177181/33-16_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR121/ERR12177172/2MA22_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR121/ERR12177181/33-16_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR121/ERR12177171/2369_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR121/ERR12177184/793_2.fastq.gz"
)

for url in "${URLS[@]}"; do
  wget -nc -P "${OUT_DIR}" "${url}"
done

echo "Download complete. FASTQ files in ${OUT_DIR}"
