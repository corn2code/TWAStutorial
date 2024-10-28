# Tutorial for Transcriptome-Wide Association Study (TWAS)
**Part 1: Obtaining Gene Expression Levels from RNA-Seq Data**

This tutorial guides you through obtaining gene expression data to perform TWAS. It is based on data used in the paper [Torres‐Rodríguez, J. Vladimir, et al., 2024, "Population‐level gene expression can repeatedly link genes to functions in maize"](https://onlinelibrary.wiley.com/doi/full/10.1111/tpj.16801).

---

## Prerequisites
The tutorial assumes familiarity with **Bash**, **Python**, and **R**. You will need access to `fastqc`, `multiqc`, `trimmomatic`, and `kallisto`.  
Large datasets, including trimmed FASTA files and the `.idx` file for Kallisto, are available on Figshare: [Figshare Link](https://figshare.com/articles/dataset/TWAS_tutorial/27312822).

---

## Step 1: Setup
Create and navigate to the working directory:

```bash
cd $WORK
mkdir -p TWASTutorial/{scripts,output,input,log.out}
cd TWASTutorial
```

## Step 2: Download Data
We will download a subset of RNA-Seq data from the [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena/browser/home) under study accession **PRJEB67964**, which contains RNA-Seq data from 750 maize samples. For this tutorial, we’ll download only 12 samples from six individuals.

1. Download `downloadSamples.sh` to the `scripts` folder.
2. Make it executable and run it in the background.

```bash
chmod +x scripts/downloadSamples.sh
./scripts/downloadSamples.sh > download_log.txt 2>&1

# Check file count
ls ../fasta | wc -l  # Should be 12 (six paired samples)
zcat ../fasta/2369_1.fastq.gz | head  # View file structure
```

## Step 3: Quality Check with FastQC and MultiQC
We’ll use [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [MultiQC](https://github.com/MultiQC/MultiQC) for data quality assessment.

### Run FastQC on a single sample
```bash
ml fastqc/0.12
fastqc ../fasta/2369_1.fastq.gz -o ../fasta/fastqc_reports
```

### Run FastQC on multiple samples using SLURM
Create a `fastqc.slurm` file with the following code:

```bash
#!/bin/sh
#SBATCH --array=1-12
#SBATCH --job-name=fastqc
#SBATCH --time=3-00:00:00
#SBATCH --mem-per-cpu=10GB
#SBATCH --output=../log.out/%x_%a.out
#SBATCH --partition=schnablelab,batch
#SBATCH --mail-type=ALL

ml fastqc/0.12

SAMPLE_FILE="samples.txt"
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLE_FILE)
fastqc ${SAMPLE} -o ../fasta/fastqc_reports
```

Submit the job:
```bash
sbatch fastqc.slurm
```

Run MultiQC to aggregate FastQC reports:
```bash
ml multiqc
multiqc ../fasta/fastqc_reports -o ../fasta/multiqc_report
```

## Step 4: Adapter and Quality Trimming with Trimmomatic
Use [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) to remove adapters and low-quality bases. Trimming order matters, so perform `ILLUMINACLIP` first.

1. Generate a file with fastq paths:
   ```bash
   find ../fasta -name "*_1.fastq.gz" | sed 's/_1.fastq.gz//g' > files.path.txt
   ```

2. Create the SLURM script `trimmomatic.slurm`:
   ```bash
   #!/bin/sh
   #SBATCH --array=1-6
   #SBATCH --job-name=trimm
   #SBATCH --time=1-00:00:00
   #SBATCH --mem-per-cpu=10GB
   #SBATCH --output=../log.out/%x_%a.out
   #SBATCH --partition=schnablelab,batch
   #SBATCH --mail-type=ALL

   ml trimmomatic/0.33 java/12

   samplesheet="files.path.txt"
   f=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}'`
   o=`echo ${f} | cut -d'/' -f3-`

   mkdir ../fasta.trimmed/${o}

   java -jar $TM_HOME/trimmomatic.jar PE -phred33 \
   ${f}_1.fastq.gz ${f}_2.fastq.gz \
   ../fasta.trimmed/${o}/${o}_1_paired.fastq.gz ../fasta.trimmed/${o}/${o}_1_unpaired.fastq.gz \
   ../fasta.trimmed/${o}/${o}_2_paired.fastq.gz ../fasta.trimmed/${o}/${o}_2_unpaired.fastq.gz \
   ILLUMINACLIP:../TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35
   ```

Submit the job:
```bash
sbatch trimmomatic.slurm
```

## Step 5: Quality Check After Trimming
Repeat FastQC and MultiQC on trimmed files to ensure quality.

```bash
mkdir ../fasta.trimmed/fastqc_reports

nano fastqc2.slurm  # Create a new SLURM file for trimmed samples
```

SLURM script `fastqc2.slurm`:
```bash
#!/bin/sh
#SBATCH --array=1-12
#SBATCH --job-name=fastqc2
#SBATCH --time=3-00:00:00
#SBATCH --mem-per-cpu=10GB
#SBATCH --output=../log.out/%x_%a.out
#SBATCH --partition=schnablelab,batch
#SBATCH --mail-type=ALL

ml fastqc
SAMPLE_FILE="samples.trimmed.txt"
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLE_FILE)
fastqc ${SAMPLE} -o ../fasta.trimmed/fastqc_reports
```

Submit the job:
```bash
sbatch fastqc2.slurm
```

Aggregate reports with MultiQC:
```bash
ml multiqc
multiqc ../fasta.trimmed/fastqc_reports -o ../fasta.trimmed/multiqc_report
```

## Step 6: Quantify Gene Expression with Kallisto
Kallisto’s pseudo-alignment method provides accurate gene quantification.

1. Create the Kallisto index:
   ```bash
   ml kallisto/0.46
   kallisto index -i Zm-B73.idx Zmays_transcript.fa.gz
   ```

2. Map paired reads with Kallisto:
   ```bash
   mkdir ../input/out.kallisto

   nano kallisto.slurm  # SLURM script for Kallisto quantification
   ```

SLURM script `kallisto.slurm`:
```bash
#!/bin/sh
#SBATCH --array=1-6
#SBATCH --job-name=quant
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=10GB
#SBATCH --output=../log.out/%x_%a.out
#SBATCH --partition=schnablelab,batch
#SBATCH --mail-type=ALL

ml kallisto/0.46

samplesheet="files.path.txt"
f=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}'`
o=`echo ${f} | cut -d'/' -f3-`

kallisto quant --threads=10 -i ../input/Zm-B73.idx -o ../input/out.kallisto/${o} \
 ../fasta.trimmed/${o}/${o}_1_paired.fastq.gz ../fasta.trimmed/${o}/${o}_2_paired.fastq.gz
```

Submit the job:
```bash
sbatch kallisto.slurm
```

## Step 7: Generate the Final Gene Expression Table
The following Python script merges all `abundance.tsv` files into a single gene expression matrix.

1. Save as `make_rna_ss.py`:
   ```python
   import os
   import sys

   mydir = sys.argv[1]
   if not os.path.exists(mydir):
       sys.exit("{0} is not a valid directory".format(mydir))

   gene_exp_dict = {}
   sample_list = []
   for asample in os.listdir(mydir):
       if not os.path.exists(f"{mydir}/{asample}/abundance.tsv"): continue
       with open(f"{mydir}/{asample}/abundance.tsv") as fh:
           sample_list.append(asample)
           fh.readline()  # Skip header
           for line in fh:
               gene, _, _, tpm = line.strip().split()
               if gene not in gene_exp_dict: gene_exp_dict[gene] = {}
               gene_exp_dict[gene][asample] = float(tpm)

   with open("merged_gene_tpms.csv", 'w') as outfile:
       header = "Gene," + ",".join(sample_list)
       print(header, file=outfile)
       for gene, expdict in gene_exp_dict.items():
           exps = [str(expdict.get(x, 0.0)) for x in sample_list]
           print(gene + "," + ",".join(exps), file=outfile)
   ```

2. Run the script:
   ```bash
   python3 make_rna_ss.py ../input/out.kallisto
   ```

The resulting file `merged_gene_tpms.csv` will contain the TPM values for each gene across all samples.
