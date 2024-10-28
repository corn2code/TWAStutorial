# Tutorial TWAS (part 1: get expression data from RNA-Seq data)

The tutorial asumes the student has basic coding knowledge in Bash, Python, and R.

This tutorial is a guide to perform TWAS with the data used in [Torres‐Rodríguez, J. Vladimir, et al. "Population‐level gene expression can repeatedly link genes to functions in maize." The Plant Journal (2024)](https://onlinelibrary.wiley.com/doi/full/10.1111/tpj.16801).

This first part of the tutorial aims to show the process of obtaining the expression data needed to perform TWAS.

## Step 1: Setup
Go to your work directory and create the respective folders to work on

```bash
cd $WORK
mkdir TWASTutorial
mkdir TWASTutorial/scripts
mkdir TWASTutorial/output
mkdir TWASTutorial/input
mkdir TWASTutorial/log.out
cd TWASTutorial
```

## Step 2: Download data
We are downloading raw data from the [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena/browser/home) under the study accession number: PRJEB67964. This data includes RNA-Seq from 750 maize individuals (1,500 files) collected from leaves at the mature stage in Nebraska in 2020.

Download the script "downloadSamples.sh", move it to folder 'scripts' and run it. The code provided will run in the Background

```bash
#Run this command to give the script executable permissions:
chmod +x downloadFiveSamples.sh

#Execute the Script
./downloadFiveSamples.sh > download_log.txt 2>&1

#check if the files exist
ls ../fasta
ls ../fasta -lh # -lh displays more information
ls ../fasta | wc -l # displays the number of files, for this example should be 12 belonging to six individuals (for example 2369_1.fastq.gz & 2369_2.fastq.gz)

#if you want to open the faste file you can use, you can se multiple lines for each
#for more information check: https://en.wikipedia.org/wiki/FASTQ_format
zcat ../fasta/2369_1.fastq.gz | head

```
## Step 3: Check quality of the data
We can use [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [MultiQC](https://github.com/MultiQC/MultiQC)

First, run and example for "2369_1.fastq.gz". Please measure the time with a stopwatch (around 3 min per sample).
```bash
ml fastqc #if the module is not loaded
fastqc ../fasta/2369_1.fastq.gz -o ../fasta/fastqc_reports
```
If we are running few samples we can adjust the code to run like (see below), but the time can increase massively (i.e for 12 samples = 36 min but what if we use the 1,500 samples)
```bash
fastqc ../fasta/*.fastq.gz -o ../fasta/fastqc_reports
```
We will run jobs using slurm files and the array option:

```bash
nano fastqc.slurm # to create a new slurm file
#copy the following code inside

#!/bin/sh
#SBATCH --array=1-12
#SBATCH --job-name=fastqc
#SBATCH --time=3-00:00:00
#SBATCH --mem-per-cpu=10GB
#SBATCH --output=../log.out/%x_%a.out
#SBATCH --partition=schnablelab,batch
#SBATCH --mail-type=ALL

# Load FastQC module
module load fastqc

# Define the list of samples and get the current sample
SAMPLE_FILE="samples.txt"
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLE_FILE)
echo ${SAMPLE}

# Run FastQC on the selected sample
fastqc ${SAMPLE} -o ../fasta/fastqc_reports
```
Then we can use MultiQC to aggregate all the reports

```bash
ml multiqc #if the module is not loaded
multiqc ../fasta/fastqc_reports -o ../fasta/multiqc_report
```
you can submit the job with:

```bash
sbatch fastqc.slurm # to run the file
```

## Step 4: Remove adapters and low quality bases
This step uses [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)

The order is important! Trimming occurs in the order in which the steps are specified on the command line. It is recommended that adapter clipping (ILLUMINACLIP) is done as early as possible.

Since we know now how to parallelize, let's use arrays (but still each job will take ~20 minutes. We suggest to download data from xx and placed as folder "fasta.trimmed")

We can find multiple [adapter](https://github.com/usadellab/Trimmomatic/tree/main/adapters) libraries in the Trimmomatic GitHub page.

Create a file with the location of the fastq files:
```bash
#to create a file with the name of the individuals (six in this example))
find ../fasta -name "*_1.fastq.gz" | sed 's/_1.fastq.gz//g' > files.path.txt
```

This step will take ~20-25 min per sample (I suggest to use the samples provided in the folder to download to save some time but you are welcome to run this anytime)

```bash
# first create this folder:
mkdir ../fasta.trimmed

nano trimmomatic.slurm
#copy pase the following code in side

#!/bin/sh
#SBATCH --array=1-6
#SBATCH --job-name=trimm
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=10GB
#SBATCH --output=../log.out/%x_%a.out
#SBATCH --partition=schnablelab,batch
#SBATCH --mail-type=ALL

ml load trimmomatic/0.33
ml load java/12

samplesheet="files.path.txt"
f=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet |  awk '{print $1}'`
o=`echo ${f} | cut -d'/' -f3-`

mkdir ../fasta.trimmed/${o}

java -jar $TM_HOME/trimmomatic.jar PE -phred33 ${f}_R1_001.fastq.gz ${f}_R2_001.fastq.gz ../fasta.trimmed/${o}/${o}_1_paired.fastq.gz ../fasta.trimmed/${o}/${o}_1_unpaired.fastq.gz output/${o}/${o}_2_paired.fastq.gz ../fasta.trimmed/${o}/${o}_2_unpaired.fastq.gz ILLUMINACLIP:../TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35

```
you can submit the job with:

```bash
sbatch trimmomatic.slurm
```

## Step 5: check quality after removing adapters and low quality bases

We will run jobs using slurm files and the array option:

```bash
mkdir ../fasta.trimmed/fastqc_reports
nano fastqc2.slurm # to create a new slurm file

#the file will include

#!/bin/sh
#SBATCH --array=1-12
#SBATCH --job-name=fastqc2
#SBATCH --time=3-00:00:00
#SBATCH --mem-per-cpu=10GB
#SBATCH --output=../log.out/%x_%a.out
#SBATCH --partition=schnablelab,batch
#SBATCH --mail-type=ALL

# Load FastQC module
module load fastqc

# find ../fasta.trimmed -name "*_paired.fastq.gz" > samples.trimmed.txt

# Define the list of samples and get the current sample
SAMPLE_FILE="samples.trimmed.txt"
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLE_FILE)
echo ${SAMPLE}

# Run FastQC on the selected sample
fastqc ${SAMPLE} -o ../fasta.trimmed/fastqc_reports
```
you can submit the job with:

```bash
sbatch fastqc2.slurm
```

Then we can use MultiQC to aggregate all the reports

```bash
ml multiqc #if the module is not loaded
multiqc ../fasta.trimmed/fastqc_reports -o ../fasta.trimmed/multiqc_report
```

## Step 6: Quantify gene expression
This step will use [Kallisto](https://pachterlab.github.io/kallisto/manual). 

Kallisto’s pseudo-alignment approach makes it exceptionally fast and memory-efficient, suitable for analyzing large RNA-Seq datasets. It’s widely used in gene expression studies due to its accuracy and speed, though it may lack the detailed read-level information that traditional aligners provide.

First we need to create an index of the transcript. Download the file "Zmays_833_Zm-B73-REFERENCE-NAM-5.0.55.transcript_primaryTranscriptOnly.fa" and load it to your working directory. In the folder where you placed this file run: (In this example I am placing it in the folder "input")

```bash
ml kallisto/0.46
kallisto index -i Zm-B73-REFERENCE-NAM-5.0.55.transcript_primaryTranscriptOnly.idx Zmays_833_Zm-B73-REFERENCE-NAM-5.0.55.transcript_primaryTranscriptOnly.fa.gz
```

We then can map all the reads from "*paired.gz" files.

```bash
# Ensure the output directory exists:
mkdir ../input/out.kallisto

nano kallisto.slurm
#copy paste the following code

#!/bin/sh
#SBATCH --array=1-6
#SBATCH --job-name=quant
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=10GB
#SBATCH --output=../log.out/%x_%a.out
#SBATCH --partition=schnablelab,batch
#SBATCH --mail-type=ALL

ml load kallisto/0.46

samplesheet="files.path.txt"
f=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet |  awk '{print $1}'`
o=`echo ${f} | cut -d'/' -f3-`

kallisto quant --threads=10 -i ../input/Zm-B73-REFERENCE-NAM-5.0.55.transcript_primaryTranscriptOnly.idx -o ../input/out.kallisto/${o} ../fasta.trimmed/${o}/${o}_1_paired.fastq.gz ../fasta.trimmed/${o}/${o}_2_paired.fastq.gz
```

This will create a folder for each individual with the following files:

1. **abundance.tsv**: This is the primary output file, which contains transcript-level expression estimates. It includes the following columns:
   - `target_id`: The transcript ID (usually matching the transcript IDs in the reference transcriptome).
   - `length`: The length of each transcript.
   - `eff_length`: The effective length of the transcript, adjusted for read length.
   - `est_counts`: The estimated number of reads assigned to each transcript.
   - `tpm`: Transcripts per million, a normalized measure of expression.

2. **abundance.h5**: A binary file containing all abundance information in HDF5 format. This file is used by downstream tools, such as Sleuth, to perform differential expression analysis and provides more efficient data storage and access than a plain text file.

3. **run_info.json**: A JSON file containing metadata about the run, including the Kallisto version, the command used, the number of processed reads, and the total runtime. This file is useful for tracking parameters and ensuring reproducibility.


## Step 7: final gene expression table

Create a new pyhton file:
```bash
nano make_rna_ss.py
```
paste inside the code:

```python
import os
import sys

mydir = sys.argv[1]
if not os.path.exists(mydir):
    sys.exit("{0} is not a valid directory".format(mydir))

gene_exp_dict = {}
sample_list = []
mysamples = os.listdir(mydir)
for asample in mysamples:
    if not os.path.isdir(mydir + "/" + asample): continue
    if not os.path.exists(mydir + "/" + asample + "/" + "abundance.tsv"): continue
    fh = open(mydir + "/" + asample + "/" + "abundance.tsv")
    sample_list.append(asample)
    fh.readline()
    for x in fh:
        y = x.strip().split('\t')
        mygene = y[0]
        mytpm = float(y[-1])
        if not mygene in gene_exp_dict: gene_exp_dict[mygene] = {}
        gene_exp_dict[mygene][asample] = mytpm
    fh.close()
#print(sample_list)
#print(gene_exp_dict[list(gene_exp_dict)[0]])
fh = open("merged_gene_tpms.csv",'w') #change with the desire name for output
myheader = ["GeneID"] + sorted(sample_list)
fh.write(",".join(myheader)+"\n")
for agene in sorted(list(gene_exp_dict)):
    plist = [agene]
    for asample in sorted(sample_list):
        plist.append(gene_exp_dict[agene][asample])
    fh.write(",".join(map(str,plist))+"\n")
```

Run it:

```bash
python3 make_rna_ss.py ../input/out.kallisto
```
...




