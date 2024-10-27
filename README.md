# TWAStutorial
This tutorial is a guide to perform TWAS with the data used in [Torres‐Rodríguez, J. Vladimir, et al. "Population‐level gene expression can repeatedly link genes to functions in maize." The Plant Journal (2024)](https://onlinelibrary.wiley.com/doi/full/10.1111/tpj.16801).

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
We are downloading raw data from the [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena/browser/home) under the study accession number: PRJEB67964. This data includes RNA-Seq from 750 maize individuals collected from leaf at mature stage in Nebraska 2020.

Download the script "downloadFiveSamples.sh", move it to folder 'scripts' and run it. The code provided will run in the Background

```bash
#Run this command to give the script executable permissions:
chmod +x downloadFiveSamples.sh

#Execute the Script
./downloadFiveSamples.sh > download_log.txt 2>&1
```
...

