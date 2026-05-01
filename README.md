

# CS39620: Stop Codon Frequency

### Author: Tom Janvier
### Scheme: Data Science 

## Project Description

DNA is converted into protein using a triplet code, where each group of three nucleotides (a codon) represents an amino acid. Three specific codons TAG, TAA, and TGA act as stop codons that do not code any amino acid. DNA stores the genetic information of the cell, and within genes, stop codons mark the end of a protein.

The aim of this project is to analyse stop codon usage across a large number of bacterial genomes using data from the AllTheBacteria database. By extracting coding sequences from these genomes, the project will measure the proportions of each stop codon.

The main goal is to build a computational system that can process large genomic datasets, identify unusual stop codon usage, and produce visualisations that highlight differences between groups of bacteria. I will then use my analysis to check if any genomes are unusual within their species or genus and produce suitable visualisations of the results.

## Proposed Tasks

The first task will be to carry out background research into genetic coding, stop codons, and codon usage bias. This will involve reading introductory material on genomics and reviewing existing studies on stop codon usage in bacteria in order to understand what patterns are expected.

The second task will be to learn how to access and process genome data from the AllTheBacteria database. This will include understanding how coding regions are represented, and extracting relevant sequence data.

The third task will be to design and implement software to analyse the genomes. This will include writing code (likely in Python)

Parse genome files, Identify coding sequences, Count occurrences of each stop codon, Calculate stop codon proportions for each genome

The fourth task will be to perform comparative analysis across species and genera. This will involve grouping genomes taxonomically and identifying outliers or unusual cases.

The final task will be to produce visualisations and interpret results. This will potentially involve generating plots such as bar charts, box plots, or scatter plots to show how stop codon usage varies across different groups of bacteria.

## Project Deliverables

The main deliverables for this project are expected to include:

-A working software for analysing and counting stop codon usage in bacterial genomes.

-A processed dataset containing stop codon counts and proportions for all analysed genomes.

-Visualisations showing variation in stop codon usage across species and genera.

-A short literature review on stop codon usage and codon bias.

-Documentation explaining system design, data sources, and methods used.

-A final project report presenting methodology, results, analysis, and conclusions.

## Bibliography

“AllTheBacteria Documentation — AllTheBacteria Documentation.” Readthedocs.io, 2024, allthebacteria.readthedocs.io/en/latest/. Accessed 9 Feb. 2026.

Reddy, Michael K. “Amino Acid | Definition, Structure, & Facts.” Encyclopædia Britannica, 11 Jan. 2024, www.britannica.com/science/amino-acid
.

“What Are the 3 Stop Codons and How Do They Work?” Biology Insights, 9 Jan. 2026, biologyinsights.com/what-are-the-3-stop-codons-and-how-do-they-work/. Accessed 9 Feb. 2026.

“Documentation · Biopython.” Biopython.org, biopython.org/wiki/Documentation. Accessed 12 Feb. 2026.




Stop Codon Frequency Analyser
A Python tool for analysing terminal stop codon usage (TAA, TAG, TGA) across bacterial genomes, using FASTA assembly files and Bakta JSON annotation files.

Requirements

Python 3.10+
Biopython — pip install biopython
Requests — pip install requests


Project Structure
project/
├── analysis.py            # Core stop codon extraction logic
├── batch_pipeline.py      # Large-scale download and processing pipeline
├── DataFiles/
│   ├── FAFiles/           # FASTA genome files (.fa)
│   ├── JsonFiles/         # Bakta annotation files (.bakta.json)
│   └── Lists/
│       ├── file_list.r0.2.v2.tsv        # AllTheBacteria FA manifest
│       └── atb.bakta.r0.2.status.tsv    # AllTheBacteria Bakta manifest
├── stop_codon_results.csv         # Output: per-genome stop codon counts
├── codon_usage_results.csv        # Output: full codon usage table
└── pipeline_progress.json         # Auto-generated: pipeline resume file

Quick Start — Single Genome (analysis.py)
Use this if you have a FASTA file and its Bakta JSON annotation and just want to analyse one genome locally.
1. Place your files in the correct folders:
DataFiles/FAFiles/myGenome.fa
DataFiles/JsonFiles/myGenome.bakta.json
The filename stem (the part before the first .) must match between the two files.
2. Enable test mode in analysis.py:
Open analysis.py and set:
pythonUSE_TEST_FILES = True
Then change "testFile1" in match_genome_files() to match your filename stem, or set USE_TEST_FILES = False to process all matched files in the folders.
3. Run:
bashpython analysis.py
4. Check outputs:

stop_codon_results.csv — stop codon counts and proportions per genome
codon_usage_results.csv — full codon frequency table per genome


Large-Scale Pipeline — AllTheBacteria (batch_pipeline.py)
Use this to process thousands of genomes from the AllTheBacteria dataset. The pipeline downloads one batch archive at a time, analyses each genome, writes results, and deletes temporary files before moving to the next batch.
1. Download the manifest files from AllTheBacteria and place them at:
DataFiles/Lists/file_list.r0.2.v2.tsv
DataFiles/Lists/atb.bakta.r0.2.status.tsv
2. Run the pipeline:
bashpython batch_pipeline.py
The pipeline will:

Match samples that appear in both manifest files with a Bakta status of PASS
Discover Bakta archive download URLs automatically from the AllTheBacteria OSF project
Download, extract, analyse, and clean up one batch at a time
Append results to stop_codon_results.csv as each genome is processed
Save progress to pipeline_progress.json after every sample

3. Resuming an interrupted run:
Simply run python batch_pipeline.py again. Already-processed samples listed in pipeline_progress.json are skipped automatically.

Output Files
stop_codon_results.csv
One row per genome with the following columns:
ColumnDescriptionGenomeSample IDSpeciesSpecies assignment (batch pipeline only)TAARaw count of TAA stop codonsTAGRaw count of TAG stop codonsTGARaw count of TGA stop codonsTotal_Valid_Stop_CodonsTotal of TAA + TAG + TGACDS_CountNumber of CDS features processedInvalid_CountCDS features where stop codon could not be extractedTAA_ProportionTAA / Total_ValidTAG_ProportionTAG / Total_ValidTGA_ProportionTGA / Total_ValidAvg_CDS_LengthMean CDS length in base pairsGC_ContentGenomic GC content (%)Genome_Size_bpTotal genome size in base pairsContig_CountNumber of contigs in the assembly
codon_usage_results.csv
One row per genome, one column per unique codon triplet observed across the dataset. Values are raw counts of each codon across all CDS features in that genome.

Configuration
Key settings at the top of analysis.py:
VariableDefaultDescriptionFA_DIRDataFiles/FAFilesFolder containing FASTA filesJSON_DIRDataFiles/JsonFilesFolder containing Bakta JSON filesOUT_FILEstop_codon_results.csvStop codon output pathCODON_OUT_FILEcodon_usage_results.csvCodon usage output pathUSE_TEST_FILESTrueIf True, only processes testFile1

Notes

File pairing is done by filename stem: SAMEA1234.fa is matched with SAMEA1234.bakta.json.
Bakta uses 1-based inclusive coordinates. The code handles conversion to Python's 0-based indexing automatically.
On the positive strand, the stop codon is the last 3 bases of the annotated CDS region.
On the negative strand, the stop codon is the reverse complement of the first 3 bases of the annotated region.
Any extracted triplet that is not TAA, TAG, or TGA is recorded as Invalid_Count and excluded from proportions.
