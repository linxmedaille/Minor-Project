

# CS39620: Stop Codon Frequency

### Author: Tom Janvier (230128875)
### Scheme: Data Science 7G73
### Date: 09/02/2026

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
