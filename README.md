# CRISPR-TN5
Code for characterizing complex genomic alterations caused by CRISPR gene editing using Illumina Nextera. Code is provided in MATLAB and should run in multiple operating systems but only tested on Windows. 

# Purpose
The code is provided for reproducibility of the published report (not published yet). This is not intended as a research tool but if you use it and have questions, please don't hesitate to ask and I'll do my best to respond. The code analyzes genomic results of dual nucleases for targeted deletion of a region of DNA. The analysis also includes room for two AAV vector genomes to queary for targeted integration, either intentional or unintentional.

# Function
- Read in fastq file in either .fastq or .fastq.gz
- Remove misprimnig events
- Bin reads as one of the following:
1. Intact unedited allele
2. Indel
3. Targeted deletion
4. AAV vector genome integration
5. Targeted inversion
6. Other (mystery) - these were then analyzed manually or with other alignment software. These could include chromosome translocations, large unexpected deletions, and sequencing/polymerase artificats. 

# Experimental  design
1. Nuclease (or CRISPR gRNA) design
2. Primer design
3. Optimization and troubleshooting
4. Analysis

![Figure 1](https://github.com/chrisnelsonlab/CRISPR-TN5/blob/master/images/Nextera_Fig1.png)
Figure 1 - Possible outcomes of multiplexed genome editing

![Figure 2](https://github.com/chrisnelsonlab/CRISPR-TN5/blob/master/images/Nextera_Fig2.png)
Figure 2 - Overview of Nextera enrichment for genome editing analysis

![Figure 3](https://github.com/chrisnelsonlab/CRISPR-TN5/blob/master/images/Nextera_Fig3.png)

## Parameter file
- Line 1 - Amplicon sequence if no editing takes place
- Line 2 - Amplicon sequence if a targeted deletion occurs
- Line 3 - Amplicon sequence if a targeted inversion occurs
- Line 4 - DSB location in bp
- Line 5 - BP window to check for indels on each side of DSB
- Line 6 - Minimum score for aligning 
- Line 7 - Minimum score for aligning short portion
- Line 8 - Minimum alignment score for the AAV genome
- Line 9 - Large portion of target gene with inversion for unexpected modifications (e.g. chewback)
- Line 10 -  Large portion of target gene for unexpected modifications (e.g. chewback)
- Line 11 - minimum amplicon size. Reject shorter amplicons
- Line 12 - AAV Genome 1
- Line 13 - AAV genome 2
- Line 14 - Value for gap opening when aligning to the dystrophin gene
- Line 15 - Value for gap opening when aligning to the AAV genome
- Line 16 - Value for gap extending for AAV genome

## Example Parameter file
> See linked file here

## Example data set
> See linked .fastq files here

## Data to recreate figures from (PAPER) can be found here
> SRA ########
