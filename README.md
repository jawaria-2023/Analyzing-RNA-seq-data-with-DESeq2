# Analyzing-RNA-seq-data-with-DESeq2
Analyzing count data from RNA-seq involves a fundamental task: identifying differentially expressed genes. The data is structured as a table, detailing the number of sequence fragments assigned to each gene for every sample. Similar data formats emerge in various assay types such as comparative ChIP-Seq, HiC, shRNA screening, and mass spectrometry. A key analytical objective is to quantify and statistically infer systematic changes between conditions, relative to within-condition variability.

DESeq2 is a package that addresses this challenge, employing negative binomial generalized linear models to test for differential expression. The estimates of dispersion and logarithmic fold changes incorporate data-driven prior distributions. This vignette serves to elucidate the usage of the DESeq2 package, providing insights into typical workflows. While a Bioconductor website RNA-seq workflow covers similar content at a more deliberate pace, including count matrix generation from FASTQ files, this document specifically focuses on DESeq2 package.

RNA-seq Analysis in R
Differential Expression Analysis
30 November 2020
Authors: Belinda Phipson, Anna Trigos, Matt Ritchie, Shian Su, Maria Doyle, Harriet Dashnow, Charity Law

Data files are available from: https://figshare.com/s/1d788fd384d33e913a2a. You should download the files listed below and place them into a folder called data in your working directory.

Data files:

GSE60450_Lactation-GenewiseCounts.txt
SampleInfo.txt
SampleInfo_Corrected.txt
Data files were originally obtained from:
ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60450/suppl/GSE60450_Lactation-GenewiseCounts.txt.gz
http://bioinf.wehi.edu.au/MSigDB/v7.1/Mm.c2.all.v7.1.entrez.rds
http://bioinf.wehi.edu.au/MSigDB/v7.1/Mm.h.all.v7.1.entrez.rds

This material has been inspired by the following resources:
http://www.statsci.org/smyth/pubs/QLedgeRPreprint.pdf (Lun, Chen, and Smyth 2016)
http://monashbioinformaticsplatform.github.io/RNAseq-DE-analysis-with-R/RNAseq_DE_analysis_with_R.html
Inspired by:

QLedgeRPreprint.pdf (Lun, Chen, and Smyth 2016)
RNAseq_DE_analysis_with_R
Bioconductor/R Packages
Packages used:

limma
edgeR
Glimma
org.Mm.eg.db
gplots
RColorBrewer
To install the packages, use the following commands in R:

R
Copy code
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR", "Glimma", "org.Mm.eg.db", "gplots", "RColorBrewer", "NMF", "BiasedUrn"))

____________Overview____________

Reading in Table of Counts
Adding Annotation
Filtering Lowly Expressed Genes
Quality Control
Normalization for Composition Bias
Differential Expression Analysis
Testing Relative to a Threshold
Visualization
Gene Set Testing
Introduction and Data Import
Measuring gene expression on a genome-wide scale has become common practice over the last two decades or so, with microarrays predominantly used pre-2008. With the advent of next-generation sequencing technology in 2008, an increasing number of scientists use this technology to measure and understand changes in gene expression in often complex systems. As sequencing costs have decreased, using RNA-Seq to simultaneously measure the expression of tens of thousands of genes for multiple samples has never been easier. The cost of these experiments has now moved from generating the data to storing and analyzing it.

There are many steps involved in analyzing an RNA-Seq experiment. Analyzing an RNAseq experiment begins with sequencing reads. These are aligned to a reference genome, then the number of reads mapped to each gene can be counted. This results in a table of counts, which is what we perform statistical analyses on in R. While mapping and counting are important and necessary tasks, today we will be starting from the count data and getting stuck into analysis.
