# Assessment_2

**Author**

200394457


**Process**


Practical8 run using R studio (version 4.3.1).

**Data** 


The GSE116583 dataset: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116583.

Download reference mouse transcriptome M31: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M31/gencode.vM31.transcripts.fa.gz.

Download quant.sf files: https://github.com/sjcockell/mmb8052/raw/main/practicals/practical_08/results/counts.zip.

Download the count data: https://github.com/sjcockell/mmb8052/raw/main/practicals/practical_08/results/counts.zip.

**Software** 


FASTQ data were downloaded (September 2023) with sra-tools (fastq-dump v2.8.0).

`Salmon` (v1.9.0) used for quantification. 

Statistical analysis was performed in R studio (version 4.3.1).

R studio add-on packages used from `bioconductor` (BiocManager 1.30.19): `tximport`, `DESeq2`, `biomaRt`, `pheatmap`, and `tidyverse`.
