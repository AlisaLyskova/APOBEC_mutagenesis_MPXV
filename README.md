# APOBEC mutagenesis in MPXV

We analyzed activity of APOBEC3 family proteins in long-read Oxford Nanopore and short-read ILLUMINA transcriptomic and genomic RNA-seq data with MPXV.

## SCRIPTS

We downloaded samples (sratoolkit-3.1.1, fastq-dump) from projects PRJEB56841, PRJEB60728, PRJNA1183318, PRJNA906618, PRJNA980137, PRJNA845087, PRJNA981509 (see data/samples_description.csv).  

In order to correctly determine the largest number of virus reads, we considered three options for the alignment of the RNA-seq data and performed the procedure for generating reads using InSilicoSeq-2.0.1 (https://github.com/HadrienG/InSilicoSeq), see scripts/iss_monkeys.sh and scripts/iss_human.sh.

We created pipeline for processing and analyzing RNA-seq data (see pipeline.sh). The pipeline consists of two parts. First part is processing RNA-seq data including fastq files quality control (FASTQC v0.12.1), alignment (minimap2-2.28-r1209 and STAR_2.7.11b), alignment quality control (mosdepth-0.3.9 and samtools-1.20 stats), variant calling (clair3-v1.0.10), variants quality control (using VAF distribution). Second part is searching for APOBEC substitutions in samples and analyzing activity of APOBEC3 family proteins.

Second part of the pipeline contains the following scripts:
- Creating pictures for found SNPs in samples (heatmap, LOGO, number of reads for 3 nucleotide context of substitutions)
- Filtering APOBEC subsitutions
- Random modulation procedure (to access the distribution of detected APOBEC substitutions)
- Translating substitutions C>T and G>A to the coding chain (file with annotation is required, see data/GCF_014621545.1_ASM1462154v1_genomic.gff); to additionally assess the presence of antisense transcription, files data/TSS.xlsx and data/TES.xlsx are required (here from https://journals.asm.org/doi/10.1128/msphere.00356-24)
- Analyzing the dependence of the positions of substitutions C>T and G>A and genomic structures (circos-0.69-8, files with type of genes (see data/ucsc_early_late.txt), with repeats (see data/repeats_MPXV.tsv) and with annotation are required)
- Analyzing the secondary structure in the detected positions (RNAsselem, https://github.com/KazanovLab/RNAsselem)
- Comparison of the nucleotide context of substitutions C>T and G>A with known specificity of proteins of the APOBEC family