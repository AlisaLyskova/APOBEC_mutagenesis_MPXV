# APOBEC mutagenesis in MPXV

We analyzed activity of APOBEC3 family proteins in long-read Oxford Nanopore and ILLUMINA RNA-seq data with MPXV.

## SCRIPTS

We downloaded samples (sratoolkit-3.1.1, fastq-dump), indexed genome, aligned reads to the genome and found SNPs using bwa-0.7.17 and freebayes-1.3.6 for ILLUMINA samples (scripts/illumina_samples_process.sh) and minimap2-2.28-r1209 and clair3-v1.0.10 for long-read Oxford Nanopore samples (scripts/nanopore_samples_process.sh).

We analyzed context of APOBEC3-like mutations. We annotated these positions and build genome map using circos-0.69-8.

We approved APOBEC3 mutagenesis on MPXV by performing 100 generations. 

Program description: We found potential targets with 1 of 4 motifs (TCT, TCG, TCA, TCC) in the genome. For each sample we created dataframe with read ID and potential targets positions. For each row of this table we randomly placed 1 (mutation) or 0 (no mutation) with the total number of 1 was equal to the total number of reads with mutations in this motif.