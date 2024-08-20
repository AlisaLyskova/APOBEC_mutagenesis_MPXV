# APOBEC mutagenesis in MPXV

Programs for processing samples from PRJEB56841 project

## SCRIPTS:

1. RNAseq_MPXV_downl_process.sh - download samples, index genome, map reads from samples to the genome using minimap2 and parameters from the article DOI: 10.1038/s41597-023-02149-4
    
    OUTPUT: BAM files and indexed BAM files 
    
2. clair3.sh - search for SNP in BAM files, program description: https://github.com/HKU-BAL/Clair3

    OUTPUT: VCF files

3. Visualization of SNP from vcf files
    
    3.1. MPXV_APOBEC.ipynb - number of reads and allele frequences for each sample with SNP by allele
        
	Archive with plots is available https://drive.google.com/file/d/16i0P8gNPb4JL6YCDeZRfEDB77ym-wwaf/view?usp=sharing
    
    3.2. MPXV_vcf.Rmd - number of reads and allele frequences for C>T mutations for replicas by time and allele
        
	Pictures are located in output_files/pictures/CTmutations_by_replica_time/

4. Simulation 

	Program description: We find potential targets in the genome. For each sample we create dataframe with read ID and potential targets positions. For each allele (TCT, TCG, TCA, TCC) we know number of reads with mutation C>T in each sample. We randomly place mutations by dataframe and count numer of targets that were mutated in the genome for simulation and count number of mutated reads for each potential target position.

	Scripts:
    
	4.1. simulation.sh runs features_samples.py, features_samples_dRNA.py and simulation.py
    
	4.2. features_samples.py and features_samples_dRNA.py extract samples description from 41597_2023_2149_MOESM2_ESM.xlsx
    
	4.3. simulation.py runs simulation for each sample
        
	OUTPUT:
        
	1) simulations_shares.csv - dataframe with real and simulated number of positions that were mutated in the genome and the share of these positions from potential ones
            
	Plots were made using MPXV_APOBEC.ipynb, are located in output_files/pictures/simulation_shares/
        
	2) simulations_Nreads_positions.csv - file with potential targets coordinate and real and simulated number of reads for each position in sample
            
	Plots were made using MPXV_APOBEC.ipynb, are located in output_files/pictures/simulation_Nreads/ 

	File simulations_Nreads_positions.csv was processed in MPXV_APOBEC.ipynb and positions with the greatest and the least number of reads in samples were found.
        
	3) position_with_the_greatest_Nreads.csv
        
	4) positions_with_the_greatest_simulated_Nreads.csv