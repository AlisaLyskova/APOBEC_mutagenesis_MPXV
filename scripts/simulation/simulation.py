import pandas as pd
import pysam
import re
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from random import shuffle
        
# OUTPUT
OUT_SHARES = "simulations_shares.csv"
OUT_NREADS = "simulations_Nreads_positions.csv"
OUT_LOG = "simulations.log"

# INPUT
## min and max for variant frequency
min_VF = 0.001
max_VF = 0.9

## genome file
genome_file = open("NC_063383.1.fa")
genome_dict = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))

## paths to vcf and bam files
vcf_dir = ""
bam_dir = ""

# FUNCTIONS

## find coordinates of potential targets with motif in sample
## output: list with potential targets in the genome for sample
def find_potential_targets(bam_path, motif, genome=str(genome_dict['NC_063383.1'].seq)):
    motif_complement = str(Seq(motif).reverse_complement())
    targets_coordinates = [i.start() for i in re.finditer(motif, genome, flags=0)] + [i.start() for i in re.finditer(motif_complement, genome, flags=0)]

    ### filter targets that are covered by at leats 1 read
    targets_coordinates_filtered = []
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    for i in targets_coordinates:
        #### write i only if reads number > 0
        #### bam file start with 1, i coordinates start with 0
        for read in bamfile.fetch("NC_063383.1", i+2, i+3):
            if i+2 in read.get_reference_positions():
                targets_coordinates_filtered.append(i)
                break
    return targets_coordinates_filtered


## extract N reads with mutation for each position for allele from vcf file
## output: list of lists with allele coordinates and N reads for each coordinate
def Nreads_mutated_sample(vcf_file, motif, min_vf=min_VF, max_vf=max_VF, genome=str(genome_dict['NC_063383.1'].seq)):

    ###read vcf file
    ###skip description in vcf file
    with open(vcf_file, 'r') as f:
        reader=f.readlines()
        row = 0
        while reader[row].startswith('##') == True:
            row += 1
    df_vcf = pd.read_csv(vcf_file, skiprows=row, sep='\t')
    df_vcf = df_vcf.dropna()

    ### count N reads for positions with motif and C>T or G>A
    motif_positions = []
    Nreads_pos = []
    for row in range(0, len(df_vcf)):
        mutation_pos = df_vcf.loc[row, 'POS'] - 1
        allele = genome[mutation_pos-1:mutation_pos+2]
        alt_nucl = df_vcf.loc[row, 'ALT']
        if (len(alt_nucl) != 1) or (len(df_vcf.loc[row, 'REF']) != 1) or (len(allele) != 3):
            continue
        if ((allele == motif) and (alt_nucl == 'T')) or ((allele == str(Seq(motif).reverse_complement())) and (alt_nucl == 'A')):
            FORMAT = df_vcf.iloc[row, -2].split(':')
            AD_ind = FORMAT.index('AD')
            AF_ind = FORMAT.index('AF')
            AD = df_vcf.iloc[row, -1].split(':')[AD_ind]
            Nreads_alt = int(AD.split(',')[1])
            AF = df_vcf.iloc[row, -1].split(':')[AF_ind]
            freq = float(AF)
            if freq >= min_vf and freq < max_vf:
                motif_positions.append(mutation_pos-1)
                Nreads_pos.append(Nreads_alt)      
    return [motif_positions, Nreads_pos]
    
   
## for each read find potential target positions then randomly place 0 and 1 (1 means that position is mutated) with the sum of units = Nmutated_reads (real number of reads with mutation and allele) and count positions in the genome that were mutated and number of reads with mutation for each potential target position
def simulation(bam_path, sample_id, time, group, file_type, allele, targets_genome, Nmutated_reads, outfile_shares=OUT_SHARES, outfile_nReads=OUT_NREADS, outfile_log=OUT_LOG):
    
    ### for each read in bam file find potential positions for mutation
    read_positions_dict = dict()
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    
    for read in bamfile.fetch():
        ### list with positions that were aligned to reference
        aligned_pos_list = read.get_reference_positions(full_length=False)
        read_target_coordinates = list(set(aligned_pos_list) & set(targets_genome))
        read_positions_dict[read.query_name] = read_target_coordinates

    ### create dataframe with id read and potential positions for mutation for each read
    reads_list = []
    positions_list = []
    for key, value in read_positions_dict.items():
        reads_list += [key]*len(value)
        positions_list += value
    d = {'id_read': reads_list, 'motif_coordinates': positions_list}
    df = pd.DataFrame(d)

    ### generate list with 0 - there isn't mutation in this position and 1 - there is mutation in this position
    new_list = [1]*Nmutated_reads+[0]*(len(df)-Nmutated_reads)

    ### add shuffled list to df with reads and target positions for each simulation
    with open(outfile_shares, "a") as fout1, open(outfile_nReads, "a") as fout2, open(outfile_log, "a") as fout_log:
    
        for i in range(1, 11):
            shuffle(new_list)
            df['generation_'+str(i)] = new_list
            
            #### count N reads with mutation for each position in targets_genome
            Nreads_simulated_for_targets_genome = []
            #### count genome positions that were mutated at least 1 time
            genome_n_mutations = 0
            for x in targets_genome:
                df2 = df[df['motif_coordinates'] == x]
                Nreads_x = df2['generation_'+str(i)].sum()
                Nreads_simulated_for_targets_genome.append(Nreads_x)
                if Nreads_x != 0:
                    genome_n_mutations += 1
            
            fout1.write(sample_id+"\t"+time+"\t"+group+"\t"+file_type+"\t"+allele+"\t"+'generation_'+str(i)+"\t"+str(genome_n_mutations) +"\t"+str(genome_n_mutations/len(targets_genome))+"\n")
            fout2.write(sample_id+"\t"+allele+"\t"+"Nreads_generation_"+str(i)+"\t"+'\t'.join(list(map(str, Nreads_simulated_for_targets_genome)))+'\n')
            fout_log.write(sample_id+"\t"+allele+"\t"+"Nreads_generation_"+str(i)+"\t"+str(sum(Nreads_simulated_for_targets_genome))+"\n")
    
    
# BODY
## samples features
features = sys.argv[1].split('_')
SAMPLE_ID, TIME, GROUP, FILE_TYPE, ALLELE = features[0], features[1], features[2], features[3], features[4]

## find potential targets in genome
BAM_FILE = os.path.join(bam_dir, SAMPLE_ID+'.bam')
potential_targets_genome = find_potential_targets(BAM_FILE, ALLELE)


## count reads with mutations for each position in real sample
VCF_FILE = os.path.join(vcf_dir, SAMPLE_ID+'.vcf')
### list of lists with allele coordinates and N reads for each coordinate in sample
mutated_pos_Nreads_sample = Nreads_mutated_sample(VCF_FILE, ALLELE)

### write proportion mutated positions from all potential positions
share_real_mutations = len(mutated_pos_Nreads_sample[0]) / len(potential_targets_genome)
with open(OUT_SHARES, "a") as fout1:
    fout1.write(SAMPLE_ID+"\t"+TIME+"\t"+GROUP+"\t"+FILE_TYPE+"\t"+ALLELE+"\t"+"real\t"+str(len(mutated_pos_Nreads_sample[0]))+"\t"+str(share_real_mutations)+"\n")
    
### write to OUT_NREADS file coordinates with potential targets and N reads mutated for each this coordinate
Nreads_mutated_for_potential_targets = []
for x in potential_targets_genome:
    #### if position is mutated
    if x in mutated_pos_Nreads_sample[0]:
        Nreads_x = mutated_pos_Nreads_sample[1][mutated_pos_Nreads_sample[0].index(x)]
        Nreads_mutated_for_potential_targets.append(Nreads_x)
    else:
        Nreads_mutated_for_potential_targets.append(0)
with open(OUT_NREADS, "a") as fout2:
    fout2.write(SAMPLE_ID+"\t"+ALLELE+"\t"+"targets_coordinates\t"+'\t'.join(list(map(str, potential_targets_genome)))+'\n')
    fout2.write(SAMPLE_ID+"\t"+ALLELE+"\t"+"Nreads_mutated_sample\t"+'\t'.join(list(map(str, Nreads_mutated_for_potential_targets)))+'\n')
with open(OUT_LOG, "a") as fout_log:
    fout_log.write(SAMPLE_ID+"\t"+ALLELE+"\t"+"Nreads_mutated_real\t"+str(sum(Nreads_mutated_for_potential_targets))+"\n")
    

## count reads with mutations for each position in simulation
simulation(BAM_FILE, SAMPLE_ID, TIME, GROUP, FILE_TYPE, ALLELE, potential_targets_genome, sum(Nreads_mutated_for_potential_targets))

