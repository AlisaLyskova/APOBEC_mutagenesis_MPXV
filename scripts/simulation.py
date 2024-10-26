import pandas as pd
import numpy as np
import pysam
import re
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import random

# INPUT
## min for variant frequency
min_VF = 0.001

## genome file
genome_file = open("data/NC_063383.1.fa")
genome_dict = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
GENOME=str(genome_dict['NC_063383.1'].seq)

## paths to vcf and bam files
vcf_dir = ""
bam_dir = ""


# FUNCTIONS
## find coordinates of potential targets with motif in sample
## output: dataframe with potential targets coordinates in the genome for sample
def find_potential_targets(bam_path, motif, genome=GENOME):
    
    CG_coordinates = [i.start() for i in re.finditer('C', genome, flags=0)] + [i.start() for i in re.finditer('G', genome, flags=0)]
    CG_coordinates.sort()
    targets_coordinates = [i-1 for i in CG_coordinates if (genome[i-1:i+2] == motif or genome[i-1:i+2] == str(Seq(motif).reverse_complement())) ]

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
    out_df = pd.DataFrame({'target_coordinate':targets_coordinates_filtered})
    return out_df


## for all targets with motif find number of reads with mutation in vcf file 
## output: dataframe from input with potantial targets and number of reads with mutation in real sample
def Nreads_mutated_sample(vcf_file, targets_df, min_vf=min_VF, genome=GENOME):

    ###read vcf file
    ###skip description in vcf file
    with open(vcf_file, 'r') as f:
        reader=f.readlines()
        row = 0
        while reader[row].startswith('##') == True:
            row += 1
    df_vcf = pd.read_csv(vcf_file, skiprows=row, sep='\t')
    df_vcf = df_vcf.dropna()
    
    targets_df = pd.concat([targets_df, pd.DataFrame({'Nreads_mutated_real':[0]*len(targets_df)})], axis=1)
    for row in range(0, len(targets_df)):
        pos_target = targets_df.loc[row, 'target_coordinate']
        pos_snp = pos_target + 2
        df_vcf_pos = df_vcf[df_vcf['POS'] == pos_snp]
        if len(df_vcf_pos) != 0:
            df_vcf_pos = df_vcf_pos.reset_index(drop=True)
            allele = genome[pos_target:pos_target+3]
            if len(df_vcf_pos) != 1:
                print(len(df_vcf_pos))
                print(pos_target)
            alt_nucl = df_vcf_pos.loc[0, 'ALT']
            if (len(alt_nucl) != 1) or (len(df_vcf_pos.loc[0, 'REF']) != 1):
                continue
            if (allele[1] == 'C' and alt_nucl == 'T') or (allele[1] == 'G' and alt_nucl == 'A'):
                FORMAT = df_vcf_pos.iloc[0, -2].split(':')
                AD_ind = FORMAT.index('AD')
                AF_ind = FORMAT.index('AF')
                AD = df_vcf_pos.iloc[0, -1].split(':')[AD_ind]
                Nreads_alt = int(AD.split(',')[1])
                AF = df_vcf_pos.iloc[0, -1].split(':')[AF_ind]
                freq = float(AF)
                if freq >= min_vf:
                    targets_df.at[row, 'Nreads_mutated_real'] = Nreads_alt
                    
    return targets_df


# find targets for all reads from bam file to simulation
def simulation(bam_path, targets_df, outfile_shares, outfile_nReads):
    
    ### for each read in bam file find potential positions for mutation
    read_positions_dict = dict()
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    
    for read in bamfile.fetch():
        ### list with positions that were aligned to reference
        aligned_pos_list = read.get_reference_positions(full_length=False)
        read_target_coordinates = list(set(aligned_pos_list) & set(targets_df.target_coordinate.tolist()))
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
    Nmutated_reads = targets_df.Nreads_mutated_real.sum() #сумма всех ридов с мутацией
    new_list = [1]*Nmutated_reads+[0]*(len(df)-Nmutated_reads)

    ### add shuffled list to df with reads and target positions for each simulation
    #### count shares
    shares = []
    N_mutated_pos_list = []
    new_df = pd.DataFrame(0, index=np.arange(len(targets_df)), columns=["generation_"+str(i) for i in range(1, 101)])
    targets_df = pd.concat([targets_df, new_df], axis=1)
    
    for i in range(1, 101):
        shuffled_list = sorted(new_list, key=lambda x: random.random())
        df['generation'] = shuffled_list

        #### count number of reads with mutation for each position in targets_genome
        for row in range(0, len(targets_df)):
            x = targets_df.loc[row, 'target_coordinate']
            df2 = df[df['motif_coordinates'] == x]
            Nreads_x = df2['generation'].sum()
            targets_df.at[row, 'generation_' + str(i)] = Nreads_x
            
        N_mutated_pos = len(targets_df[targets_df['generation_' + str(i)] != 0])
        share = N_mutated_pos / len(targets_df)
        shares.append(share)
        N_mutated_pos_list.append(N_mutated_pos)
    
    ### write proportion of mutated positions from all potential positions
    with open(outfile_shares, 'w') as fout:
        N_mutated_pos = len(targets_df[targets_df['Nreads_mutated_real'] != 0])
        share_real = N_mutated_pos / len(targets_df)
        fout.write("Number of real positions: " + str(N_mutated_pos) + "\n")
        fout.write("proportion: " + str(share_real) + "\n")
        fout.write("Mean number of mutated positions in simulations: " + str(sum(N_mutated_pos_list)/len(N_mutated_pos_list)) + "\n")
        fout.write("proportion: " + str(sum(shares)/len(shares)) + "\n")
    
    targets_df.to_csv(outfile_nReads, sep='\t', index=False)



# BODY
## input: sample id, allele
SAMPLE_ID = sys.argv[1]
ALLELE = sys.argv[2]

## find potential targets in genome
BAM_FILE = os.path.join(bam_dir, SAMPLE_ID+'.bam')
print("Searching for potential targets for " + ALLELE)
RES_DF = find_potential_targets(BAM_FILE, ALLELE)
print("Searching for potential targets for " + ALLELE + ": done")

## find number of reads with mutation in vcf file for all targets
VCF_FILE = os.path.join(vcf_dir, SAMPLE_ID+'.vcf')
print("Counting number of reads with mutation in sample for " + ALLELE)
RES_DF = Nreads_mutated_sample(VCF_FILE, RES_DF)
print("Counting number of reads with mutation in sample for " + ALLELE + ": done")

# OUTPUT
OUT_SHARES = "data/simulation/simulations_shares_" + ALLELE + ".txt"
OUT_NREADS = "data/simulation/simulations_Nreads_" + ALLELE + ".csv"

print("Simulation for " + ALLELE)
simulation(BAM_FILE, RES_DF, OUT_SHARES, OUT_NREADS)
print("Simulation for " + ALLELE + ": done")
