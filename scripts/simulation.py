import pandas as pd
import numpy as np
import pysam
import re
import os
import sys
import cyvcf2
from Bio import SeqIO
from Bio.Seq import Seq
import random
import matplotlib.pyplot as plt
import seaborn as sns

sys.stdout.flush()
# INPUT
## min for variant frequency
min_VF = 0.01
plot_color1 = '#00246B'
plot_color2 = '#993955'

##directory with bam files
SAMPLE_ID = sys.argv[1]
BAM_FILE = sys.argv[2]
VCF_FILE = sys.argv[3]
TARGETS_FILE = sys.argv[4]

## genome file
genome_file = open(sys.argv[5])
genome_dict = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
##MPXV
CHR = 'NC_063383.1'
GENOME = str(genome_dict[CHR].seq)

MOTIF = sys.argv[6]

OUTDIR = sys.argv[7]

N_READS_MAX = int(sys.argv[8])

# FUNCTIONS
## find coordinates of potential targets with motif in sample
## output: dataframe with potential targets coordinates in the genome for sample
def find_potential_targets(bam_path, motif, genome=GENOME):
    
    #CG_coordinates = [i.start() for i in re.finditer('C', genome, flags=0)] + [i.start() for i in re.finditer('G', genome, flags=0)]
    #CG_coordinates.sort()
    #targets_coordinates = [i-1 for i in CG_coordinates if (genome[i-1:i+2] == motif or genome[i-1:i+2] == str(Seq(motif).reverse_complement())) ]
    
    targets_coordinates = [i.start() for i in re.finditer(motif, genome, flags=0)] + [i.start() for i in re.finditer(str(Seq(motif).reverse_complement()), genome, flags=0)]
    targets_coordinates.sort()
    print(len(targets_coordinates))

    ### find max number of reads
    global nreads_max_sample
    nreads_max_sample = 1

    ### filter targets that are covered by at leats 1 read
    #targets_coordinates_filtered = []
    #bamfile = pysam.AlignmentFile(bam_path, "rb")
    #for i in targets_coordinates:
        #### write i only if reads number > 0
        #### bam file start with 1, i coordinates start with 0
    #    counter = 0
    #    for read in bamfile.fetch("NC_063383.1", i+2, i+3):
    #        counter += 1
    #    if counter >= 1:
    #        targets_coordinates_filtered.append(i)
    #    if counter > nreads_max_sample:
    #        nreads_max_sample = counter

    #out_df = pd.DataFrame({'target_coordinate':targets_coordinates_filtered})

    out_df = pd.DataFrame({'target_coordinate':targets_coordinates})

    print(out_df.head())
    return out_df


## for all targets with motif find number of reads with mutation in vcf file 
## output: dataframe from input with potantial targets and number of reads with mutation in real sample
def Nreads_mutated_sample(input_file, targets_file, min_vf=min_VF, genome=GENOME, virus_chr=CHR, motif=MOTIF):

    #filter 3d position in motif
    nucl3 = motif[2]

    vcf_file = cyvcf2.VCF(input_file, gts012=True)
    
    for record in vcf_file:
        VAF = record.format('VAF').item()
        ad = record.format('AD')
        Nreads_alt = ad[0][1]
    
    targets_df = pd.read_csv(targets_file, names=["target_coordinate"])
    targets_df = pd.concat([targets_df, pd.DataFrame({'Nreads_mutated_real':[0]*len(targets_df)})], axis=1)
    for row in range(0, len(targets_df)):
        pos_target = targets_df.loc[row, 'target_coordinate']
        region = "{0}:{1}-{2}"
        new_region = region.format(virus_chr, pos_target-1, pos_target+5)
        for v in vcf_file(new_region):
            REF, ALT, POS = v.REF, v.ALT[0], v.POS
            if POS == pos_target + 2:
                if REF == "C":
                    prev_nucl = genome[POS-2]
                    next_nucl = genome[POS]
                    if prev_nucl != "T":
                        continue
                    if next_nucl != nucl3:
                        continue
                if REF == "G":
                    prev_nucl = genome[POS]
                    next_nucl = genome[POS-2]
                    if prev_nucl != "A":
                        continue
                    if next_nucl != str(Seq(nucl3).reverse_complement()):
                        continue

                Nreads_alt = v.format('AD')[0][1]
                targets_df.at[row, 'Nreads_mutated_real'] = Nreads_alt
                              
    return targets_df


# find targets for all reads from bam file to simulation
def simulation(bam_path, targets_df, outfile_shares, outfile_plot, nreads_max_sample=N_READS_MAX):
    
    ### for each read in bam file find potential positions for mutation
    read_positions_dict = dict()
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    
    print("Extract reads ids")
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
    
    print("Generation")
    ### generate list with 0 - there isn't mutation in this position and 1 - there is mutation in this position
    Nmutated_reads = int(targets_df.Nreads_mutated_real.sum()) #сумма всех ридов с мутацией
    try:
        new_list = [1]*Nmutated_reads+[0]*(len(df)-Nmutated_reads)
    except TypeError:
        print(Nmutated_reads)
        return

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
    
    print("Calculation shares")
    ### write proportion of mutated positions from all potential positions
    with open(outfile_shares, 'w') as fout:
        N_mutated_pos = len(targets_df[targets_df['Nreads_mutated_real'] != 0])
        share_real = N_mutated_pos / len(targets_df)
        fout.write("Number of real positions: " + str(N_mutated_pos) + "\n")
        fout.write("proportion: " + str(share_real) + "\n")
        fout.write("Mean number of mutated positions in simulations: " + str(sum(N_mutated_pos_list)/len(N_mutated_pos_list)) + "\n")
        fout.write("proportion: " + str(sum(shares)/len(shares)) + "\n")


    print("Create plot")
    ### plot with N reads
    targets_df['generation_mean'] = targets_df.drop('Nreads_mutated_real', axis=1).drop('target_coordinate', axis=1).mean(axis=1)
    targets_df = targets_df.rename(columns={'Nreads_mutated_real': 'sample', 'target_coordinate': 'coordinate'})
    #targets_df.to_csv("/data1/lyskovaa/mpxv_pipeline/simulation/ERR10513574/targets_df_simulation.csv", sep=",", index=False)
    
    #normalize number of reads
    xmin = 1
    xmax = nreads_max_sample
    normalized_df = (targets_df.drop('coordinate', axis=1) - xmin)/xmax
    normalized_df['coordinate'] = targets_df['coordinate']

    ylim = normalized_df['sample'].max(axis=0)
    print(ylim)

    dfm = normalized_df[['coordinate', 'sample', 'generation_mean']].melt('coordinate', var_name='generation', value_name='vals')

    df1 = dfm[dfm['coordinate'] < 50000]
    df2 = dfm[(dfm['coordinate'] >= 50000) & (dfm['coordinate'] <= 100000)]
    df3 = dfm[(dfm['coordinate'] > 100000) & (dfm['coordinate'] < 150000)]
    df4 = dfm[dfm['coordinate'] >= 150000]

    sns.set_style("whitegrid")
    fig = plt.figure(figsize=(16, 8))
    ax1 = plt.subplot2grid(shape=(4, 1), loc=(0, 0))
    ax2 = plt.subplot2grid(shape=(4, 1), loc=(1, 0))
    ax3 = plt.subplot2grid(shape=(4, 1), loc=(2, 0))
    ax4 = plt.subplot2grid(shape=(4, 1), loc=(3, 0))

    ax1 = sns.lineplot(data=df1, x="coordinate", y="vals", hue="generation", palette=[plot_color1, plot_color2], ax=ax1)
    ax2 = sns.lineplot(data=df2, x="coordinate", y="vals", hue="generation", palette=[plot_color1, plot_color2], legend=False, ax=ax2)
    ax3 = sns.lineplot(data=df3, x="coordinate", y="vals", hue="generation", palette=[plot_color1, plot_color2], legend=False, ax=ax3)
    ax4 = sns.lineplot(data=df4, x="coordinate", y="vals", hue="generation", palette=[plot_color1, plot_color2], legend=False, ax=ax4)
    ax1.set_ylim(0, ylim+0.05)
    ax2.set_ylim(0, ylim+0.05)
    ax3.set_ylim(0, ylim+0.05)
    ax4.set_ylim(0, ylim+0.05)

    fig.savefig(outfile_plot)
    #targets_df.to_csv(outfile_nReads, sep='\t', index=False)



# BODY
## find potential targets in genome
#print("Searching for potential targets for " + MOTIF)
#RES_DF = find_potential_targets(BAM_FILE, MOTIF)
#print("Searching for potential targets for " + MOTIF + ": done")

## find number of reads with mutation in vcf file for all targets
print("Counting number of reads with mutation in sample for " + MOTIF)
RES_DF = Nreads_mutated_sample(VCF_FILE, TARGETS_FILE)
print("Counting number of reads with mutation in sample for " + MOTIF + ": done")

## result tables
OUT_SHARES = os.path.join(OUTDIR, SAMPLE_ID+"_simulations_shares_" + MOTIF + ".csv")
OUT_PLOT = os.path.join(OUTDIR, SAMPLE_ID+"_simulations_Nreads_" + MOTIF + ".png")

print("Simulation for " + MOTIF)
simulation(BAM_FILE, RES_DF, OUT_SHARES, OUT_PLOT)
print("Simulation for " + MOTIF + ": done")


