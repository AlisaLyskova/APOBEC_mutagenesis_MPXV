import pandas as pd
import numpy as np
import os
import sys
import itertools
import re
import math
from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import seaborn as sns
import cyvcf2

##min variant frequency
min_VF = 0.01
CALLER = "clair3"

##vcf file
VCF = sys.argv[1]
OUTDIR = sys.argv[2]

##genome file
genome_file = open(sys.argv[3])
genome_dict = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
##MPXV
CHR = 'NC_063383.1'
GENOME = str(genome_dict[CHR].seq)

##palette for pictures
colors_list = ['#7F3C8D', '#11A579', '#3969AC', '#E68310', '#008695', '#CF1C90']


#function to rounding up top position
def round_up(n, decimals=0):
    multiplier = 10**decimals
    return math.ceil(n * multiplier) / multiplier

def get_plot_vcf(input_file, ID, variant_caller=CALLER, min_vf=min_VF, genome=GENOME, colors=colors_list):
    
    ##make dataframe from vcf file
    ###list with all allele variants of 3 nucleotides
    A = ['A'.join(i) for i in itertools.product('ATCG', repeat=2)]
    T = ['T'.join(i) for i in itertools.product('ATCG', repeat=2)]
    C = ['C'.join(i) for i in itertools.product('ATCG', repeat=2)]
    G = ['G'.join(i) for i in itertools.product('ATCG', repeat=2)]
    allele_list = A + T + C + G

    ###dataframe for counting reads number for mutations
    df_Nreads = pd.DataFrame(columns=['ref_nucl', 'allele', 'A', 'T', 'C', 'G'])
    df_Nreads['ref_nucl'] = ['A']*16 + ['T']*16 + ['C']*16 + ['G']*16
    df_Nreads['allele'] = allele_list
    df_Nreads = df_Nreads.set_index('allele')
    df_Nreads.fillna(0, inplace=True)
    ###replace cells with no mutation (like A>A) with NaN
    for nucl in ['A', 'T', 'C', 'G']:
        df_Nreads.loc[df_Nreads['ref_nucl'] == nucl, nucl] = np.nan
    
    ###dataframe for counting allele's frequencies for mutations
    df_freq = pd.DataFrame(columns=['ref_nucl', 'allele', 'A', 'T', 'C', 'G'])
    df_freq['ref_nucl'] = ['A']*16 + ['T']*16 + ['C']*16 + ['G']*16
    df_freq['allele'] = allele_list
    df_freq = df_freq.set_index('allele')
    df_freq.fillna(0, inplace=True)
    ###replace cells with no mutation (like A>A) with NaN
    for nucl in ['A', 'T', 'C', 'G']:
        df_freq.loc[df_freq['ref_nucl'] == nucl, nucl] = np.nan
        

    vcf_file = cyvcf2.VCF(input_file, gts012=True)
    for record in vcf_file:
        REF, ALT, POS = record.REF, record.ALT[0], record.POS-1
        allele = genome[POS-1:POS+2]
        
        if variant_caller == "clair3":
            ad = record.format('AD')
            Nreads_alt = ad[0][1]
            freq = record.format('VAF').item()
            df_Nreads.at[allele, ALT] += Nreads_alt
            df_freq.at[allele, ALT] += freq        
    
    df_Nreads.reset_index(inplace=True)
    df_Nreads = df_Nreads.set_index(['allele'])
    df_Nreads = df_Nreads[['A', 'T', 'C', 'G']].stack()
    df_Nreads = df_Nreads.to_frame().reset_index()
    df_Nreads = df_Nreads.dropna()
    df_Nreads = df_Nreads.rename(columns={"level_1": "ALT", 0:'Nreads'})
    df_Nreads['REF'] = df_Nreads['allele'].apply(lambda x: x[1])
    df_Nreads['mutation'] = df_Nreads['REF'] + '>' + df_Nreads['ALT']
    df_Nreads = df_Nreads[['allele', 'REF', 'ALT', 'mutation', 'Nreads']]
    df_Nreads = df_Nreads.sort_values(by=['mutation', 'allele'])
    df_Nreads = df_Nreads.reset_index(drop=True)
    
    df_freq.reset_index(inplace=True)
    df_freq = df_freq.set_index(['allele'])
    df_freq = df_freq[['A', 'T', 'C', 'G']].stack()
    df_freq = df_freq.to_frame().reset_index()
    df_freq = df_freq.dropna()
    df_freq = df_freq.rename(columns={"level_1": "ALT", 0:'freq'})
    df_freq['REF'] = df_freq['allele'].apply(lambda x: x[1])
    df_freq['mutation'] = df_freq['REF'] + '>' + df_freq['ALT']
    df_freq = df_freq[['allele', 'REF', 'ALT', 'mutation', 'freq']]
    df_freq = df_freq.sort_values(by=['mutation', 'allele'])
    df_freq = df_freq.reset_index(drop=True)

    # picture
    custom_palette = []
    for x in colors:
        custom_palette.extend([x]*16)
    ##divide dataframe into two parts
    df1_Nreads = df_Nreads[df_Nreads['REF'].isin(['A', 'C'])]
    df2_Nreads = df_Nreads[df_Nreads['REF'].isin(['G', 'T'])]
    df1_freq = df_freq[df_freq['REF'].isin(['A', 'C'])]
    df2_freq = df_freq[df_freq['REF'].isin(['G', 'T'])]

    mutation_list1 = df1_Nreads.mutation.to_list()
    mutation_list1 = list(set(mutation_list1))
    mutation_list1.sort()
    mutation_list2 = df2_Nreads.mutation.to_list()
    mutation_list2 = list(set(mutation_list2))
    mutation_list2.sort()

    ###picture for number of reads
    fig1 = plt.figure(figsize=(16, 7))
    ax11 = plt.subplot2grid(shape=(2, 1), loc=(0, 0))
    ax12 = plt.subplot2grid(shape=(2, 1), loc=(1, 0))

    ax11 = sns.barplot(x=df1_Nreads.index, y=df1_Nreads.Nreads, palette=custom_palette, ax=ax11)
    ax12 = sns.barplot(x=df2_Nreads.index, y=df2_Nreads.Nreads, palette=custom_palette, ax=ax12)
    ax11.set_xticklabels(df1_Nreads.allele, rotation=90, horizontalalignment='center')
    ax12.set_xticklabels(df2_Nreads.allele, rotation=90, horizontalalignment='center')

    ####change lim of picture
    top_position = int(df_Nreads.Nreads.max())
    round_pos = 0
    while top_position >= 10:
        top_position = int(top_position / 10)
        round_pos += 1
    new_top_position = round_up(df_Nreads.Nreads.max(), -round_pos)

    ####rectangle must place 10% of space
    part_ten = int(new_top_position/10)
    ax11.set_ylim(0, new_top_position+part_ten)
    ax12.set_ylim(0, new_top_position+part_ten)

    ####add rectangles
    for i in range(0, len(colors)):
        ax11.add_patch(Rectangle((i*16, new_top_position), 16, part_ten, facecolor = colors[i]))
        ax11.text(i*16+8, new_top_position+int(part_ten/2), mutation_list1[i], horizontalalignment='center', verticalalignment='center', size='x-large', color='white', weight='semibold')
        ax12.add_patch(Rectangle((i*16, new_top_position), 16, part_ten, facecolor = colors[i]))
        ax12.text(i*16+8, new_top_position+int(part_ten/2), mutation_list2[i], horizontalalignment='center', verticalalignment='center', size='x-large', color='white', weight='semibold')

    ####save picture
    #fig1.suptitle('Number of reads by mutation type for sample '+ID, fontsize=18)
    fig1.suptitle('Количество прочтений по типу мутации для образца '+ID, fontsize=18)
    plt.tight_layout()
    plt.close(fig1)


    ###picture for freq
    fig2 = plt.figure(figsize=(16, 7))
    ax21 = plt.subplot2grid(shape=(2, 1), loc=(0, 0))
    ax22 = plt.subplot2grid(shape=(2, 1), loc=(1, 0))

    ax21 = sns.barplot(x=df1_freq.index, y=df1_freq.freq, palette=custom_palette, ax=ax21)
    ax22 = sns.barplot(x=df2_freq.index, y=df2_freq.freq, palette=custom_palette, ax=ax22)
    ax21.set_xticklabels(df1_freq.allele, rotation=90, horizontalalignment='center')
    ax22.set_xticklabels(df2_freq.allele, rotation=90, horizontalalignment='center')

    ####change lim of picture
    top_position = df_freq.freq.max()
    if top_position >= 1:
        new_top_position = top_position + 0.1
    else:
        round_pos = 0
        while top_position < 1:
            top_position = top_position * 10
            round_pos += 1
        new_top_position = round_up(df_freq.freq.max(), round_pos)
        

    ####rectangle must place 10% of space
    part_ten = new_top_position/10
    ax21.set_ylim(0, new_top_position+part_ten)
    ax22.set_ylim(0, new_top_position+part_ten)
    

    ####add rectangles
    for i in range(0, len(colors)):
        ax21.add_patch(Rectangle((i*16, new_top_position), 16, part_ten, facecolor = colors[i]))
        ax21.text(i*16+8, new_top_position+part_ten/2, mutation_list1[i], horizontalalignment='center', verticalalignment='center', size='x-large', color='white', weight='semibold')
        ax22.add_patch(Rectangle((i*16, new_top_position), 16, part_ten, facecolor = colors[i]))
        ax22.text(i*16+8, new_top_position+part_ten/2, mutation_list2[i], horizontalalignment='center', verticalalignment='center', size='x-large', color='white', weight='semibold')

    #fig2.suptitle('Allele frequency by mutation type for sample '+ID, fontsize=18)
    fig2.suptitle('Аллельная частота по типу мутации для образца '+ID, fontsize=18)
    plt.tight_layout()
    plt.close(fig2)
    
    return [fig1, fig2]


fig_list = ['Nreads', 'freq']
sample_id = os.path.basename(VCF).split('_')[0]
plot_list = get_plot_vcf(VCF, sample_id)
for i in range(0, len(plot_list)):
            plot_list[i].savefig(os.path.join(OUTDIR, sample_id+'_'+fig_list[i]+'.png'))
