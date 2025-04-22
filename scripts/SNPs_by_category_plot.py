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

##directory with vcf files
VCF = sys.argv[1]
OUTDIR = sys.argv[2]

##genome file
genome_file = open(sys.argv[3])
genome_dict = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
##MPXV
CHR = 'NC_063383.1'
GENOME = str(genome_dict[CHR].seq)

##Variant Caller: clair3, bcftools
CALLER = "clair3"

##mutations categories
category_list = ['TCT', 'TCA', 'TCC', 'TCG', '*']
category_complement_list = ['AGA', 'TGA', 'GGA', 'CGA', '*']


#function to rounding up top position
def round_up(n, decimals=0):
    multiplier = 10**decimals
    return math.ceil(n * multiplier) / multiplier

def Nread_freq_by_category(input_file, category = category_list, category_complement = category_complement_list, variant_caller=CALLER, min_vf=min_VF, genome=GENOME):
    
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
    
    #count for category
    Nreads = []
    for i in range(0, len(category[:-1])):
        motif = category[i]
        motif_complement = category_complement[i]
        Nr1 = df_Nreads[(df_Nreads['allele'] == motif) & (df_Nreads['mutation'] == 'C>T')]['Nreads'].sum()
        Nr2 = df_Nreads[(df_Nreads['allele'] == motif_complement) & (df_Nreads['mutation'] == 'G>A')]['Nreads'].sum()
        Nreads.append(Nr1+Nr2)
    Nreads.append(df_Nreads['Nreads'].sum() - sum(Nreads))
    
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
    
    #count for category
    freq = []
    for i in range(0, len(category[:-1])):
        motif = category[i]
        motif_complement = category_complement[i]
        freq1 = df_freq[(df_freq['allele'] == motif) & (df_freq['mutation'] == 'C>T')]['freq'].sum()
        freq2 = df_freq[(df_freq['allele'] == motif_complement) & (df_freq['mutation'] == 'G>A')]['freq'].sum()
        freq.append(freq1+freq2)
    freq.append(df_freq['freq'].sum() - sum(freq))
    
    new_list=[]
    for i in range(0, len(category)):
        new_list.append([category[i]]+[Nreads[i]]+[freq[i]])
    
    return new_list


sample_id = os.path.basename(VCF).split('_')[0]
samples = []
samples.append(sample_id)

df = pd.DataFrame(index=samples)
df['allele'] = ''

for sample in samples:
    df.loc[sample]['allele'] = Nread_freq_by_category(VCF)

#transform table
transformed_df = df.explode('allele')
df2 = transformed_df.allele.apply(pd.Series)
df2.columns = ['allele', 'Nreads', 'freq']
result = pd.concat([transformed_df[transformed_df.columns[:-1]], df2], axis=1)
result = result.reset_index()
result = result.rename(columns={"index": "sample", "Nreads":"Number of reads", "freq":"Allele frequency"})


#picture
colors_list = ['#7F3C8D', '#3969AC', '#E68310', '#008695', '#CF1C90']

##Number of reads
plot = sns.catplot(
    data=result, kind="bar",
    x="allele", y="Number of reads", col="sample", hue="allele",
    height=4, aspect=.9, palette=colors_list, dodge=False
)
sns.set_style("whitegrid")
plot.savefig(os.path.join(OUTDIR, "Nreads_by_category.png"), dpi=800)

##Allele frequency
plot = sns.catplot(
    data=result, kind="bar",
    x="allele", y="Allele frequency", col="sample", hue="allele",
    height=4, aspect=.9, palette=colors_list, dodge=False
)
sns.set_style("whitegrid")
plot.savefig(os.path.join(OUTDIR, "freq_by_category.png"), dpi=800)
