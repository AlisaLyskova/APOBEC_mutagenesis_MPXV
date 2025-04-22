import pandas as pd
import numpy as np
import os
import sys
import itertools
import re
from Bio import SeqIO
from Bio.Seq import Seq
import cyvcf2

##min variant frequency
min_VF = 0.01

##directory with vcf files
VCF_DIR = sys.argv[1]

##genome file
genome_file = open(sys.argv[2])
genome_dict = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
##MPXV
CHR = 'NC_063383.1'
GENOME = str(genome_dict[CHR].seq)

STAT_FILE = sys.argv[3]

##Variant Caller: clair3, bcftools
CALLER = "clair3"

##mutations categories
category_list = ['TCT', 'TCA', 'TCC', 'TCG', 'C>T', '*']
category_complement_list = ['AGA', 'TGA', 'GGA', 'CGA', 'C>T', '*']


#function to rounding up top position
def round_up(n, decimals=0):
    multiplier = 10**decimals
    return math.ceil(n * multiplier) / multiplier

def Nread_freq_by_category(filepath, category = category_list, category_complement = category_complement_list, variant_caller=CALLER, min_vf=min_VF, genome=GENOME):
    
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
        
                
    vcf_file = cyvcf2.VCF(filepath, gts012=True)
    for record in vcf_file:
        REF, ALT, POS = record.REF, record.ALT[0], record.POS-1
        if variant_caller == "clair3":
            VAF = record.format('VAF').item()
            ad = record.format('AD')
            Nreads_alt = ad[0][1]
            allele = genome[POS-1:POS+2]
            if VAF >= min_vf:
                df_Nreads.at[allele, ALT] += Nreads_alt
        
    
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
    Nreads = {}
    Nreads_APOBEC_sum = 0
    for i in range(0, len(category[:-2])):
        motif = category[i]
        motif_complement = category_complement[i]
        Nr1 = df_Nreads[(df_Nreads['allele'] == motif) & (df_Nreads['mutation'] == 'C>T')]['Nreads'].sum()
        Nr2 = df_Nreads[(df_Nreads['allele'] == motif_complement) & (df_Nreads['mutation'] == 'G>A')]['Nreads'].sum()
        Nreads[motif] = Nr1+Nr2
        Nreads_APOBEC_sum += Nreads[motif]
    Nreads["C>T"] = Nreads_APOBEC_sum
    Nreads["*"] = df_Nreads['Nreads'].sum() - Nreads_APOBEC_sum
    
    #new_value = "{0} ({1}%)"
    for k, v in Nreads.items():
        pc = round((v / df_Nreads['Nreads'].sum()) * 100, 2)
        #value = new_value.format(v, pc)
        #Nreads[k] = value
        Nreads[k] = pc
    
    return Nreads



def Npos_by_category(filepath, category = category_list, category_complement = category_complement_list, variant_caller=CALLER, min_vf=min_VF, genome=GENOME):
    
    ##make dataframe from vcf file
    ###list with all allele variants of 3 nucleotides
    A = ['A'.join(i) for i in itertools.product('ATCG', repeat=2)]
    T = ['T'.join(i) for i in itertools.product('ATCG', repeat=2)]
    C = ['C'.join(i) for i in itertools.product('ATCG', repeat=2)]
    G = ['G'.join(i) for i in itertools.product('ATCG', repeat=2)]
    allele_list = A + T + C + G

    ###dataframe for counting reads number for mutations
    df_pos = pd.DataFrame(columns=['ref_nucl', 'allele', 'A', 'T', 'C', 'G'])
    df_pos['ref_nucl'] = ['A']*16 + ['T']*16 + ['C']*16 + ['G']*16
    df_pos['allele'] = allele_list
    df_pos = df_pos.set_index('allele')
    df_pos.fillna(0, inplace=True)
    ###replace cells with no mutation (like A>A) with NaN
    for nucl in ['A', 'T', 'C', 'G']:
        df_pos.loc[df_pos['ref_nucl'] == nucl, nucl] = np.nan
        
                
    vcf_file = cyvcf2.VCF(filepath, gts012=True)
    for record in vcf_file:
        REF, ALT, POS = record.REF, record.ALT[0], record.POS-1
        if variant_caller == "clair3":
            VAF = record.format('VAF').item()
            allele = genome[POS-1:POS+2]
            if VAF >= min_vf:
                df_pos.at[allele, ALT] += 1
        
    
    df_pos.reset_index(inplace=True)
    df_pos = df_pos.set_index(['allele'])
    df_pos = df_pos[['A', 'T', 'C', 'G']].stack()
    df_pos = df_pos.to_frame().reset_index()
    df_pos = df_pos.dropna()
    df_pos = df_pos.rename(columns={"level_1": "ALT", 0:'Npos'})
    df_pos['REF'] = df_pos['allele'].apply(lambda x: x[1])
    df_pos['mutation'] = df_pos['REF'] + '>' + df_pos['ALT']
    df_pos = df_pos[['allele', 'REF', 'ALT', 'mutation', 'Npos']]
    df_pos = df_pos.sort_values(by=['mutation', 'allele'])
    df_pos = df_pos.reset_index(drop=True)
    
    #count for category
    pos = {}
    pos_APOBEC_sum = 0
    for i in range(0, len(category[:-2])):
        motif = category[i]
        motif_complement = category_complement[i]
        Nr1 = df_pos[(df_pos['allele'] == motif) & (df_pos['mutation'] == 'C>T')]['Npos'].sum()
        Nr2 = df_pos[(df_pos['allele'] == motif_complement) & (df_pos['mutation'] == 'G>A')]['Npos'].sum()
        pos[motif] = int(Nr1+Nr2)
        pos_APOBEC_sum += pos[motif]
    pos["C>T"] = int(pos_APOBEC_sum)
    pos["*"] = int(df_pos['Npos'].sum() - pos_APOBEC_sum)
    
    return pos




samples = []
for filename in os.listdir(VCF_DIR):
    f = os.path.join(VCF_DIR, filename)
    # checking if it is a vcf file
    if f.endswith(".vcf.gz"):
        sample_id = filename.split('_')[0]
        samples.append(sample_id)

df = pd.DataFrame(index=samples, columns=category_list)

for filename in os.listdir(VCF_DIR):
    f = os.path.join(VCF_DIR, filename)
    # checking if it is a vcf file
    if f.endswith(".vcf.gz"):
        sample_id = filename.split('_')[0]
        print("Process file: ", filename)

        #for number of positions
        new_dict = Npos_by_category(f)

        #for number of reads
        #new_dict = Nread_freq_by_category(f)
        for v, k in new_dict.items():
            df.loc[sample_id][v] = k

print(df)
df = df.reset_index()
df.to_csv(STAT_FILE, sep=",", index=False)

