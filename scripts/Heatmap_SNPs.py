import pandas as pd
import numpy as np
import os
import sys
import itertools
import re
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
import cyvcf2

##min variant frequency
min_VF = 0.01
CALLER = "bcftools"
QUALITY = "MQ10BQ10"
colors = ['#CADCFC', '#8AB6F9', '#00246B']
my_cmap = LinearSegmentedColormap.from_list("custom", colors, N=100)

##vcf file
VCF_DIR = sys.argv[1]
OUTDIR = sys.argv[2]
STAT_FILE = sys.argv[3]


#function to rounding up top position
def round_up(n, decimals=0):
    multiplier = 10**decimals
    return math.ceil(n * multiplier) / multiplier

def get_plot_vcf(directory, outdir, stat_file, min_vf=min_VF, colors=my_cmap, variant_caller=CALLER, quality=QUALITY):
    
    #file with statistics
    #stat_df = pd.read_csv(stat_file, sep="\t")
    #stat_df = stat_df.set_index('Sample')
    
    print(directory)
    samples_list = []
    for filename in os.listdir(directory):
        filepath = os.path.join(directory, filename)
        # checking if it is a file
        if filepath.endswith(".vcf.gz"):
            sample = filename.split('_')[0]
            samples_list.append(sample)
    print(samples_list)
    subst_list = []
    for ref in ['A', 'T', 'C', 'G']:
        for alt in ['A', 'T', 'C', 'G']:
            if ref != alt:
                subst_list.append(ref+">"+alt) 
                
    ##make dataframe from vcf file
    ###dataframe for counting reads number for mutations
    df_Nreads = pd.DataFrame(columns=['substitution']+samples_list)
    df_Nreads['substitution'] = subst_list
    df_Nreads = df_Nreads.set_index('substitution')
    df_Nreads.fillna(0, inplace=True)
    
                
    for filename in os.listdir(directory):
        filepath = os.path.join(directory, filename)
        # checking if it is a file
        if filepath.endswith(".vcf.gz"):
            print("Process: ", filepath)
            sample = filename.split('_')[0]
            
            #n reads with any substitutions
            Nreads_alt_all = 0
            
            # for ERR10513574 there isn't statistics in stat_file
            #if sample == "ERR10513574":
            #    Nreads = stat_df.loc['ERR10963128', 'N_reads_MPXV']
            #else:
            #    Nreads = stat_df.loc[sample, 'N_reads_MPXV']
                
            vcf_file = cyvcf2.VCF(filepath, gts012=True)
            for record in vcf_file:
                REF, ALT = record.REF, record.ALT[0]
                subst=REF+">"+ALT
        
                if variant_caller == "clair3":
                    VAF = record.format('VAF').item()
                    ad = record.format('AD')
                    Nreads_alt = ad[0][1]
                    
                if variant_caller == "bcftools":
                    ad = record.format('AD')
                    Nreads_alt = ad[0][1]
                    VAF = record.format('VAF').item()
                    
                Nreads_alt_all += Nreads_alt
                    
                #fraction of reads with sunstitution (from all number of reads)
                #Nreads_alt_fract = Nreads_alt/Nreads
            
                freq = float(VAF)
                if freq >= min_vf:
                    df_Nreads.at[subst, sample] += Nreads_alt
                        
            #each cell devide to Nreads_alt_all
            df_Nreads[sample] = df_Nreads[sample]/Nreads_alt_all
                        
    df_Nreads.reset_index(inplace=True)
    df_Nreads = df_Nreads.set_index(['substitution'])
    df_Nreads = df_Nreads[samples_list].stack()

    df_Nreads = df_Nreads.to_frame().reset_index()
    df_Nreads = df_Nreads.dropna()
    df_Nreads = df_Nreads.rename(columns={"level_1": "sample", 0:'Nreads'})
    df_Nreads = df_Nreads.reset_index(drop=True)
    df_Nreads = df_Nreads.pivot(index="substitution", columns="sample", values="Nreads")
    print(df_Nreads)

    plt.figure(figsize=(15, 10))
    sns.heatmap(df_Nreads, cmap=colors)
    outfile1 = os.path.join(outdir, variant_caller+quality+'Nreads_all_SNPs.png')
    plt.savefig(outfile1, dpi=800)
    plt.close()

    return
        
get_plot_vcf(VCF_DIR, OUTDIR, STAT_FILE)
