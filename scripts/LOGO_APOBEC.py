import pandas as pd
import numpy as np
import os
import sys
import itertools
import re
from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import seaborn as sns
import cyvcf2

##min variant frequency
min_VF = 0.01
CALLER = "clair3"

##vcf file
VCF_DIR = sys.argv[1]
OUTDIR = sys.argv[2]

##genome file
genome_file = open(sys.argv[3])
genome_dict = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
##MPXV
CHR = 'NC_063383.1'
GENOME = genome_dict[CHR].seq

##palette for pictures
colors_list = ['#CADCFC', '#00246B', '#993955', '#8AB6F9']


def LOGO_plot(directory=VCF_DIR, outdir=OUTDIR, variant_caller=CALLER, min_vf=min_VF, genome=GENOME, my_colors=colors_list):

    ###dataframe for counting reads number for mutations
    position = ["-1", "0", "1"]
    nucl = ['A', 'T', 'C', 'G']
    df_Nreads = pd.DataFrame(columns=['position', 'A', 'T', 'C', 'G'])
    df_Nreads['position'] = position
    df_Nreads = df_Nreads.set_index('position')
    df_Nreads.fillna(0, inplace=True)
        
    for filename in os.listdir(directory):
        filepath = os.path.join(directory, filename)
        # checking if it is a file
        if filepath.endswith(".vcf.gz"):
            print("Process: ", filepath)
            sample = filename.split('_')[0]
            
            vcf_file = cyvcf2.VCF(filepath, gts012=True)
            for record in vcf_file:
                REF, ALT, POS = record.REF, record.ALT[0], record.POS-1
        
                if variant_caller == "clair3":
                    VAF = record.format('VAF').item()
                    ad = record.format('AD')
                    Nreads_alt = ad[0][1]

                if REF == "C":
                    nucl_prev = str(genome[POS-1])
                    nucl_forw = str(genome[POS+1])
                    
                if REF == "G":
                    nucl_forw = str(Seq(genome[POS-1]).reverse_complement())
                    nucl_prev = str(Seq(genome[POS+1]).reverse_complement())
                    REF = "C"
                
                if REF != "C":
                    continue

                freq = float(VAF)
                if freq >= min_vf:
                    df_Nreads.at["-1", nucl_prev] += Nreads_alt
                    df_Nreads.at["0", REF] += Nreads_alt
                    df_Nreads.at["1", nucl_forw] += Nreads_alt

    df_Nreads = df_Nreads.div(df_Nreads.sum(axis=1), axis=0)
    print(df_Nreads)
    
    plt.figure(figsize=(5, 5))
    df_Nreads.plot(kind='bar', stacked=True, color=my_colors, width=1.0, edgecolor='white')
    plt.xticks(rotation=0)
    outfile = os.path.join(outdir, 'LOGO.png')
    plt.savefig(outfile, dpi=800)
    
    return

LOGO_plot()
