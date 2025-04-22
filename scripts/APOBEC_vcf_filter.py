import pandas as pd
import numpy as np
import sys
import cyvcf2
import os
from Bio import SeqIO
from Bio.Seq import Seq

##min variant frequency
min_VF = 0.01

##vcf file
VCF = sys.argv[1]
OUT = sys.argv[2]

##genome file
genome_file = open(sys.argv[3])
genome_dict = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
##MPXV
CHR = 'NC_063383.1'
GENOME = str(genome_dict[CHR].seq)

def filter_APOBEC(input_file, output_file, genome=GENOME, min_vaf=min_VF):
    file_name = os.path.basename(input_file)
    caller = file_name.split("_")[1]
    feature = file_name.split("_")[2]
    vcf_type = caller + "_" + feature

    vcf_file = cyvcf2.VCF(input_file, gts012=True)
    
    vcf_out = output_file
    w = cyvcf2.Writer(vcf_out, vcf_file)
    
    for record in vcf_file:
        REF, ALT, POS = record.REF, record.ALT[0], record.POS
        
        if (REF == "C" and ALT == "T") or (REF == "G" and ALT == "A"):
        
            if caller == "bcftools":
                vaf = record.format('VAF').item()
                
            if caller == "lofreq":
                vaf = record.INFO.get('AF')
                
            if caller == "clair3":
                vaf = record.format('VAF').item()

            #motif filter
            if REF == "C":
                prev_nucl = genome[POS-2]
                if prev_nucl != "T":
                    continue
            if REF == "G":
                prev_nucl = genome[POS]
                if prev_nucl != "A":
                    continue

            if vaf >= min_vaf:
                w.write_record(record)

                #print(REF)
                #if REF == "C":
                #    print(genome[POS-2:POS])
                #if REF == "G":
                #    print(genome[POS-1:POS+1])
                
    w.close(); vcf_file.close()

filter_APOBEC(VCF, OUT)

