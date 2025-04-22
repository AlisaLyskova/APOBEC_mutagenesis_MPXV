import sys
import pandas as pd
import cyvcf2
import os
import re
from Bio import SeqIO
from Bio.Seq import Seq

sys.path.insert(1, '/data1/biosoft/RNAsselem/RNAsselem/rnasselem')
import rnasselem as r

VCF = sys.argv[1]
FASTA = sys.argv[2]
OUTDIR = sys.argv[3]
feature_for_name = sys.argv[4]

def load_genome(fasta_file):
    genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    new_genome = genome["NC_063383.1"].seq
    return str(new_genome)


def vcf_reader(filepath=VCF, genome_file=FASTA, outdir=OUTDIR):
    vcf_file = cyvcf2.VCF(filepath, gts012=True)

    column_names=["POS", "structure_type"]
    df = pd.DataFrame(columns=column_names)

    for record in vcf_file:
        POS = record.POS - 1
        if POS < 100 or (POS > 197209-105):
            continue
        new_dict = r.get_predicted_structure_type_full(outdir, POS, genome_file)
        print(new_dict)
        new_row = [POS, new_dict["structure_type"]]
        df.loc[len(df)] = new_row

    df.to_csv(os.path.join(outdir, "APOBEC_subst_secondary_structure_"+feature_for_name+".txt"), sep='\t', index=False)
    print(df.groupby(by=["structure_type"]).count())

def find_targets(motif="TC", genome_file=FASTA, Chr="NC_063383.1", outdir=OUTDIR):

    genome = load_genome(genome_file)

    motif2 = str(Seq(motif).reverse_complement())
    #write for TC and GA position of substitution (+1 and 0 accordingly)
    targets_list = [i.start()+1 for i in re.finditer('(?={0})'.format(re.escape(motif)), genome)] + [i.start() for i in re.finditer('(?={0})'.format(re.escape(motif2)), genome)]
    targets_list.sort()

    column_names=["POS", "structure_type"]
    df = pd.DataFrame(columns=column_names)

    for pos in targets_list:
        if pos < 100 or pos > len(genome)-105:
            continue
        print(pos)
        new_dict = r.get_predicted_structure_type_full(outdir, pos, genome_file)
        #print(new_dict)
        new_row = [pos, new_dict["structure_type"]]
        df.loc[len(df)] = new_row

    df.to_csv(os.path.join(outdir, "APOBEC_targets_secondary_structure.txt"), sep='\t', index=False)
    print(df.groupby(by=["structure_type"]).count())

vcf_reader()

#find_targets()
