import pandas as pd
import numpy as np
import sys
import cyvcf2
import os
from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns


DIR = sys.argv[1]
OUTDIR = sys.argv[2]
GENOME = sys.argv[3]
GTF_FILE = sys.argv[4]

colors = ['#CADCFC', '#8AB6F9', '#00246B']
my_cmap = LinearSegmentedColormap.from_list("custom", colors, N=100)


def load_genome(fasta_file):
    genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    new_genome = genome["NC_063383.1"].seq
    return new_genome

def load_gtf(gtf_file):
    """Load the GTF file into a pandas DataFrame."""
    gtf_columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    gtf = pd.read_csv(gtf_file, sep="\t", comment='#', header=None, names=gtf_columns)
    return gtf

def get_cds(gtf, chromosome, position):
    """Find the CDS where the mutation is located."""
    cds_rows = gtf[(gtf['seqname'] == chromosome) & 
                   (gtf['feature'] == 'CDS') & 
                   (gtf['start'] <= position) & 
                   (gtf['end'] >= position)]
    return cds_rows

def translate_strand(genome, gtf, position, ref, alt, chromosome="NC_063383.1"):
    """Detect the amino acid change caused by a mutation."""
    cds = get_cds(gtf, chromosome, position)
    
    if cds.empty:
        ref_norm = ref
        alt_norm = alt
    else:
        strand = cds.iloc[0]['strand']
        if strand == '-':
            ref_norm = str(Seq(ref).reverse_complement())
            alt_norm = str(Seq(alt).reverse_complement())
        else:
            ref_norm = ref
            alt_norm = alt

    return ref_norm+">"+alt_norm



def main(directory=DIR, outdir=OUTDIR, genome_file=GENOME, gtf_file=GTF_FILE, chromosome="NC_063383.1", colors=my_cmap):

    genome = load_genome(genome_file)
    gtf = load_gtf(gtf_file)

    samples_list = []
    for filename in os.listdir(directory):
        filepath = os.path.join(directory, filename)
        # checking if it is a file
        if filepath.endswith(".vcf.gz"):
            sample = filename.split('_')[0]
            samples_list.append(sample)

    df = pd.DataFrame(columns=['substitution']+samples_list)
    df['substitution'] = ["C>T", "G>A"]
    df = df.set_index('substitution')
    df.fillna(0, inplace=True)

    df_norm = df.copy()

    for filename in os.listdir(directory):
        filepath = os.path.join(directory, filename)
        # checking if it is a file
        if filepath.endswith(".vcf.gz"):
            print("Process: ", filepath)
            sample = filename.split('_')[0]

            vcf_file = cyvcf2.VCF(filepath, gts012=True)
            for record in vcf_file:
                REF, ALT, POS = record.REF, record.ALT[0], record.POS - 1
                subst=REF+">"+ALT

                subst_norm = translate_strand(genome, gtf, POS, REF, ALT)

                df.at[subst, sample] += 1
                df_norm.at[subst_norm, sample] += 1


    df.reset_index(inplace=True)
    df = df.set_index(['substitution'])
    df = df[samples_list].stack()
    df = df.to_frame().reset_index()
    df = df.dropna()
    df = df.rename(columns={"level_1": "sample", 0:'Nsubst'})
    df = df.reset_index(drop=True)
    df = df.pivot(columns="substitution", index="sample", values="Nsubst")
    print(df)
    #sns.set_context("notebook", font_scale=0.5, rc={"figure.figsize": (8, 10)})
    #sns.heatmap(df, cmap=colors, annot=True, fmt='g')
    #outfile1 = os.path.join(outdir, 'APOBEC_subst_not_translated.png')
    #plt.savefig(outfile1, dpi=800)
    #plt.close()

    df_norm.reset_index(inplace=True)
    df_norm = df_norm.set_index(['substitution'])
    df_norm = df_norm[samples_list].stack()
    df_norm = df_norm.to_frame().reset_index()
    df_norm = df_norm.dropna()
    df_norm = df_norm.rename(columns={"level_1": "sample", 0:'Nsubst'})
    df_norm = df_norm.reset_index(drop=True)
    df_norm = df_norm.pivot(columns="substitution", index="sample", values="Nsubst")
    print(df_norm)
    #sns.set_context("notebook", font_scale=0.5, rc={"figure.figsize": (8, 10)})
    #sns.heatmap(df_norm, cmap=colors, annot=True, fmt='g')
    #outfile2 = os.path.join(outdir, 'APOBEC_subst_translated.png')
    #plt.savefig(outfile2, dpi=800)
    #plt.close()

    return

main()
