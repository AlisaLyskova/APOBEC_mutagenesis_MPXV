import pandas as pd
import numpy as np
import sys
import cyvcf2
from Bio import SeqIO
from Bio.Seq import Seq


VCF = sys.argv[1]
GENOME = sys.argv[2]
GTF_FILE = sys.argv[3]
TSS_FILE = sys.argv[4]
TES_FILE = sys.argv[5]


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

def translate_strand(genome, gtf, position, ref, chromosome="NC_063383.1"):
    """Detect the amino acid change caused by a mutation."""
    cds = get_cds(gtf, chromosome, position)

    if cds.empty:
        ref_norm = ref
    else:
        strand = cds.iloc[0]['strand']
        if strand == '-':
            ref_norm = str(Seq(ref).reverse_complement())
        else:
            ref_norm = ref

    return ref_norm



def main(filepath=VCF, genome_file=GENOME, gtf_file=GTF_FILE, tss_table=TSS_FILE, tes_table=TES_FILE, chromosome="NC_063383.1"):

    genome = load_genome(genome_file)
    gtf = load_gtf(gtf_file)

    df_GA = pd.DataFrame(columns=['POS', 'REF', 'REF_NORM'])

    vcf_file = cyvcf2.VCF(filepath, gts012=True)
    for record in vcf_file:
        REF, POS = record.REF, record.POS - 1

        ref_norm = translate_strand(genome, gtf, POS, REF)
        if ref_norm == "G":
            df_GA.loc[len(df_GA)] = [POS, REF, ref_norm]

    TSS = pd.read_excel(tss_table)
    TES = pd.read_excel(tes_table)

    #filter like in the article
    TSS = TSS[(TSS['dRNA_in_25']==True) & (TSS['count']>=10) & (TSS['promoter_in_40']==True)]
    TSS = TSS.reset_index(drop=True)
    TES = TES[(TES['PAS_in_50']==True) & (TES['count']>=3)]
    TES = TES.reset_index(drop=True)


    df_GA['transcript'] = '-'
    for row in range(0, len(df_GA)):
        pos = int(df_GA.loc[row, 'POS'])
        ref = df_GA.loc[row, 'REF']
        if ref == 'C':
            TSS_filtered = TSS[(TSS['strand']=='+') & (TSS['end'] <= pos)]
            TSS_filtered = TSS_filtered.reset_index(drop=True)
            if len(TSS_filtered) != 0:
                TSS_coordinate = TSS_filtered['start'].iloc[-1]
                #if there is TES in [TSS:pos] than pos is not transcribed
                if len(TES[(TES['strand']=='+') & (TES['start'] > TSS_coordinate) & (TES['end'] < pos)]) == 0:
                    TES_filtered = TES[(TES['strand']=='+') & (TES['end'] >= pos)]
                    TES_filtered = TES_filtered.reset_index(drop=True)
                    TES_coordinate = TES_filtered.loc[0, 'start']
                    df_GA.at[row, 'transcript'] = str(TSS_coordinate)+","+str(TES_coordinate)
        if ref == 'G':
            TSS_filtered = TSS[(TSS['strand']=='-') & (TSS['end'] >= pos)]
            TSS_filtered = TSS_filtered.reset_index(drop=True)
            if len(TSS_filtered) != 0:
                TSS_coordinate = TSS_filtered.loc[0, 'end']
                #if there is TES in [pos:TSS] than pos is not transcribed
                if len(TES[(TES['strand']=='-') & (TES['end'] < TSS_coordinate) & (TES['start'] > pos)]) == 0:
                    TES_filtered = TES[(TES['strand']=='-') & (TES['start'] <= pos)]
                    TES_filtered = TES_filtered.reset_index(drop=True)
                    TES_coordinate = TSS_filtered['end'].iloc[-1]
                    df_GA.at[row, 'transcript'] = str(TSS_coordinate)+","+str(TES_coordinate)

    print(len(df_GA))
    print("Transcription is possible in {} positions".format(len(df_GA[df_GA['transcript']!= '-'])))
    return

main()
