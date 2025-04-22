import pandas as pd
import numpy as np
import sys
import cyvcf2
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable

##min variant frequency
min_VF = 0.01

##vcf file
VCF = sys.argv[1]
OUT = sys.argv[2]
GENOME = sys.argv[3]
GTF_FILE = sys.argv[4]


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

def translate_codon(codon):
    """Translate a codon to its corresponding amino acid."""
    return str(Seq(codon).translate(table=CodonTable.unambiguous_dna_by_id[11]))

def detect_mutation_impact(genome, gtf, chromosome, position, ref, alt):
    """Detect the amino acid change caused by a mutation."""
    cds = get_cds(gtf, chromosome, position)
    
    motif3 = genome[position-1:position+2]
    if cds.empty:
        motif3_norm = motif3
        ref_norm = ref
        alt_norm = alt
    else:    
        strand = cds.iloc[0]['strand']
        if strand == '-':
            motif3_norm = motif3.reverse_complement()
            ref_norm = str(Seq(ref).reverse_complement())
            alt_norm = str(Seq(alt).reverse_complement())
        else:
            motif3_norm = motif3
            ref_norm = ref
            alt_norm = alt

    return motif3, motif3_norm, ref_norm, alt_norm




def preprocessing(input_file=VCF, output_file=OUT, genome_file=GENOME, gtf_file=GTF_FILE, chromosome="NC_063383.1"):

    column_names=['pos', 'ref', 'alt', 'motif3', 'isAPOBEC', 'isAPOBECextended', 'vaf', 'altcnt', 'ref_norm', 'alt_norm', 'motif3_norm']
    df = pd.DataFrame(columns=column_names)

    vcf_file = cyvcf2.VCF(input_file, gts012=True)

    genome = load_genome(genome_file)
    gtf = load_gtf(gtf_file)
    
    for record in vcf_file:
        REF, ALT, POS = record.REF, record.ALT[0], record.POS - 1
        vaf = record.format('VAF').item()
        ad = record.format('AD')
        Nreads_alt = ad[0][1]

        motif3, motif3_norm, ref_norm, alt_norm = detect_mutation_impact(genome, gtf, chromosome, POS, REF, ALT)

        isAPOBEC = 0
        isAPOBECextended = 0
        if REF == "C":
            prev_nucl = genome[POS-1]
            if prev_nucl != "T":
                isAPOBECextended = 2
            else:
                isAPOBEC = 1
                isAPOBECextended = 1
        if REF == "G":
            prev_nucl = genome[POS+1]
            if prev_nucl != "A":
                isAPOBECextended = 2
            else:
                isAPOBEC = 1
                isAPOBECextended = 1

        new_row = [POS, REF, ALT, motif3, isAPOBEC, isAPOBECextended, vaf, Nreads_alt, ref_norm, alt_norm, motif3_norm]
        df.loc[len(df)] = new_row
    df[['pos', 'ref_norm', 'alt_norm', 'motif3_norm', 'isAPOBEC', 'isAPOBECextended', 'vaf', 'altcnt']].to_csv(output_file, sep='\t', index=False, header=False)
    #df.to_csv(output_file, sep='\t', index=False)


preprocessing()

