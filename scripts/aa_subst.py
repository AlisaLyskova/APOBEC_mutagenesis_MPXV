import pandas as pd
import numpy as np
import sys
import cyvcf2
import os
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable


VCF = sys.argv[1]
OUTFILE_APOBEC = sys.argv[2]
OUTFILE_TARGETS = sys.argv[3]
GENOME = sys.argv[4]
GTF_FILE = sys.argv[5]


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

def detect_mutation_impact(genome, gtf, position, ref, alt, df_out, chromosome="NC_063383.1"):
    """Detect the amino acid change caused by a mutation."""
    cds = get_cds(gtf, chromosome, position)
    
    if cds.empty:
        new_row = [position, ref, alt, "-", "-", "-", "-", "-", "intergenic", "-"]
        df_out.loc[len(df_out)] = new_row
        return df_out
    for row in range(0, len(cds)):
        
        strand = cds.iloc[row]['strand']
        gene_start = cds.iloc[row]['start']

        # Calculate the codon position
        codon_start = gene_start + ((position - gene_start) // 3) * 3

        # Extract the codon
        codon_seq = genome[codon_start-1:codon_start+2]
        if strand == '-':
            codon_seq = codon_seq.reverse_complement()
            ref = str(Seq(ref).reverse_complement())
            alt = str(Seq(alt).reverse_complement())

        # Translate the original codon
        original_amino_acid = translate_codon(codon_seq)

        # Create the mutated codon
        mutation_index = (position - codon_start) % 3
        if strand == '-':
            if mutation_index == 0:
                mutation_index = 2
            elif mutation_index == 2:
                mutation_index = 0

        mutated_codon = list(codon_seq)
        mutated_codon[mutation_index] = alt
        mutated_codon = Seq(''.join(mutated_codon))

        # Translate the mutated codon
        mutated_amino_acid = translate_codon(mutated_codon)

        nucl_pos = mutation_index+1

        if mutated_amino_acid == original_amino_acid:
            category = "synonymous"
        else:
            if mutated_amino_acid == "*":
                category = "nonsense"
            else:
                category = "nonsynonymous"

        new_row = [position, ref, alt, nucl_pos, codon_seq, original_amino_acid, mutated_codon, mutated_amino_acid, category, strand]
        df_out.loc[len(df_out)] = new_row
    
    return df_out


def find_targets(outfile=OUTFILE_TARGETS, genome_file=GENOME, gtf_file=GTF_FILE, motif="TC", Chr="NC_063383.1"):

    genome = load_genome(genome_file)
    gtf = load_gtf(gtf_file)

    motif2 = str(Seq(motif).reverse_complement())
    #write for TC and GA position of substitution (+1 and 0 accordingly, starting with 1 like in vcf file)
    targets_list = [i.start()+2 for i in re.finditer('(?={0})'.format(re.escape(motif)), str(genome))] + [i.start()+1 for i in re.finditer('(?={0})'.format(re.escape(motif2)), str(genome))]
    targets_list.sort()

    column_names=['POS', 'REF', 'ALT', 'nucl_pos_codon', 'codon', "aa", "mutated_codon", "mutated_aa", "mutation_category", "strand"]
    targets = pd.DataFrame(columns=column_names)

    for pos in targets_list:
        ref = genome[pos-1]
        if ref == "C":
            alt = "T"
        if ref == "G":
            alt = "A"
        targets = detect_mutation_impact(genome, gtf, pos, ref, alt, targets)

    targets.to_csv(outfile, sep="\t", index=False)
    return


def vcf_reader(input_file=VCF, outfile=OUTFILE_APOBEC, genome_file=GENOME, gtf_file=GTF_FILE, chromosome="NC_063383.1"):

    genome = load_genome(genome_file)
    gtf = load_gtf(gtf_file)

    #table with amino acids changes
    column_names=['POS', 'REF', 'ALT', 'nucl_pos_codon', 'codon', "aa", "mutated_codon", "mutated_aa", "mutation_category", "strand"]
    df = pd.DataFrame(columns=column_names)

    vcf_file = cyvcf2.VCF(input_file, gts012=True)

    for record in vcf_file:
        REF, POS, ALT = record.REF, record.POS, record.ALT[0]
        df = detect_mutation_impact(genome, gtf, POS, REF, ALT, df)

    df.to_csv(outfile, sep="\t", index=False)
    return


find_targets()
vcf_reader()
