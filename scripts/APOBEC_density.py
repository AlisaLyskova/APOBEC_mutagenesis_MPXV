from scipy.stats import wilcoxon
import pandas as pd
import sys
import cyvcf2
import re
from Bio import SeqIO
from Bio.Seq import Seq

def load_genome(fasta_file):
    genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    new_genome = genome["NC_063383.1"].seq
    return str(new_genome)

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

def find_genes(genome, gtf, position, ref, df_out, df_genes, df_repeats, chromosome="NC_063383.1"):
    """Detect the amino acid change caused by a mutation."""
    cds = get_cds(gtf, chromosome, position)
    repeats = df_repeats[(df_repeats['chromStart']<=position) & (df_repeats['chromEnd']>=position)]
    if len(repeats) != 0:
        inSTR = 1
    else:
        inSTR = 0
    
    if cds.empty:
        new_row = [position, ref, "-", "intergenic", inSTR]
        df_out.loc[len(df_out)] = new_row
        return df_out
        
    for row in range(0, len(cds)):
        strand = cds.iloc[row]['strand']
        info = cds.iloc[row]['attribute'].split(';')
        for x in info:
            if x.startswith("gene="):
                gene_name = x.split("=")[1]
                gene_type_list = df_genes[(df_genes["name"] == gene_name) & (df_genes["strand"] == strand)].stage.tolist()
                if len(gene_type_list) != 0:
                    gene_type = gene_type_list[0]
                else:
                    gene_type = "unknown"
                new_row = [position, ref, gene_name, gene_type, inSTR]
                df_out.loc[len(df_out)] = new_row
    return df_out


def find_targets(motif, genome, gtf, genes, repeats, Chr="NC_063383.1"):
    motif2 = str(Seq(motif).reverse_complement())
    #write for TC and GA position of substitution (+1 and 0 accordingly)
    targets_list = [i.start()+1 for i in re.finditer('(?={0})'.format(re.escape(motif)), genome)] + [i.start() for i in re.finditer('(?={0})'.format(re.escape(motif2)), genome)]
    targets_list.sort()
    
    column_names=['pos', 'ref', 'gene_name', 'gene_type', "inRepeat"]
    targets = pd.DataFrame(columns=column_names)

    for pos in targets_list:
        ref = genome[pos]
        targets = find_genes(genome, gtf, pos, ref, targets, genes, repeats)
        
    return targets


def vcf_reader(input_file, genome, gtf, genes, repeats):
    
    column_names=['pos', 'ref', 'gene_name', 'gene_type', "inRepeat"]
    subst = pd.DataFrame(columns=column_names)
    
    vcf_file = cyvcf2.VCF(input_file, gts012=True)
    
    for record in vcf_file:
        REF, POS = record.REF, record.POS - 1
        subst = find_genes(genome, gtf, POS, REF, subst, genes, repeats)
        
    return subst


def APOBEC_density(df_targets, df_APOBEC):
    
    new_data = {'genes_type': ['early', 'intermediate', 'late', 'intergenic'], 'D_APOBEC_in': '-', 'D_APOBEC_out': '-'}
    D_APOBEC_df = pd.DataFrame(new_data)
    D_APOBEC_df = D_APOBEC_df.set_index('genes_type')
    
    for gene_stage in ['early', 'late', 'intermediate', 'intergenic']:
        N_TC_in = len(df_targets[df_targets.gene_type.str.match(gene_stage)])
        N_APOBEC_in = len(df_APOBEC[df_APOBEC.gene_type.str.match(gene_stage)])
        D_APOBEC_in = N_APOBEC_in/N_TC_in
        D_APOBEC_df.at[gene_stage, 'D_APOBEC_in'] = D_APOBEC_in
    
        N_TC_out = len(df_targets) - N_TC_in
        N_APOBEC_out = len(df_APOBEC) - N_APOBEC_in
        D_APOBEC_out = N_APOBEC_out/N_TC_out
        D_APOBEC_df.at[gene_stage, 'D_APOBEC_out'] = D_APOBEC_out

    D_APOBEC_df = D_APOBEC_df.reset_index()
    print(D_APOBEC_df)
    print(wilcoxon(D_APOBEC_df.D_APOBEC_in.tolist(), D_APOBEC_df.D_APOBEC_out.tolist()))
    
    #for repeats
    N_TC_in = len(df_targets[df_targets["inRepeat"] == 1])
    N_APOBEC_in = len(df_APOBEC[df_APOBEC["inRepeat"] == 1])
    D_APOBEC_in = N_APOBEC_in/N_TC_in
    
    N_TC_out = len(df_targets) - N_TC_in
    N_APOBEC_out = len(df_APOBEC) - N_APOBEC_in
    D_APOBEC_out = N_APOBEC_out/N_TC_out
    print("D_APOBEC_in = ", D_APOBEC_in)
    print("D_APOBEC_out = ", D_APOBEC_out)
    print(wilcoxon(D_APOBEC_in, D_APOBEC_out))


VCF1 = sys.argv[1]
VCF_common_pos = sys.argv[2]
VCF_all_pos = sys.argv[3]
GENES = sys.argv[4]
GENOME = sys.argv[5]
GTF_FILE = sys.argv[6]
REPEATS_FILE = sys.argv[7]

genome = load_genome(GENOME)
gtf = load_gtf(GTF_FILE)
genes = pd.read_csv(GENES, sep="\t")
repeats = pd.read_csv(REPEATS_FILE, sep="\t")

df_targets = find_targets("TC", genome, gtf, genes, repeats)
df_APOBEC_sample = vcf_reader(VCF1, genome, gtf, genes, repeats)
APOBEC_density(df_targets, df_APOBEC_sample)
print(len(df_APOBEC_sample[df_APOBEC_sample["inRepeat"] == 1]))
print(df_APOBEC_sample[df_APOBEC_sample["inRepeat"] == 1])

#df_APOBEC_common = vcf_reader(VCF_common_pos, genome, gtf, genes, repeats)
#APOBEC_density(df_targets, df_APOBEC_common)
#print(len(df_APOBEC_common[df_APOBEC_common["inRepeat"] == 1]))

#df_APOBEC_all = vcf_reader(VCF_all_pos, genome, gtf, genes, repeats)
#APOBEC_density(df_targets, df_APOBEC_all)
#print(len(df_APOBEC_all[df_APOBEC_all["inRepeat"] == 1]))
#print(df_APOBEC_all[df_APOBEC_all["inRepeat"] == 1])
