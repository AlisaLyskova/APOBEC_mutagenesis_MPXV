import pandas as pd
from scipy.stats import pearsonr
import numpy as np
import math
import os
import sys
import argparse


parser = argparse.ArgumentParser(prog="'script_for_motives_specificity.py'", description="Python script for specificity matrix construction and comparison.")
parser.add_argument('-i', '--input', help="Enter the whole path to input TXT file with list of positions", type=str)
parser.add_argument('-g', '--genome', help="Enter the whole path to genome FASTA file.", type=str)
parser.add_argument('-m', '--matrix', help='Enter the whole path to literature matrix TXT file.', type=str)
#parser.add_argument('-lf', '--len_frame', help="Choose the length of the frame: 11 (-5 to +5), 7 (-3 to 3), 5 (-2 to 2), 3 (-1 to 1).", type=str)
parser.add_argument('-o', '--output', help="Enter the whole path to output TXT file with user matrix.", type=str)

args = parser.parse_args()
input_file = args.input
genome_file = args.genome
matrix_file = args.matrix
#frame = args.len_frame
output_file = args.output

# Get list of positions #
df = pd.read_csv(input_file, dtype={'pos':int})
list_of_pos = df['pos'].tolist()

# Get sequence of genome #
with open(genome_file, 'r') as file:
    data = file.readlines()[1:]
    genome_seq = ''.join([s.strip() for s in data])

# Define frame #
frame = ['-5', '-4', '-3', '-2', '-1', '0', '+1', '+2', '+3', '+4', '+5']
'''
if frame == '11':
    frame = ['-5', '-4', '-3', '-2', '-1', '0', '+1', '+2', '+3', '+4', '+5']
elif frame == '7':
    frame = ['-3', '-2', '-1', '0', '+1', '+2', '+3']
elif frame == '5':
    frame = ['-2', '-1', '0', '+1', '+2']
elif frame == '3':
    frame = ['-1', '0', '1']
else:
    print('Enter the length of the frame correctly!')
    sys.exit()
'''
### STEP 1 - USER MATRIX CONSTRUCTION ###
# Function to create user matrix #
def create_matrix(genome, list_pos, frame):
    def get_reverse_complement(seq):
        new_seq = ''
        for n in seq:
            if n == 'A':
                new_seq += 'T'
            elif n == 'C':
                new_seq += 'G'
            elif n == 'G':
                new_seq += 'C'
            else:
                new_seq += 'A'
        return new_seq[::-1]
    
    
    # Matrix construction #
    data = pd.DataFrame(index=['A', 'C', 'G', 'T'], columns=frame, data=0)
    len_frame = len(frame)
    for p in list_pos:
        subseq = genome[(p-(math.floor(len_frame / 2)+1)):(p + math.floor(len_frame / 2))]
        mutation_pos = genome[p - 1]
        if mutation_pos in ['A', 'T']:
            print(f"Position {p} is {mutation_pos} nucleotide in genome! It will not be taken into account when constructing the matrix.")
        elif mutation_pos == 'G':
            subseq = get_reverse_complement(subseq)
        
        for s, i in zip(subseq, frame):
            data.loc[s, i] += 1

    # Normalize matrix #    
    for i in frame:
        data[i] = round(data[i] / data[i].sum(), 3)
    
    return data

user_matrix = create_matrix(genome_seq, list_of_pos, frame)
user_matrix.to_csv(output_file, sep="\t")

### STEP 2 - MATRIX COMPARISON ###
literature_matrix = pd.read_csv(matrix_file, sep='\t', index_col=0)
def compare_two_matrix_by_Pearson(m1, m2):
    m1 = m1[m2.columns]
    vm1 = m1.values.flatten(order='F') # to vector by column to column 
    vm2 = m2.values.flatten(order='F') # to vector by column to column
    comparison_result = pearsonr(vm1, vm2)
    
    return comparison_result

print(f"Answer:\n{compare_two_matrix_by_Pearson(user_matrix, literature_matrix)}")
