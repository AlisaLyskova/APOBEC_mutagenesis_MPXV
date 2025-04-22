import pandas as pd
import numpy as np
import sys
import os
import plotly.express as px
import cyvcf2

#caller_name: clair3, lofreq, bcftools
##change color_list
def VAF_distribution(input_file):
    file_name = os.path.basename(input_file)
    caller = file_name.split("_")[1]
    feature = file_name.split("_")[2]
    vcf_type = caller + "_" + feature

    vcf_file = cyvcf2.VCF(input_file, gts012=True)
    
    VAF_list = []
    
    for record in vcf_file:
        if caller == "bcftools":
            ad = record.format('AD')
            ad1 = ad[0][0]
            ad2 = ad[0][1]
            #vaf = record.format('VAF').item()
            #dp = ad1 + ad2

        if caller == "lofreq":
            dp4 = record.INFO.get('DP4')
            ad1 = dp4[0]+dp4[1]
            ad2 = dp4[2]+dp4[3]
            #dp = record.INFO.get('DP')
            #vaf = record.INFO.get('AF')

        if caller == "clair3":
            ad = record.format('AD')
            ad1 = ad[0][0]
            ad2 = ad[0][1]
            #vaf = record.format('VAF').item()
            #dp = record.format('DP').item()

        ref_depth, alt_depth = int(ad1), int(ad2)
        VAF_value = alt_depth / (ref_depth + alt_depth) if (ref_depth + alt_depth) != 0 else 0
        VAF_list.append(VAF_value)


    # Define bin edges using linspace
    bin_edges = np.linspace(0, 1, 50)

    # Bin the data using digitize
    bin_indices = np.digitize(VAF_list, bin_edges, right=True)

    # Calculate histogram counts using bin count
    hist = np.bincount(bin_indices)
    hist = [(i/len(VAF_list))*100 for i in hist]

    return [(round(bin, 4),round(val, 4)) for bin,val in zip(bin_edges, hist)]


directory=sys.argv[1]
outfile=sys.argv[2]
#caller=sys.argv[3]

x_list = []
y_list = []
samples = []
n_samples = 0
for filename in os.listdir(directory):
    f = os.path.join(directory, filename)
    if f.endswith(".vcf.gz"):
        res = VAF_distribution(f)
        file_name = os.path.basename(f)
        sample = file_name.split("_")[0]
        caller = file_name.split("_")[1]
        feature = file_name.split("_")[2]
        #vcf is empty
        if not res:
            continue
        x, y = zip(*res)
        x_list += x
        y_list += y
        #samples += [sample+'_'+caller+'_'+feature]*len(x)
        samples += [sample]*len(x)
        n_samples +=1

print(n_samples)
data = {'sample':samples, 'x':x_list, 'y':y_list}
#data = {'sample':samples, 'x':y_list, 'y':x_list}
df = pd.DataFrame.from_dict(data)

#color_list=[]
#color_discrete_sequence=color_list,

fig = px.bar(df, x='x', y='y', color='sample', category_orders={'sample':sorted(list(set(samples)))}, barmode='overlay', color_discrete_sequence=px.colors.qualitative.Bold, labels={'x': 'Аллельная частота (VAF)', 'y':'Процент'})
fig.update_layout(
    plot_bgcolor='white',
    font_size=18
)
fig.update_xaxes(
    mirror=True,
    ticks='outside',
    showline=True,
    linecolor='black',
    gridcolor='lightgrey'
)
fig.update_yaxes(
    mirror=True,
    ticks='outside',
    showline=True,
    linecolor='black',
    gridcolor='lightgrey'
)
fig.show()
fig.write_html(outfile)
