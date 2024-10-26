import pandas as pd
import numpy as np
import sys
import os
import gzip
from dash import Dash, dcc, html
from plotly import graph_objects as go

#input - directory with samples and mosdepth results
directory = sys.argv[1]

count = 0
color_list = ['#C6878F', '#3D405B', '#81B29A', '#E07A5F', '#5F797B', '#F2CC8F']

app = Dash(__name__)
fig = go.Figure()

for filename in os.listdir(directory):
    f = os.path.join(directory, filename)
    if f.endswith(".regions.bed.gz"):
        print(filename)
        df = pd.read_csv(f, compression='gzip', names=["#chrom", 'start', "end", "coverage"], sep="\t")
        fig.add_trace(go.Scatter(
            name=filename.split("_")[0], 
            x=df.start.tolist(), 
            y=df.coverage.tolist(), 
            offsetgroup=count, 
            mode='lines',
            opacity=0.8,
            marker=dict(color=color_list[count])
        ))
        count += 1


fig.update_layout(title="Coverage",
            xaxis_title="Position",
            yaxis_title="Coverage")

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

fig.update_yaxes(range=[0, 10000])

fig.write_html(os.path.join(directory, "coverage.html"))
