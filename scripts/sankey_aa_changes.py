from dash import Dash, dcc, html, Input, Output

import json, urllib
import pandas as pd
import sys
import plotly.graph_objects as go

def sankey_plot(INFILE, OUTFILE):
    df = pd.read_csv(INFILE, sep=',')
    df = df.reset_index()
    
    df_grouped = df.groupby(['grantham_rank'])['index'].count()
    print(df_grouped)
    df['aa_v2'] = [x+'_observed' for x in df.aa.tolist()]
    df['mutated_aa_v2'] = [x+'_editing' for x in df.mutated_aa.tolist()]
    df1 = df.groupby(['aa_v2','mutated_aa_v2', 'grantham_rank_color'])['index'].count().reset_index()
    df1.columns = ['source','target', 'grantham_rank_color', 'value']
    df1 = df1.sort_values(by='value', ascending=True)

    unique_source_target = list(pd.unique(df1[['source','target']].values.ravel('k')))

    unique_source_target_label = [list(x)[0] for x in unique_source_target]

    mapping_dict = {k: v for v, k in enumerate(unique_source_target)}

    df1['source'] = df1['source'].map(mapping_dict)
    df1['target'] = df1['target'].map(mapping_dict)

    new_dict = df1.to_dict(orient='list')


    fig = go.Figure(data=[go.Sankey(
        node = dict(
            pad = 20,
            thickness=30,
            line=dict(color='black', width=0.5),
            label = unique_source_target_label,
            color='#BABBB6'
        ),
        link = dict(
            source= new_dict['source'],
            target = new_dict['target'],
            value = new_dict['value'],
            color = new_dict['grantham_rank_color']
        )

    )
    ])
    #legend={"yellow": "conservative", "pink": "moderately conservative", "dark blue": "moderately radical", "brown": "radical", "green": "synonymous"}
    fig.update_layout(
        font_family="Arial",
        font_color="black",
        font_size=18,
        width=700,
        height=1000
    )

    fig.write_html(OUTFILE)


sankey_plot(sys.argv[1], sys.argv[2])
