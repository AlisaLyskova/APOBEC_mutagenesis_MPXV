#create cumulitive plot for bam files based on mosdepth statistics
#https://github.com/brentp/mosdepth
##INPUT: <sample>.mosdepth.global.dist.txt files 
##OUTPUT: html plot with coverage on x axes and Proportion of bases at coverage on y axes

import sys
import os
import string
import json
import itertools as it
from operator import itemgetter
import collections
from argparse import ArgumentParser


def main():
    args = get_args()
    traces = collections.defaultdict(list)
    chroms = collections.OrderedDict()
    chroms["total"] = True

    #dictionary with samples features
    ##for PRJEB56841
    #samples_features = {'ERR10550019': '1hA', 'ERR10963107': '1hA', 'ERR10550020': '1hB', 'ERR10963108': '1hB', 'ERR10550021': '1hC', 'ERR10963109': '1hC', 'ERR10550022': '2hA', 'ERR10963110': '2hA', 'ERR10550023': '2hB', 'ERR10963111': '2hB', 'ERR10550024': '2hC', 'ERR10963112': '2hC', 'ERR10550025': '4hA', 'ERR10963113': '4hA', 'ERR10550026': '4hB', 'ERR10963114': '4hB', 'ERR10550027': '4hC', 'ERR10963115': '4hC', 'ERR10550028': '6hA', 'ERR10963116': '6hA', 'ERR10550029': '6hB', 'ERR10963117': '6hB', 'ERR10550030': '6hC', 'ERR10963118': '6hC', 'ERR10543479': '12hA', 'ERR10963119': '12hA', 'ERR10543480': '12hB', 'ERR10963120': '12hB', 'ERR10543481': '12hC', 'ERR10963121': '12hC', 'ERR10550031': '24hA', 'ERR10963122': '24hA', 'ERR10550032': '24hB', 'ERR10963123': '24hB', 'ERR10550033': '24hC', 'ERR10963124': '24hC', 'ERR10550034': 'mockA', 'ERR10963125': 'mockA', 'ERR10550035': 'mockB', 'ERR10963126': 'mockB', 'ERR10550036': 'mockC', 'ERR10963127': 'mockC', 'ERR10513574': 'dRNA', 'ERR10963128': 'dRNA'}
    #samples_features={}
    ##PRJEB60728
    #samples_features={'ERR11026612': '12h', 'ERR11026613': '6h', 'ERR11026615': '12h', 'ERR11026616': '24h', 'ERR11026617': '12h', 'ERR11026625': '24h', 'ERR11026630': '12h', 'ERR11026631': '6h', 'ERR11026635': '24h', 'ERR11026636': '6h', 'ERR11026637': '24h', 'ERR11026638': '6h', 'ERR11030164': '12h', 'ERR11030165': '6h', 'ERR11030167': '12h', 'ERR11030168': '24h', 'ERR11030169': '12h', 'ERR11030177': '24h', 'ERR11030182': '12h', 'ERR11030183': '6h', 'ERR11030187': '24h', 'ERR11030188': '6h', 'ERR11030189': '24h', 'ERR11030190': '6h'}

    for f in args.input:
        file_name = os.path.basename(f)
        sample = file_name.replace(".mosdepth.global.dist.txt", "")
        
        #sample = sample+"_"+samples_features[sample]
        
        gen = (x.rstrip().split("\t") for x in open(f))
        for chrom, data in it.groupby(gen, itemgetter(0)):
            if chrom.startswith("GL"):
                continue
            if "Un" in chrom: continue
            if "random" in chrom or "HLA" in chrom: continue
            if chrom.endswith("alt"): continue
            chroms[chrom] = True
            xs, ys = [], []
            v50 = 0
            found = False
            for _, x, y in data:
                y = float(y)
                if y < 0.01:
                    continue
                if not found and y > 0.5:
                    v50 = x
                    found = True
                    print("{}\t{}\t{}\t{:.3f}".format(sample, chrom, x, y))

                xs.append(float(x))
                ys.append(y)

            if len(xs) > 100:
                xs = [x for i, x in enumerate(xs) if ys[i] > 0.02]
                ys = [y for y in ys if y > 0.02]
                if len(xs) > 100:
                    xs = xs[::2]
                    ys = ys[::2]

            traces[chrom].append({
                'x': [round(x, 3) for x in xs],
                'y': [round(y, 3) for y in ys],
                'mode': 'lines',
                'name': sample + (" (%.1f)" % float(v50))
            })
    #select one chromosome
    chroms = ['NC_063383.1']

    tmpl = """<html>
    <head>
      <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
      <style>
    div {
      width: 1200px;
      height: 700px;
    }
    </style
    </head>
    <body>$plot_divs</div>
    <script>
    var layout = {
        hovermode: 'closest',
        xaxis: {title: 'Coverage'},
        yaxis: {title: 'Proportion of bases at coverage', domain: [0, 1], dtick: 0.25},
        showlegend: $showlegend,
        autosize: true,
        legend: {
            x: 1,
            y: 1
        },
    }
    """
    footer = """
    </script>
    </body>
    </html>"""

    chr_tmpl = """
    Plotly.newPlot('plot-div-$chrom', $data, layout, {displayModeBar: false, displaylogo: false, fillFrame: false, autosizeable: true});
    """

    tmpl = string.Template(tmpl)
    try:
        with open(args.output, "w") as html:
            divs = "\n".join("<{div}>{chrom}</{div}><div id='plot-div-{chrom}'></div><hr/>".format(
                chrom=c, div="h2" if c == "total" else "b") for c in chroms)
            html.write(tmpl.substitute(showlegend="true" if len(
                sys.argv[1:]) < 50 else "false", plot_divs=divs))
            for chrom in chroms:
                html.write(string.Template(chr_tmpl).substitute(
                    chrom=chrom, data=json.dumps(traces[chrom])))
            html.write(footer)
    except FileNotFoundError:
        sys.exit("ERROR: failed creating output file, does the path exist?")


def get_args():
    parser = ArgumentParser(description="Creates html plots from mosdepth results.")
    parser.add_argument("-o", "--output",
                        default="dist.html",
                        help="path and name of output file. Directories must exist.")
    parser.add_argument("input",
                        nargs='+',
                        help="the dist file(s) to use for plotting")
    return parser.parse_args()


if __name__ == '__main__':
    main()

