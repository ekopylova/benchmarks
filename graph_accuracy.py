#!/usr/bin/python

import click
import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import  matplotlib.pyplot as plt
import matplotlib.cm as cm
from os.path import join, exists
import glob
import brewer2mpl


def graph_accuracy(accuracy_fp,
                   output_acc_fp,
                   offset=0,
                   platform="Illumina"):
    """
    """
    mpl.rcParams['legend.fontsize'] = 10
    mpl.rcParams['xtick.major.pad']='10'
    mpl.rcParams['ytick.major.pad']='10'
    mpl.rcParams['grid.linestyle'] = "-"
    mpl.rcParams['grid.color'] = "0.7"
    fig = plt.figure()
    ax = fig.add_subplot(111, axisbg='0.97')

    if offset == 0:
        acc_ind = 10
    elif offset == 5:
        acc_ind = 11
    elif offset == 10:
        acc_ind = 12
    else:
        raise ValueError("%s offset is not allowed" % offset)

    # count number of tools
    tools = []
    with open(accuracy_fp, 'U') as accuracy_f:
        for line in accuracy_f:
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')
            tool = line[0]
            if tool not in tools:
                tools.append(tool)            

    num_colors = len(tools)
    colors = brewer2mpl.get_map('Dark2', 
                                'qualitative', 
                                num_colors).mpl_colors

    # plot F-measure vs. accuracy
    tools = []
    fmeasure_list = []
    accuracy_list = []
    proxy_list = []
    tools_list = []

    with open(accuracy_fp, 'U') as accuracy_f:
        ind = 0
        for line in accuracy_f:
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')
            tool = line[0]
            fmeasure = float(line[9])
            accuracy = round(float(line[acc_ind])/100.0, 3)
            usertime = round(float(line[11]),0)
            walltime = round(float(line[12]),0)
            if tool not in tools:
                # plot existing tool
                if len(fmeasure_list) > 0:
                    ax.scatter(accuracy_list, fmeasure_list, marker='o', s=60, label=tools[ind], color=colors[ind])
                    proxy = plt.Rectangle((0,0), 1, 1, fc=colors[ind])
                    proxy_list.append(proxy)
                    tools_list.append(tools[ind])
                    ind += 1
                tools.append(tool)
                # reset points for next tool
                fmeasure_list = [fmeasure]
                accuracy_list = [accuracy]
            else:
                fmeasure_list.append(fmeasure)
                accuracy_list.append(accuracy)

        # plot last tool
        if len(fmeasure_list) > 0:
            ax.scatter(accuracy_list, fmeasure_list, marker='o', s=60, label=tools[ind], color=colors[ind])
            proxy = plt.Rectangle((0,0), 1, 1, fc=colors[ind])
            proxy_list.append(proxy)
            tools_list.append(tools[ind])

        #ax.xaxis._axinfo['label']['space_factor'] = 1.9
        #ax.yaxis._axinfo['label']['space_factor'] = 1.9
        ax.set_xlabel('Normalized Weighted True-Positive Alignment score [0,1]')
        ax.set_ylabel('F-measure [0,1]')
        zed = [tick.label.set_fontsize(10) for tick in ax.yaxis.get_major_ticks()]
        zed = [tick.label.set_fontsize(10) for tick in ax.xaxis.get_major_ticks()]

        #legend
        #lgd = ax.legend(proxy_list, tools_list, bbox_to_anchor=(1.3, 0.9), borderaxespad=0.5)
        #plt.savefig(output_acc_fp, bbox_extra_artists=(lgd,), bbox_inches='tight', pad_inches=0.8)
        plt.grid()
        plt.savefig(output_acc_fp)



@click.command()
@click.argument('accuracy_fp', required=True,
                type=click.Path(resolve_path=True, readable=True, exists=True,
                                file_okay=True))
@click.argument('output_acc_fp', required=True,
                type=click.Path(resolve_path=True, readable=True, exists=False,
                                file_okay=True))
@click.option('--offset', required=False, type=int, default=0, show_default=True,
              help="Maximum difference between expected alignment position and observed")
@click.option('--platform', required=False, type=str, default='Illumina', show_default=True,
              help="platform can be Illumina, Roche 454 or Ion Torrent PGM")
def _main(accuracy_fp, output_acc_fp, offset, platform):
    """
    """

    if platform not in ["Illumina", "Roche 454", "Ion Torrent PGM"]:
        raise ValueError("%s can only be one of Illumina, Roche 454 or Ion Torrent PGM" % platform)

    graph_accuracy(accuracy_fp=accuracy_fp,
                   output_acc_fp=output_acc_fp,
                   offset=offset,
                   platform=platform)


if __name__ == "__main__":
    _main()
                
                
