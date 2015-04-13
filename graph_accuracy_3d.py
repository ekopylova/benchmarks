#!/usr/bin/python

import click
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
mpl.use('Agg')
import numpy as np
import  matplotlib.pyplot as plt
import matplotlib.cm as cm
from os.path import join, exists
import glob
import brewer2mpl


def graph_accuracy(accuracy_fp,
                   output_acc_fp,
                   output_walltime_fp,
                   output_usertime_fp,
                   offset=0,
                   time="Walltime",
                   platform="Illumina"):
    """
    """
    mpl.rcParams['legend.fontsize'] = 10
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

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

    # plot F-measure vs. accuracy vs. usertime
    tools = []
    fmeasure_list = []
    accuracy_list = []
    usertime_list = []
    walltime_list = []
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
            # 10: offset 0; 11: offset 5; 12: offset 10
            accuracy = round(float(line[acc_ind])/100.0, 3)
            usertime = round(float(line[13]),0)
            walltime = round(float(line[14]),0)

            if tool not in tools:
                # plot existing tool
                if len(fmeasure_list) > 0:
                    if time == "Walltime":
                        time_list = walltime_list
                    else:
                        time_list = usertime_list
                    ax.scatter(accuracy_list, fmeasure_list, time_list, marker='o', s=60, label=tools[ind], color=colors[ind])
                    proxy = plt.Rectangle((0,0), 1, 1, fc=colors[ind])
                    proxy_list.append(proxy)
                    tools_list.append(tools[ind])
                    ind += 1
                tools.append(tool)
                # reset points for next tool
                fmeasure_list = [fmeasure]
                accuracy_list = [accuracy]
                walltime_list = [walltime]
                usertime_list = [usertime]
            else:
                fmeasure_list.append(fmeasure)
                accuracy_list.append(accuracy)
                walltime_list.append(walltime)
                usertime_list.append(usertime)

        # plot last tool
        if len(fmeasure_list) > 0:
            if time == "Walltime":
                time_list = walltime_list
            else:
                time_list = usertime_list
            ax.scatter(accuracy_list, fmeasure_list, time_list, marker='o', s=60, label=tools[ind], color=colors[ind])
            proxy = plt.Rectangle((0,0), 1, 1, fc=colors[ind])
            proxy_list.append(proxy)
            tools_list.append(tools[ind])

        ax.xaxis._axinfo['label']['space_factor'] = 1.9
        ax.yaxis._axinfo['label']['space_factor'] = 1.9
        ax.zaxis._axinfo['label']['space_factor'] = 1.9
        ax.set_xlabel('Normalized Weighted True-Positive Alignment score [0,1]')
        ax.set_ylabel('F-measure [0,1]')
        ax.set_zlabel('%s (sec)' % time)
        zed = [tick.label.set_fontsize(10) for tick in ax.yaxis.get_major_ticks()]
        zed = [tick.label.set_fontsize(10) for tick in ax.xaxis.get_major_ticks()]
        zed = [tick.label.set_fontsize(10) for tick in ax.zaxis.get_major_ticks()]
        ax.set_zlim(0)

        #legend
        lgd = ax.legend(proxy_list, tools_list, bbox_to_anchor=(0.1, 0.9), borderaxespad=0.5)
        plt.savefig(output_acc_fp, bbox_extra_artists=(lgd,), bbox_inches='tight', pad_inches=0.8)
        #plt.show()


@click.command()
@click.argument('accuracy_fp', required=True,
                type=click.Path(resolve_path=True, readable=True, exists=True,
                                file_okay=True))
@click.argument('output_acc_fp', required=True,
                type=click.Path(resolve_path=True, readable=True, exists=False,
                                file_okay=True))
@click.argument('output_walltime_fp', required=False,
                type=click.Path(resolve_path=True, readable=True, exists=False,
                                file_okay=True))
@click.argument('output_usertime_fp', required=False,
                type=click.Path(resolve_path=True, readable=True, exists=False,
                                file_okay=True))
@click.option('--offset', required=False, type=int, default=0, show_default=True,
              help="Maximum difference between expected alignment position and observed")
@click.option('--time', required=False, type=str, default='Walltime', show_default=True,
              help="z-axis is Walltime or Usertime")
@click.option('--platform', required=False, type=str, default='Illumina', show_default=True,
              help="platform can be Illumina, Roche 454 or Ion Torrent PGM")
def _main(accuracy_fp, output_acc_fp, output_walltime_fp, output_usertime_fp, offset, time, platform):
    """
    """

    if time not in ["Walltime", "Usertime"]:
        raise ValueError("%s is not an option of --time" % time)

    if platform not in ["Illumina", "Roche 454", "Ion Torrent PGM"]:
        raise ValueError("%s can only be one of Illumina, Roche 454 or Ion Torrent PGM" % platform)

    graph_accuracy(accuracy_fp=accuracy_fp,
                   output_acc_fp=output_acc_fp,
                   output_walltime_fp=output_walltime_fp,
                   output_usertime_fp=output_usertime_fp,
                   offset=offset,
                   time=time,
                   platform=platform)


if __name__ == "__main__":
    _main()
                
                
