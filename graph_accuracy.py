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
                   output_walltime_fp,
                   output_usertime_fp):
    """
    """
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

    colors = cm.rainbow(np.linspace(0,1,len(tools)))
    colors = np.random.permutation(colors)

    # plot F-measure vs. accuracy
    tools = []
    fmeasure_list = []
    accuracy_list = []
    usertime_list = []
    walltime_list = []

    with open(accuracy_fp, 'U') as accuracy_f:
        ind = 0
        for line in accuracy_f:
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')
            tool = line[0]
            fmeasure = float(line[9])
            accuracy = round(float(line[10])/100.0, 3)
            usertime = round(float(line[11]),0)
            walltime = round(float(line[12]),0)
            if tool not in tools:
                # plot existing tool
                if len(fmeasure_list) > 0:
                    plt.figure(0)
                    plt.scatter(accuracy_list, fmeasure_list, marker='o', label=tools[ind], color=colors[ind])
                    plt.figure(1)
                    plt.scatter(accuracy_list, walltime_list, marker='o', label=tools[ind], color=colors[ind])
                    plt.figure(2)
                    plt.scatter(accuracy_list, usertime_list, marker='o', label=tools[ind], color=colors[ind])
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
            plt.figure(0)
            plt.scatter(accuracy_list, fmeasure_list, marker='o', label=tools[ind], color=colors[ind])
            plt.figure(1)
            plt.scatter(accuracy_list, walltime_list, marker='o', label=tools[ind], color=colors[ind])
            plt.figure(2)
            plt.scatter(accuracy_list, usertime_list, marker='o', label=tools[ind], color=colors[ind])

        plt.figure(0)
        plt.xlabel('Accuracy [0,1]')
        plt.ylabel('F-measure [0,1]')
        plt.title('Sensitivity plot for various parameter settings of each tool (454 reads)', y=1.08)
        lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., scatterpoints=1)
        plt.savefig(output_acc_fp, bbox_extra_artists=(lgd,), bbox_inches='tight')

        plt.figure(1)
        plt.xlabel('Accuracy [0,1]')
        plt.ylabel('Wall time 64 threads (sec)')
        plt.title('Sensitivity plot for various parameter settings of each tool (454 reads)', y=1.08)
        lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., scatterpoints=1)
        plt.savefig(output_walltime_fp, bbox_extra_artists=(lgd,), bbox_inches='tight')

        plt.figure(2)
        plt.xlabel('Accuracy [0,1]')
        plt.ylabel('User time (sec)')
        plt.title('Sensitivity plot for various parameter settings of each tool (454 reads)', y=1.08)
        lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., scatterpoints=1)
        plt.savefig(output_usertime_fp, bbox_extra_artists=(lgd,), bbox_inches='tight')


@click.command()
@click.argument('accuracy_fp', required=True,
                type=click.Path(resolve_path=True, readable=True, exists=True,
                                file_okay=True))
@click.argument('output_acc_fp', required=True,
                type=click.Path(resolve_path=True, readable=True, exists=False,
                                file_okay=True))
@click.argument('output_walltime_fp', required=True,
                type=click.Path(resolve_path=True, readable=True, exists=False,
                                file_okay=True))
@click.argument('output_usertime_fp', required=True,
                type=click.Path(resolve_path=True, readable=True, exists=False,
                                file_okay=True))
def _main(accuracy_fp, output_acc_fp, output_walltime_fp, output_usertime_fp):
    """
    """

    graph_accuracy(accuracy_fp=accuracy_fp,
                   output_acc_fp=output_acc_fp,
                   output_walltime_fp=output_walltime_fp,
                   output_usertime_fp=output_usertime_fp)


if __name__ == "__main__":
    _main()
                
                
