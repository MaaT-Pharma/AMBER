#!/usr/bin/env python

# Modifications Copyright 2020 MaaT Pharma
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import argparse
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.ticker as ticker
import numpy as np
import math
import re
import os, sys, inspect
import plotly
import plotly.graph_objs as go
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)
from src.utils import load_data
from src.utils import labels
import pandas as pd 
import pickle

LEGEND2 = False


def create_colors_list():
    colors_list = []
    for color in plt.cm.tab10(np.linspace(0, 1, 10))[:-1]:
        colors_list.append(tuple(color))
    colors_list.append("black")
    for color in plt.cm.Set2(np.linspace(0, 1, 8)):
        colors_list.append(tuple(color))
    for color in plt.cm.Set3(np.linspace(0, 1, 12)):
        colors_list.append(tuple(color))
    return colors_list


def load_results(files):
    table = []
    for file in files:
        f = open(file, 'r')
        field_names = f.readline().rstrip('\n').split('\t')
        field_names[field_names.index('tool')] = 'binning_label'
        for line in f:
            values = line.rstrip('\n').split('\t')
            results_dict = dict(zip(field_names, values))
            table.append(results_dict)
        f.close()
    return table


def scan_dir(output_dir):
    p_order = re.compile('^#\([0-9]+\)')
    order = []
    data_list = []
    binning_labels = []
    for path in [d for d in (os.path.join(output_dir, d1) for d1 in os.listdir(output_dir)) if os.path.isdir(d)]:
        f = open(os.path.join(path, 'purity_completeness.tsv'), 'r')
        data_list.append(load_data.load_tsv_table(f))

        # load label and order
        f = open(os.path.join(path, 'label.txt'), 'r')
        line = f.readline().rstrip('\n')
        match = p_order.match(line)
        match_string = match.group()
        order.append(int(match_string[2:match.end() - 1]))
        binning_labels.append(line[match.end():])
        f.close()

    return data_list, binning_labels, order

def plot_heatmap_prop(df_confusion_prop, df_confusion_nb, binning_label, output_dir, separate_bar=False):
    fig, axs = plt.subplots(figsize=(10, 8))

    sns_plot = sns.heatmap(df_confusion_prop, ax=axs, annot=False, linewidths=.0, cmap="YlGnBu_r", xticklabels=True, yticklabels=True, cbar=True)

    sns_plot.set_xlabel("Expected bins", fontsize=20)
    sns_plot.set_ylabel("Predicted bins", fontsize=20)
    # plt.yticks(fontsize=8, rotation=0)
    # plt.xticks(fontsize=8)

    fig.savefig(os.path.normpath(output_dir + '/heatmap.eps'), dpi=100, format='eps', bbox_inches='tight')
    fig.savefig(os.path.normpath(output_dir + '/heatmap.png'), dpi=100, format='png', bbox_inches='tight')
    plt.close(fig)

    df_nan = df_confusion_prop.values
    df_nan[df_nan==0] = np.nan 
    hover=np.around(df_confusion_prop.values, decimals=2)
    hover= np.nan_to_num(hover).astype(np.str)
    perc_str = np.full(hover.shape, ' % of expected bin <br> Shared seqs = ', np.chararray)
    hover=np.add(np.add(hover,perc_str),np.nan_to_num(df_confusion_nb.values).astype(np.str))

    # Heatmap Plotly
    trace = go.Heatmap(
        z=df_nan, 
        y=["_" + x for x in df_confusion_prop.index],
        x=["_" + x for x in df_confusion_prop.columns], 
        # y= df_confusion_prop.index,
        # x=df_confusion_prop.columns, 
        text=hover,
        colorscale="YlOrRd",
        hoverinfo="x+y+text",
        reversescale=True
    )
    
    data=[trace]

    layout = go.Layout(
        title='',
        xaxis=dict(title='Expected bins',ticks='outside', showticklabels=True,showgrid=False,showline=True,automargin=True),
        yaxis=dict(title='Predicted bins', ticks='outside', dtick=1,showticklabels=True,autorange='reversed',showgrid=False,showline=True,automargin=True), 
        font=dict(size=20),
        hoverlabel=dict(font=dict(size=30))

    )

    plotly.offline.plot(dict(data=data, layout=layout), filename=os.path.normpath(output_dir +'/heatmap_2.html'),auto_open=False)

    fig, axs = plt.subplots(figsize=(10, 8))

    df_nan_2 = pd.DataFrame(df_nan)
    df_nan_2.index = df_confusion_prop.index
    df_nan_2.columns = df_confusion_prop.columns

    sns_plot = sns.heatmap(df_nan_2, ax=axs, annot=False, linewidths=.0, cmap="YlOrRd", xticklabels=True, yticklabels=True, cbar=True)

    sns_plot.set_xlabel("Expected bins", fontsize=10)
    sns_plot.set_ylabel("Predicted bins", fontsize=10)
    plt.yticks(fontsize=8)
    plt.xticks(fontsize=8,rotation=-45, ha='left')
    plt.title(binning_label)

    with open(os.path.normpath(output_dir + '/heatmap_2.pkl'), "wb") as f: 
        pickle.dump((df_nan_2,binning_label),f,pickle.HIGHEST_PROTOCOL)        

    fig.savefig(os.path.normpath(output_dir + '/heatmap_2.eps'), dpi=300, format='eps', bbox_inches='tight')
    fig.savefig(os.path.normpath(output_dir + '/heatmap_2.png'), dpi=300, format='png', bbox_inches='tight')
    plt.close(fig)

    if not separate_bar:
        return

    # create separate figure for bar
    fig = plt.figure(figsize=(6, 6))
    mappable = sns_plot.get_children()[0]
    fmt = lambda x, pos: '{:.0f}'.format(x / 1000000)

    cbar = plt.colorbar(mappable, orientation='vertical', label='[millions]', format=ticker.FuncFormatter(fmt))

    text = cbar.ax.yaxis.label
    font = matplotlib.font_manager.FontProperties(size=16)
    text.set_font_properties(font)

    cbar.outline.set_visible(False)
    cbar.ax.tick_params(labelsize=14)

    # store separate bar figure
    plt.gca().set_visible(False)
    fig.savefig(os.path.normpath(output_dir + '/heatmap_bar.eps'), dpi=100, format='eps', bbox_inches='tight')

    plt.close(fig)



def plot_heatmap(df_confusion, output_dir, separate_bar=False):
    fig, axs = plt.subplots(figsize=(10, 8))

    sns_plot = sns.heatmap(df_confusion, ax=axs, annot=False, cmap="YlGnBu_r", xticklabels=False, yticklabels=False, cbar=True)
    sns_plot.set_xlabel("Expected bins", fontsize=20)
    sns_plot.set_ylabel("Predicted bins", fontsize=20)
    # plt.yticks(fontsize=8, rotation=0)
    # plt.xticks(fontsize=8)

    fig.savefig(os.path.normpath(output_dir + '/heatmap.eps'), dpi=100, format='eps', bbox_inches='tight')
    fig.savefig(os.path.normpath(output_dir + '/heatmap.png'), dpi=100, format='png', bbox_inches='tight')
    plt.close(fig)

    if not separate_bar:
        return

    # create separate figure for bar
    fig = plt.figure(figsize=(6, 6))
    mappable = sns_plot.get_children()[0]
    fmt = lambda x, pos: '{:.0f}'.format(x / 1000000)

    cbar = plt.colorbar(mappable, orientation='vertical', label='[millions]', format=ticker.FuncFormatter(fmt))

    text = cbar.ax.yaxis.label
    font = matplotlib.font_manager.FontProperties(size=16)
    text.set_font_properties(font)

    cbar.outline.set_visible(False)
    cbar.ax.tick_params(labelsize=14)

    # store separate bar figure
    plt.gca().set_visible(False)
    fig.savefig(os.path.normpath(output_dir + '/heatmap_bar.eps'), dpi=100, format='eps', bbox_inches='tight')

    plt.close(fig)


def plot_boxplot(data_list, binning_labels, metric_name, onlyComplete, output_dir, order=None):
    precision_all = []
    for metrics in data_list:
        precision = []
        for metric in metrics:
            if not math.isnan(metric[metric_name]):
                if onlyComplete : 
                    if  metric[metric_name] !=0 :  
                        precision.append(metric[metric_name])
                else : 
                    precision.append(metric[metric_name])
        precision_all.append(precision)

    if order:
        # sort binning_labels and precision_all by order
        enum_order = [(v, k) for k, v in enumerate(order)]
        enum_order = sorted(enum_order, key=lambda x: x[0])
        binning_labels = [binning_labels[i[1]] for i in enum_order]
        precision_all = [precision_all[i[1]] for i in enum_order]

    df_precision_all=pd.DataFrame(precision_all)
    una_include=str(not onlyComplete)
    df_precision_all.to_csv(os.path.normpath(output_dir + '/df_boxplot_' + metric_name + '_unassigned' +  una_include + '.csv'), sep='\t')

    fig, axs = plt.subplots(figsize=(6, 5))

    medianprops = dict(linewidth=2.5, color='gold')
    bplot = axs.boxplot(precision_all, notch=0, vert=0, patch_artist=True, labels=binning_labels, medianprops=medianprops, sym='k.')
    colors_iter = iter(create_colors_list())

    # turn on grid
    axs.grid(which='major', linestyle='-', linewidth='0.5', color='lightgrey')

    # force axes to be from 0 to 100%
    axs.set_xlim([-0.01, 1.01])

    # transform plot_labels to percentages
    vals = axs.get_xticks()
    axs.set_xticklabels(['{:3.0f}'.format(x * 100) for x in vals])

    # enable code to rotate labels
    tick_labels = axs.get_yticklabels()
    plt.setp(tick_labels, fontsize=14) ## rotation=55

    for box in bplot['boxes']:
        box.set(facecolor=next(colors_iter), linewidth=0.1)
    plt.ylim(plt.ylim()[::-1])

    if metric_name == 'purity':
        axs.set_xlabel('Purity per predicted bin $p$ (%)' if LEGEND2 else 'Purity per predicted bin (%)', fontsize=14)
    elif onlyComplete :
        axs.set_xlabel('Completeness per predicted bin $r$ (%)' if LEGEND2 else 'Completeness per predicted bin (%)', fontsize=14) 
        metric_name += '_mapped'
    else : 
        axs.set_xlabel('Completeness per predicted bin and unassigned expected bin $r$ (%)' if LEGEND2 else 'Completeness per predicted bin and unassigned expected bin (%)', fontsize=14)

    fig.savefig(os.path.normpath(output_dir + '/boxplot_' + metric_name + '.pdf'), dpi=100, format='pdf', bbox_inches='tight')
    fig.savefig(os.path.normpath(output_dir + '/boxplot_' + metric_name + '.png'), dpi=100, format='png', bbox_inches='tight')
    fig.savefig(os.path.normpath(output_dir + '/boxplot_' + metric_name + '.eps'), dpi=100, format='eps', bbox_inches='tight')

    # remove labels but keep grid
    axs.get_yaxis().set_ticklabels([])
    for tic in axs.yaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False
        tic.label1On = tic.label2On = False
    fig.savefig(os.path.normpath(output_dir + '/boxplot_' + metric_name + '_wo_legend.eps'), dpi=100, format='eps', bbox_inches='tight')
    plt.close(fig)


def plot_summary(summary_per_query, output_dir, plot_type, file_name, xlabel, ylabel):
    colors_list = create_colors_list()
    if len(summary_per_query) > len(colors_list):
        raise RuntimeError("Plot only supports 29 colors")

    fig, axs = plt.subplots(figsize=(6, 5))

    # force axes to be from 0 to 100%
    axs.set_xlim([0.0, 1.0])
    axs.set_ylim([0.0, 1.0])

    i = 0
    plot_labels = []
    if plot_type == 'e':
        for summary in summary_per_query:
            axs.errorbar(float(summary[labels.AVG_PRECISION]), float(summary[labels.AVG_RECALL]), xerr=float(summary[labels.SEM_PRECISION]), yerr=float(summary[labels.SEM_RECALL]),
                         fmt='o',
                         ecolor=colors_list[i],
                         mec=colors_list[i],
                         mfc=colors_list[i],
                         capsize=3,
                         markersize=8)
            plot_labels.append(summary[labels.TOOL])
            i += 1
    if plot_type == 'w':
        for summary in summary_per_query:
            axs.plot(float(summary[labels.PRECISION]), float(summary[labels.RECALL]), marker='o', color=colors_list[i], markersize=10)
            plot_labels.append(summary[labels.TOOL])
            i += 1
    elif plot_type == 'p':
        for summary in summary_per_query:
            axs.plot(float(summary[labels.ARI_BY_BP]), float(summary[labels.PERCENTAGE_ASSIGNED]), marker='o', color=colors_list[i], markersize=10)
            plot_labels.append(summary[labels.TOOL])
            i += 1
    elif plot_type == 'q' : 
        for summary in summary_per_query:
            axs.plot(float(summary[labels.ARI_BY_SEQ]), float(summary[labels.PERCENTAGE_ASSIGNED]), marker='o', color=colors_list[i], markersize=10)
            plot_labels.append(summary[labels.TOOL])
            i += 1

    # turn on grid
    axs.minorticks_on()
    axs.grid(which='major', linestyle='-', linewidth='0.5')
    axs.grid(which='minor', linestyle=':', linewidth='0.5')

    # transform plot_labels to percentages
    vals = axs.get_xticks()
    axs.set_xticklabels(['{:3.0f}'.format(x * 100) for x in vals])
    vals = axs.get_yticks()
    axs.set_yticklabels(['{:3.0f}'.format(x * 100) for x in vals])

    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    plt.tight_layout()
    fig.savefig(os.path.normpath(output_dir + '/' + file_name + '.eps'), dpi=100, format='eps', bbox_inches='tight')

    colors_iter = iter(colors_list)
    circles = []
    for summary in summary_per_query:
        circles.append(Line2D([], [], markeredgewidth=0.0, linestyle="None", marker="o", markersize=11, markerfacecolor=next(colors_iter)))

    lgd = plt.legend(circles, plot_labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., handlelength=0, frameon=False, fontsize=12)

    fig.savefig(os.path.normpath(output_dir + '/' + file_name + '.png'), dpi=100, format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(os.path.normpath(output_dir + '/' + file_name + '.pdf'), dpi=100, format='pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig)


def plot_weighed_precision_recall_seq(summary_per_query, output_dir):
    plot_summary(summary_per_query,
                 output_dir,
                 'w',
                 'avg_purity_completeness',
                 'Average purity per sequence $\overline{p}_{seq}$ (%)' if LEGEND2 else 'Average purity per sequence (%)',
                 'Average completeness per sequence $\overline{r}_{seq}$ (%)' if LEGEND2 else 'Average completeness per sequence (%)')


def plot_adjusted_rand_index_vs_assigned_seq(summary_per_query, output_dir):
    plot_summary(summary_per_query,
                 output_dir,
                 'q',
                 'ari_vs_assigned',
                 'Adjusted Rand Index (%)' if LEGEND2 else 'Adjusted Rand Index',
                 'Percentage of assigned sequences (%)' if LEGEND2 else 'Percentage of assigned sequences')


def plot_avg_precision_recall(summary_per_query, output_dir):
    plot_summary(summary_per_query,
                 output_dir,
                 'e',
                 'avg_purity_completeness_per_bin',
                 'Truncated average purity per bin $\overline{p}_{99}$ (%)' if LEGEND2 else 'Average purity per bin (%)',
                 'Average completeness per genome $\overline{r}$ (%)' if LEGEND2 else 'Average completeness per genome (%)')


def plot_weighed_precision_recall(summary_per_query, output_dir):
    plot_summary(summary_per_query,
                 output_dir,
                 'w',
                 'avg_purity_completeness',
                 'Average purity per base pair $\overline{p}_{bp}$ (%)' if LEGEND2 else 'Average purity per base pair (%)',
                 'Average completeness per base pair $\overline{r}_{bp}$ (%)' if LEGEND2 else 'Average completeness per base pair (%)')


def plot_adjusted_rand_index_vs_assigned_bps(summary_per_query, output_dir):
    plot_summary(summary_per_query,
                 output_dir,
                 'p',
                 'ari_vs_assigned',
                 'Adjusted Rand Index (%)' if LEGEND2 else 'Adjusted Rand Index',
                 'Percentage of assigned base pairs (%)' if LEGEND2 else 'Percentage of assigned base pairs')


# def main():
#     parser = argparse.ArgumentParser(description="Create plots from one or more tables of results")
#     parser.add_argument("files", nargs='+', help="File(s) including system path")
#     parser.add_argument('-o', '--output_dir', help="Directory to save the plots in", required=True)
#     args = parser.parse_args()
#     load_data.make_sure_path_exists(args.output_dir)
#     results = load_results(args.files)
#     plot_avg_precision_recall(results, args.output_dir)
#     plot_weighed_precision_recall(results, args.output_dir)
#     plot_adjusted_rand_index_vs_assigned_bps(results, args.output_dir)


if __name__ == "__main__":
    main()
