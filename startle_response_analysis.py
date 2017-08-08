__author__ = 'Elena Maria Daniela Hindinger'

from __future__ import division
import pandas as pd
import numpy as np
import os
import natsort as nat
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from scipy import stats
import seaborn as sns
import itertools as it

parent_folder = r'J:/Elena Hindinger/Startle Response Results/'


def form_dataframe(date_of_experiment, input_folder=parent_folder):
    input_files_folder = input_folder + date_of_experiment + '/data'
    randomization_list_folder = input_folder + date_of_experiment + '/randomization'
    sorted_filelist = nat.natsorted(os.listdir(input_files_folder))
    sorted_randomization_list = nat.natsorted(os.listdir(randomization_list_folder))
    dfAll = pd.DataFrame()
    noex = 1
    counter = 0
    for i in sorted_filelist:
        path = os.path.join(input_files_folder, i)
        matlab_output = pd.read_csv(path)
        meta = pd.read_csv(os.path.join(randomization_list_folder, sorted_randomization_list[counter]), skiprows=2)
        df = pd.concat([meta, matlab_output], axis=1)
        dfAll = pd.concat([dfAll, df])
        # ==============================================================================
        #         probability_distribution(df, date_of_experiment, output_folder, number_of_animals=24, noe=noex)
        # ==============================================================================
        noex += 1
        counter += 1
    return dfAll


def probability_graph(trial, output_folder):
    prob_count = trial.groupby(['condition', 'Probability']).count().fish.tolist()
    percentage_prob = []
    for i in prob_count:
        percentage_prob.append(i / 120 * 100)

    fig, ax = plt.subplots(facecolor='w', figsize=(15, 10))
    fig.patch.set_facecolor('w')
    sns.set_style('white')
    bar_width = 0.2
    index = np.arange(1, 12)
    ax.bar(index, percentage_prob[:11], bar_width, color='b', label='gr')
    ax.bar(index + bar_width, percentage_prob[11:], bar_width, color='g', label='het')
    plt.xlabel('Probability of Escape', fontsize=12)
    plt.ylabel('Percentage', fontsize=12)
    plt.title('Percentage of population per probability of escape class', fontsize=16)
    plt.xticks(index + bar_width / 2, ('0.0', '0.1', '0.2', '0.3', '0.4', '0.5',
                                       '0.6', '0.7', '0.8', '0.9', '1.0'))
    plt.legend()
    savename = os.path.join(output_folder, 'probability_distribution.pdf')
    fig.savefig(savename, bbox_inches='tight', dpi=600, facecolor='w')
    plt.close('all')


def HI_graph(trial, output_folder):
    df = pd.DataFrame()
    a = ['gr TL'] * 11
    b = ['het TL'] * 11
    c = list(it.chain(a, b))
    df['condition'] = c
    df['HI rounded'] = list(it.chain([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
                                      0.9, 1.0], [0.0, 0.1, 0.2, 0.3, 0.4, 0.5,
                                                  0.6, 0.7, 0.8, 0.9, 1.0]))
    df['fish'] = 0

    hi_rounded = []
    hi = trial['HI'].tolist()
    for i in hi:
        hi_rounded.append(round(i, 1))
    trial['HI rounded'] = hi_rounded
    hi_count = trial.groupby(['condition', 'HI rounded']).count()
    hi_count = hi_count.drop(['wells', 'Probability', 'Spont', 'total animals', 'HI'], axis=1).reset_index()
    built = pd.concat([df, hi_count])
    hi_count_full = built.groupby(['condition', 'HI rounded']).sum().fish.tolist()
    gr = hi_count_full[:11]
    het = hi_count_full[11:]
    gr_sum = np.sum(gr)
    het_sum = np.sum(het)
    gr_norm = []
    for i in gr:
        gr_norm.append(i / gr_sum * 100)
    het_norm = []
    for i in het:
        het_norm.append(i / het_sum * 100)

    fig, ax = plt.subplots()
    bar_width = 0.2
    index = np.arange(1, 12)
    ax.bar(index, gr_norm, bar_width, color='b', label='gr')
    ax.bar(index + bar_width, het_norm, bar_width, color='g', label='het')
    plt.xlabel('Habituation Index', fontsize=12)
    plt.ylabel('Percentage', fontsize=12)
    plt.title('Percentage of population per habituation index class', fontsize=16)
    plt.xticks(index + bar_width / 2, ('0.0', '0.1', '0.2', '0.3', '0.4', '0.5',
                                       '0.6', '0.7', '0.8', '0.9', '1.0'))
    plt.legend()
    savename = os.path.join(output_folder, 'HI_distribution.pdf')
    fig.savefig(savename, bbox_inches='tight', dpi=600, facecolor='w')
    plt.close('all')


def simple_quant(dataframe, n, outdir, name):
    df = dataframe.copy()
    fig, ax = plt.subplots(1, 1, figsize=(15, 10))
    ax = sns.barplot(x='condition', y='Probability', data=df, units='total animals', estimator=np.median, ax=ax)
    ax.set_xlabel('Genotype', fontsize=18)
    ax.set_ylabel('Probability of Escape', fontsize=18)
    plt.suptitle('Average probability of escape per genotype (n=' + str(n) + ')', fontsize='xx-large')
    savename = outdir + name + '_probability_quantified.pdf'
    fig.savefig(savename, format='pdf', bbox_inches='tight')
    plt.close('all')

    fig2, ax = plt.subplots(1, 1, figsize=(15, 10))
    ax = sns.barplot(x='condition', y='HI', data=df, units='total animals', ax=ax)
    ax.set_xlabel('Genotype', fontsize=18)
    ax.set_ylabel('Habituation Index', fontsize=18)
    plt.suptitle('Average habituation index per genotype', fontsize='xx-large')
    savename = outdir + name + '_HI_quantified.pdf'
    fig2.savefig(savename, format='pdf', bbox_inches='tight')
    plt.close('all')


def ttest(dataframe):
    gb = dataframe.groupby('condition')
    gb.groups
    gr = gb.get_group('gr TL')['Probability']
    het = gb.get_group('het TL')['Probability']
    return stats.ttest_ind(gr, het, equal_var=False)

# ==============================================================================
# date = '22-06-17'
# output_folder = parent_folder + date + '/results/'
# trial = sf.form_dataframe(date, output_folder)
# sf.probability_graph(trial, output_folder)
# ==============================================================================

# This section combines multiple experiments
date1 = '13-06-17'
date2 = '20-06-17'
date3 = '22-06-17'
trial1 = form_dataframe(date1)
trial2 = form_dataframe(date2)
trial3 = form_dataframe(date3)

gb = trial1.groupby('condition')
gb.groups
trial1_gr = gb.get_group('gr TL')
trial1_het = gb.get_group('het TL')

both = pd.concat([trial1_gr, trial1_het, trial2, trial3])
both['total animals'] = np.arange(1, len(both.index) + 1)
probability_graph(trial=both, output_folder=parent_folder)
HI_graph(both, parent_folder)
simple_quant(dataframe=both, n=240, outdir=parent_folder, name='median')
print 'Probability of Escape, t-test results:'
print ttest(both)

print 'Finished.'

# ==============================================================================
# probability_distribution(dfAll, output_folder, result='HI')
# ==============================================================================

# date = '20-06-17'
# parent_folder = r'J:/Elena Hindinger/Startle Response Results/'
# input_files_folder = parent_folder + date + '/data'
# randomization_list_folder = parent_folder + date + '/randomization'
# output_folder = parent_folder + date + '/results'
# sorted_filelist = nat.natsorted(os.listdir(input_files_folder))
# sorted_randomization_list = nat.natsorted(os.listdir(randomization_list_folder))
#
# dfAll = pd.DataFrame()
# counter = 1
# for i in sorted_filelist:
#     path = os.path.join(input_files_folder, i)
#     matlab_output = pd.read_csv(path)
#     if counter < 3:
#         meta = pd.read_csv(os.path.join(randomization_list_folder, sorted_randomization_list[0]), skiprows=2)
#     elif counter in range(3, 5):
#         meta = pd.read_csv(os.path.join(randomization_list_folder, sorted_randomization_list[1]), skiprows=2)
#     elif counter in range(5, 7):
#         meta = pd.read_csv(os.path.join(randomization_list_folder, sorted_randomization_list[2]), skiprows=2)
#     else:
#         assert False
#     df = pd.concat([meta, matlab_output], axis=1)
#     dfAll = pd.concat([dfAll, df])
#     counter += 1
#
# dfAll['total animals'] = np.arange(1, len(dfAll.index) + 1)
#
# gb = dfAll.groupby('condition')
# gb.groups
# grtl = gb.get_group('gr TL')
# hettl = gb.get_group('het TL')
# wttl = gb.get_group('wt TL')
# grtln = gb.get_group('gr TLN')
# hettln = gb.get_group('het TLN')
# wttln = gb.get_group('wt TLN')
#
#
# def probability_distribution(df, output_folder, result):
#     gb = df.groupby('condition')
#     gb.groups
#     grtl = gb.get_group('gr TL')
#     hettl = gb.get_group('het TL')
#     wttl = gb.get_group('wt TL')
#     grtln = gb.get_group('gr TLN')
#     hettln = gb.get_group('het TLN')
#     wttln = gb.get_group('wt TLN')
#
#     f, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, figsize=(15, 10), sharex=True, sharey=True, dpi=600)
#     sns.set_style('white')
#     if result == 'Probability':
#         sns.distplot(grtl[result], bins=10, color='darkred', ax=ax1)
#         sns.distplot(hettl[result], bins=10, color='darkgreen', ax=ax2)
#         sns.distplot(wttl[result], bins=10, color='darkblue', ax=ax3)
#         sns.distplot(grtln[result], bins=10, color='red', ax=ax4)
#         sns.distplot(hettln[result], bins=10, color='lightgreen', ax=ax5)
#         sns.distplot(wttln[result], bins=10, color='lightblue', ax=ax6)
#     else:
#         sns.distplot(grtl[result], color='darkred', ax=ax1)
#         sns.distplot(hettl[result], color='darkgreen', ax=ax2)
#         sns.distplot(wttl[result], color='darkblue', ax=ax3)
#         sns.distplot(grtln[result], color='red', ax=ax4)
#         sns.distplot(hettln[result], color='lightgreen', ax=ax5)
#         sns.distplot(wttln[result], color='lightblue', ax=ax6)
#     axes = [ax1, ax2, ax3, ax4, ax5, ax6]
#     for ax in axes:
#         if result == 'Probability':
#             ax.set_xlabel('Probability of Escape', fontsize=10)
#         elif result == 'HI':
#             ax.set_xlabel('Habituation Index', fontsize=10)
#         else:
#             assert False
#         ax.set_ylabel('Count', fontsize=10)
#     ax1.set_title('GR TL', fontsize=16)
#     ax2.set_title('Het TL', fontsize=16)
#     ax3.set_title('WT TL', fontsize=16)
#     ax4.set_title('GR TLN', fontsize=16)
#     ax5.set_title('Het TLN', fontsize=16)
#     ax6.set_title('WT TLN', fontsize=16)
#     if result == 'Probability':
#         plt.suptitle('Distribution of Escape Probabilities per Group (n=48)', fontsize='xx-large')
#         savename = os.path.join(output_folder, 'probability_distribution.pdf')
#     elif result == 'HI':
#         plt.suptitle('Habituation Index per Group (n=48)', fontsize='xx-large')
#         savename = os.path.join(output_folder, 'HI.pdf')
#     else:
#         assert False
#     f.savefig(savename, bbox_inches='tight', dpi=600)
#
#
# probability_distribution(dfAll, output_folder, result='Probability')
# probability_distribution(dfAll, output_folder, result='HI')

