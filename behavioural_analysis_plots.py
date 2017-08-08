from __future__ import division
import pandas as pd
import numpy as np
import os
import csv
import itertools as it
import natsort as nat
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from scipy import stats
import seaborn as sns
from behavioural_analysis_formulas import *
from behavioural_analysis_script import *

__author__ = 'Elena Maria Daniela Hindinger'

def trace_plot(dataframe, outdir, name, conc, n=48, length=25, duration=None,
               day1=None):
    df = dataframe.copy(deep=True)
    df['axis'] = (df.frame - df.frame.iat[0]) / np.timedelta64(60,
                                                               's')  # gets time series axis into right format by subtracting the date
    ptitle = 'n = ' + str(n)  # determines legend title and writes the number of animals
    sns.set_style('white')
    f, (ax) = plt.subplots(1, 1, figsize=(15, 10), dpi=600)
    ax = sns.tsplot(data=df, time='axis', value='distance', unit='animal',
                    condition='genotype', ax=ax)  # plots time series
    ax.set_xlabel('Minutes', fontsize=24)
    ax.set_ylabel('Distance (mm)', fontsize=24)
    plt.suptitle('Average total distance travelled per minute over ' + str(length) + ' hours', fontsize=36)
    lg = ax.legend(loc=1, fontsize='x-large', bbox_to_anchor=(0.97, 0.97),
                   borderaxespad=0, title=ptitle, frameon=True)
    plt.setp(lg.get_title(), fontsize='xx-large')
    setup(ax)
    if duration == 'longterm':
        x1 = 1 / length * day1
        x2 = 1 / length * (day1 + 10)
        lbo = 0
        ubo = 10
        ax.axhspan(ymin=lbo, ymax=ubo, xmin=0, xmax=x1, color='lightyellow')
        ax.axhspan(ymin=lbo, ymax=ubo, xmin=x1, xmax=x2, color='lightgrey')
        ax.axhspan(ymin=lbo, ymax=ubo, xmin=x2, xmax=1, color='lightyellow')
    else:
        pass

    ax.xaxis.set_major_locator(ticker.MultipleLocator(60))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(10))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(10))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(2))
    savename = outdir + duration + '_' + name + '_' + conc + '_in_minutes_trace.tiff'
    f.savefig(savename, format='tiff', bbox_inches='tight')
    savename = outdir + duration + '_' + name + '_' + conc + '_in_minutes_trace.pdf'
    f.savefig(savename, format='pdf', bbox_inches='tight')
    savename = outdir + duration + '_' + name + '_' + conc + '_in_minutes_trace.svg'
    f.savefig(savename, format='svg', bbox_inches='tight')
    plt.close('all')


def split_by_group(dataframe, names_list, group_number, output=None):
    df = dataframe.copy()
    gb = df.groupby('genotype')  # group by genotype
    gb.groups
    if group_number == 2:
        ''' Names list should be in format group1, group2 '''
        gr = gb.get_group(names_list[0])
        het = gb.get_group(names_list[1])
        both = pd.concat([gr, het])
        return both
    elif group_number == 4:
        ''' Names list should be in format gr, het, WIK-TL, TL '''
        gr = gb.get_group(names_list[0])
        het = gb.get_group(names_list[1])
        wiktl = gb.get_group(names_list[2])
        tl = gb.get_group(names_list[3])
        normal = pd.concat([gr, het])  # concatenate all tls
        strains = pd.concat([wiktl, tl])
        hetstrains = pd.concat([het, wiktl, tl])
        # concatenate all tlns
        if output == 'normal':
            return normal
        elif output == 'strains':
            return strains
        elif output == 'het with strains':
            return hetstrains
        else:
            print 'Sorry, this output is not supported.'
            assert False
    elif group_number == 8:
        ''' Names list should be in format gr dmso, gr con1, gr con2, gr con3,
        het dmso, het con1, het con2, het con3 '''
        gr_dmso = gb.get_group(names_list[0])
        gr_con1 = gb.get_group(names_list[1])
        gr_con2 = gb.get_group(names_list[2])
        gr_con3 = gb.get_group(names_list[3])
        het_dmso = gb.get_group(names_list[4])
        het_con1 = gb.get_group(names_list[5])
        het_con2 = gb.get_group(names_list[6])
        het_con3 = gb.get_group(names_list[7])
        dmso = pd.concat([gr_dmso, het_dmso])  # concatenate all tls
        con1 = pd.concat([gr_dmso, het_dmso, gr_con1, het_con1])
        con2 = pd.concat([gr_dmso, het_dmso, gr_con2, het_con2])
        con3 = pd.concat([gr_dmso, het_dmso, gr_con3, het_con3])  # concatenate all tlns
        if output == 'dmso':
            return dmso
        elif output == 'con1':
            return con1
        elif output == 'con2':
            return con2
        elif output == 'con3':
            return con3
        else:
            print 'Sorry, this output is not supported.'
            assert False
    else:
        print 'Sorry, this number of groups is not supported.'
        assert False


def strains(dataframe, n, length, outdir, name):
    df = dataframe.copy()
    group_mean = (df.groupby(['genotype', 'animal']).sum()) / 10 / 100
    group_mean = group_mean.reset_index()
    fig, ax = plt.subplots(1, 1, figsize=(15, 10), dpi=600)
    ax = sns.barplot(x='genotype', y='distance', data=group_mean,
                     units='animal', ax=ax)
    ax.set_xlabel('Strain', fontsize=24)
    ax.set_ylabel('Distance (m)', fontsize=24)
    plt.suptitle(
        'Strain differences in average total distance travelled over ' + str(length) + ' hours (n=' + str(n) + ')',
        fontsize=32)
    savename = outdir + name + '_quantified.tiff'
    fig.savefig(savename, format='tiff', bbox_inches='tight', dpi=600)
    savename = outdir + name + '_quantified.pdf'
    fig.savefig(savename, format='pdf', bbox_inches='tight', dpi=600)
    savename = outdir + name + '_quantified.svg'
    fig.savefig(savename, format='svg', bbox_inches='tight', dpi=600)
    plt.close('all')


def dose_response(dataframe, n, length, duration, name, outdir, conc):
    df = rearrange(dataframe.copy())
    df_m = df / 10 / 100
    fig, ax = plt.subplots(1, 1, figsize=(15, 10), dpi=600)
    ptitle = 'n = ' + str(n)  # determines legend title and writes the number of animals
    x = np.arange(4)
    ax.margins(0.1, 0)
    ax.errorbar(x, df_m['gr mean'], yerr=df_m['gr sem'], ecolor='b',
                color='b', label='gr')
    ax.errorbar(x, df_m['het mean'], yerr=df_m['het sem'], ecolor='g',
                color='g', label='het')
    plt.xticks(x, conc)
    ax.set_xlabel('Concentration (uM)', fontsize=24)
    ax.set_ylabel('Distance (m)', fontsize=24)
    plt.suptitle('Dose-response profile for average total distance travelled over ' + str(length) + ' hours',
                 fontsize=32)
    lg = ax.legend(loc=5, fontsize='x-large', bbox_to_anchor=(1.15, 0.5),
                   borderaxespad=0, title=ptitle, frameon=True)
    plt.setp(lg.get_title(), fontsize='xx-large')
    setup(ax)
    savename = outdir + name + '_dose_response_profile.tiff'
    fig.savefig(savename, format='tiff', bbox_inches='tight', dpi=600)
    savename = outdir + name + '_dose_response_profile.pdf'
    fig.savefig(savename, format='pdf', bbox_inches='tight', dpi=600)
    savename = outdir + name + '_dose_response_profile.svg'
    fig.savefig(savename, format='svg', bbox_inches='tight', dpi=600)
    plt.close('all')


def quantitative(measure, dataframe, outdir, name, n, group_names,
                 number_of_groups, day1):
    df = dataframe.copy()
    if measure == 'time':
        day_factor = (day1 + 24 - day1 - 10) * 30 * 60 * 60
        night_factor = 10 * 30 * 60 * 60
        df = (df.groupby(['genotype', 'animal', 'daytime']).sum()).reset_index()
        gb = df.groupby('daytime')
        gb.groups
        day = gb.get_group('daytime')
        night = gb.get_group('night')
        temp_day = day.loc[:, 'time'].values / day_factor * 100
        temp_night = night.loc[:, 'time'].values / night_factor * 100
        day.loc[:, 'time'] = temp_day
        night.loc[:, 'time'] = temp_night
        df = pd.concat([day, night])
        mean = df.groupby(['genotype', 'daytime']).mean().reset_index().drop(['animal', 'minute'], axis=1)
        sem = df.groupby(['genotype', 'daytime']).sem().reset_index().drop(['animal', 'minute'], axis=1)
    else:
        df = (df.groupby(['genotype', 'animal', 'daytime']).sum()) / 10 / 100
        df = df.reset_index()
        mean = df.groupby(['genotype', 'daytime']).mean().reset_index().drop(['hours', 'animal', 'minute'], axis=1)
        sem = df.groupby(['genotype', 'daytime']).sem().reset_index().drop(['hours', 'animal', 'minute'], axis=1)
    mean['sem'] = sem[measure]
    msorted = mean.sort_values(by='daytime', axis=0)
    gb = msorted.groupby('genotype')
    gb.groups
    if number_of_groups == 8:
        gr1 = gb.get_group(group_names[0])
        gr2 = gb.get_group(group_names[1])
        gr3 = gb.get_group(group_names[2])
        gr4 = gb.get_group(group_names[3])
        gr5 = gb.get_group(group_names[4])
        gr6 = gb.get_group(group_names[5])
        gr7 = gb.get_group(group_names[6])
        gr8 = gb.get_group(group_names[7])
    elif number_of_groups == 2:
        gr1 = gb.get_group(group_names[0])
        gr2 = gb.get_group(group_names[1])
    else:
        pass

    fig, ax = plt.subplots(figsize=(15, 20), dpi=600)
    setup(ax)
    timepoints = np.arange(1, 4, 2)
    width = 0.2
    labels = ['', '', 'Daytime', '', 'Nighttime', '']
    if number_of_groups == 8:
        ax.bar(timepoints - 0.7, gr1[measure], width, yerr=gr1['sem'],
               color='b', ecolor='k', capsize=10, label=group_names[0])
        ax.bar(timepoints - 0.5, gr2[measure], width, yerr=gr2['sem'],
               color='darkblue', ecolor='k', capsize=10, label=group_names[1])
        ax.bar(timepoints - 0.3, gr3[measure], width, yerr=gr3['sem'],
               color='dodgerblue', ecolor='k', capsize=10, label=group_names[2])
        ax.bar(timepoints - 0.1, gr4[measure], width, yerr=gr4['sem'],
               color='skyblue', ecolor='k', capsize=10, label=group_names[3])
        ax.bar(timepoints + 0.1, gr5[measure], width, yerr=gr5['sem'],
               color='g', ecolor='k', capsize=10, label=group_names[4])
        ax.bar(timepoints + 0.3, gr6[measure], width, yerr=gr6['sem'],
               color='darkgreen', ecolor='k', capsize=10, label=group_names[5])
        ax.bar(timepoints + 0.5, gr7[measure], width, yerr=gr7['sem'],
               color='mediumseagreen', ecolor='k', capsize=10, label=group_names[6])
        ax.bar(timepoints + 0.7, gr8[measure], width, yerr=gr8['sem'],
               color='palegreen', ecolor='k', capsize=10, label=group_names[7])
    elif number_of_groups == 2:
        ax.bar(timepoints - 0.1, gr1[measure], width, yerr=gr1['sem'],
               color='b', ecolor='k', capsize=10, label=group_names[0])
        ax.bar(timepoints + 0.1, gr2[measure], width, yerr=gr2['sem'],
               color='g', ecolor='k', capsize=10, label=group_names[4])
    else:
        pass
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(20))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(5))
    a = ax.get_xticks().tolist()
    a = labels
    ax.set_xticklabels(a)
    if measure == 'time':
        plt.suptitle('Average time in motion per condition',
                     fontsize=48)
        ax.set_ylabel('Percentage time in motion', fontsize=24)
        ax.yaxis.set_major_locator(ticker.MultipleLocator(5))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
    else:
        plt.suptitle('Average total distance moved per condition',
                     fontsize=48)
        ax.set_ylabel('Average total distance moved (m)', fontsize=24)
    ax.set_xlabel('Time of Day', fontsize=24)
    ptitle = 'n = ' + str(n)  # determines legend title and writes the number of animals
    lg = ax.legend(loc=1, fontsize='x-large', bbox_to_anchor=(0.97, 0.97),
                   borderaxespad=0, title=ptitle, frameon=True)
    plt.setp(lg.get_title(), fontsize='xx-large')
    savename = outdir + measure + '_' + name + '_quantified.tiff'
    fig.savefig(savename, format='tiff', bbox_inches='tight', dpi=600)
    savename = outdir + measure + '_' + name + '_quantified.pdf'
    fig.savefig(savename, format='pdf', bbox_inches='tight', dpi=600)
    savename = outdir + measure + '_' + name + '_quantified.svg'
    fig.savefig(savename, format='svg', bbox_inches='tight', dpi=600)
    plt.close('all')


def quantitative_3(measure, dataframe, outdir, name, n, group_names,
                   number_of_groups, day1):
    df = dataframe.copy()
    if measure == 'time':
        day1_factor = day1 * 30 * 60 * 60
        night_factor = 10 * 30 * 60 * 60
        day2_factor = (14 - day1) * 30 * 60 * 60
        df = (df.groupby(['genotype', 'animal', 'daytime']).sum()).reset_index()
        gb = df.groupby('daytime')
        gb.groups
        day1_df = gb.get_group('day1')
        night_df = gb.get_group('night')
        day2_df = gb.get_group('day2')
        temp_day1 = day1_df.loc[:, 'time'].values / day1_factor * 100
        temp_night = night_df.loc[:, 'time'].values / night_factor * 100
        temp_day2 = day2_df.loc[:, 'time'].values / day2_factor * 100
        day1_df.loc[:, 'time'] = temp_day1
        night_df.loc[:, 'time'] = temp_night
        day2_df.loc[:, 'time'] = temp_day2
        df = pd.concat([day1_df, night_df, day2_df])
        mean = df.groupby(['genotype', 'daytime']).mean().reset_index().drop(['animal', 'minute'], axis=1)
        sem = df.groupby(['genotype', 'daytime']).sem().reset_index().drop(['animal', 'minute'], axis=1)
    else:
        df = (df.groupby(['genotype', 'animal', 'daytime']).sum()) / 10 / 100
        df = df.reset_index()
        mean = df.groupby(['genotype', 'daytime']).mean().reset_index().drop(['hours', 'animal', 'minute'], axis=1)
        sem = df.groupby(['genotype', 'daytime']).sem().reset_index().drop(['hours', 'animal', 'minute'], axis=1)
    mean['sem'] = sem[measure]
    msorted = mean.sort_values(by='daytime', axis=0)
    gb = msorted.groupby('genotype')
    gb.groups
    if number_of_groups == 8:
        gr1 = gb.get_group(group_names[0]).set_index([[0, 2, 1]]).sort_index()
        gr2 = gb.get_group(group_names[1]).set_index([[0, 2, 1]]).sort_index()
        gr3 = gb.get_group(group_names[2]).set_index([[0, 2, 1]]).sort_index()
        gr4 = gb.get_group(group_names[3]).set_index([[0, 2, 1]]).sort_index()
        gr5 = gb.get_group(group_names[4]).set_index([[0, 2, 1]]).sort_index()
        gr6 = gb.get_group(group_names[5]).set_index([[0, 2, 1]]).sort_index()
        gr7 = gb.get_group(group_names[6]).set_index([[0, 2, 1]]).sort_index()
        gr8 = gb.get_group(group_names[7]).set_index([[0, 2, 1]]).sort_index()
    elif number_of_groups == 2:
        gr1 = gb.get_group(group_names[0]).set_index([[0, 2, 1]]).sort_index()
        gr2 = gb.get_group(group_names[1]).set_index([[0, 2, 1]]).sort_index()
    else:
        pass

    fig, ax = plt.subplots(figsize=(15, 20), dpi=600)
    setup(ax)
    timepoints = np.arange(1, 6, 2)
    width = 0.2
    labels = ['', '', 'Day 1', '', 'Nighttime', '', 'Day 2']
    if number_of_groups == 8:
        ax.bar(timepoints - 0.7, gr1[measure], width, yerr=gr1['sem'],
               color='b', ecolor='k', capsize=10, label=group_names[0])
        ax.bar(timepoints - 0.5, gr2[measure], width, yerr=gr2['sem'],
               color='darkblue', ecolor='k', capsize=10, label=group_names[1])
        ax.bar(timepoints - 0.3, gr3[measure], width, yerr=gr3['sem'],
               color='dodgerblue', ecolor='k', capsize=10, label=group_names[2])
        ax.bar(timepoints - 0.1, gr4[measure], width, yerr=gr4['sem'],
               color='skyblue', ecolor='k', capsize=10, label=group_names[3])
        ax.bar(timepoints + 0.1, gr5[measure], width, yerr=gr5['sem'],
               color='g', ecolor='k', capsize=10, label=group_names[4])
        ax.bar(timepoints + 0.3, gr6[measure], width, yerr=gr6['sem'],
               color='darkgreen', ecolor='k', capsize=10, label=group_names[5])
        ax.bar(timepoints + 0.5, gr7[measure], width, yerr=gr7['sem'],
               color='mediumseagreen', ecolor='k', capsize=10, label=group_names[6])
        ax.bar(timepoints + 0.7, gr8[measure], width, yerr=gr8['sem'],
               color='palegreen', ecolor='k', capsize=10, label=group_names[7])
    elif number_of_groups == 2:
        ax.bar(timepoints - 0.1, gr1[measure], width, yerr=gr1['sem'],
               color='b', ecolor='k', capsize=10, label=group_names[0])
        ax.bar(timepoints + 0.1, gr2[measure], width, yerr=gr2['sem'],
               color='g', ecolor='k', capsize=10, label=group_names[1])
    else:
        pass
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(20))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(5))
    a = ax.get_xticks().tolist()
    a = labels
    ax.set_xticklabels(a)
    if measure == 'time':
        plt.suptitle('Average time in motion per condition',
                     fontsize=48)
        ax.set_ylabel('Percentage time in motion', fontsize=32)
        ax.yaxis.set_major_locator(ticker.MultipleLocator(5))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
    else:
        plt.suptitle('Average total distance moved per condition',
                     fontsize=48)
        ax.set_ylabel('Average total distance moved (m)', fontsize=32)
    ax.set_xlabel('Time of Day', fontsize=32)
    ptitle = 'n = ' + str(n)  # determines legend title and writes the number of animals
    lg = ax.legend(loc=1, fontsize='x-large', bbox_to_anchor=(0.97, 0.97),
                   borderaxespad=0, title=ptitle, frameon=True)
    plt.setp(lg.get_title(), fontsize='xx-large')
    savename = outdir + measure + '_' + name + '_quantified.tiff'
    fig.savefig(savename, format='tiff', bbox_inches='tight', dpi=600)
    savename = outdir + measure + '_' + name + '_quantified.pdf'
    fig.savefig(savename, format='pdf', bbox_inches='tight', dpi=600)
    savename = outdir + measure + '_' + name + '_quantified.svg'
    fig.savefig(savename, format='svg', bbox_inches='tight', dpi=600)
    plt.close('all')


def binning_plot(measure, dataframe, group1, group2, n, duration, day1, outdir,
                 length, binsize):
    df = dataframe.copy(deep=True)
    mean_both = df.groupby(['genotype', 'frame']).mean().reset_index().drop('animal', axis=1)
    sem_both = df.groupby(['genotype', 'frame']).sem().reset_index().drop('animal', axis=1)
    y = length * 60 / binsize
    x_axis = np.arange(y)
    gbm = mean_both.groupby('genotype')
    gbm.groups
    gr_mean = gbm.get_group(group1).set_index(x_axis)
    het_mean = gbm.get_group(group2).set_index(x_axis)

    gbs = sem_both.groupby('genotype')
    gbs.groups
    gr_sem = gbs.get_group(group1).set_index(x_axis)
    het_sem = gbs.get_group(group2).set_index(x_axis)
    total = pd.concat([gr_mean[measure], gr_sem[measure],
                       het_mean[measure], het_sem[measure]], axis=1)
    total.columns = ['gr mean', 'gr sem', 'het mean', 'het sem']
    index = list(it.chain(np.arange((21 - day1), 24), np.arange(0, (23 - day1))))
    labels = ['12:00']
    for j in index:
        labels.append(str(j) + ':00')
    if measure == 'time':
        new = total / 108000 * 100
        shaderange = 0.5
        major_y = 1
        minor_y = 0.2
        index = list(it.chain(np.arange((22 - day1), 24),
                              np.arange(0, (23 - day1))))
        labels = ['13:00']
        for j in index:
            labels.append(str(j) + ':00')
    else:
        if binsize == 10:
            new = total / 10
            shaderange = 10
            unit = 'cm'
            major_y = 10
            minor_y = 5
        elif binsize == 5:
            new = total.copy(deep=True)
            shaderange = 40
            unit = 'mm'
            major_y = 100
            minor_y = 20
        elif binsize == 60:
            new = total / 10
            shaderange = 50
            unit = 'cm'
            major_y = 100
            minor_y = 20
            index = list(it.chain(np.arange((22 - day1), 24),
                                  np.arange(0, (23 - day1))))
            labels = ['13:00']
            for j in index:
                labels.append(str(j) + ':00')
        else:
            print 'Sorry, this binsize is not supported.'
            assert False

    fig, ax = plt.subplots(figsize=(20, 15), dpi=600)
    if measure == 'time':
        plt.suptitle('Average time in motion over 24 hours',
                     fontsize=48)
        ax.set_ylabel('Percentage of time in motion', fontsize=24)
    else:
        plt.suptitle('Average total distance moved over 24 hours',
                     fontsize=48)
        ax.set_ylabel('Distance moved (' + unit + ')', fontsize=24)

    ax.set_xlabel('Time of Day', fontsize=24)
    plt.margins(0.01, 0.01)

    ax.errorbar(x_axis, new['gr mean'], yerr=new['gr sem'], color='b',
                ecolor='lightblue', marker='s', label=group1)
    ax.errorbar(x_axis, new['het mean'], yerr=new['het sem'], color='g',
                ecolor='lightgreen', marker='s', label=group2)
    ptitle = 'n = ' + str(n)  # determines legend title and writes the number of animals
    lg = ax.legend(loc=1, fontsize='x-large', bbox_to_anchor=(0.97, 0.97),
                   borderaxespad=0, title=ptitle, frameon=True)
    plt.setp(lg.get_title(), fontsize='xx-large')
    major_freq = 60 / binsize
    minor_freq = major_freq / 2
    ax.xaxis.set_major_locator(ticker.MultipleLocator(major_freq))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(minor_freq))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(major_y))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(minor_y))
    a = ax.get_xticks().tolist()
    a = labels
    ax.set_xticklabels(a)
    setup(ax, ticksize=18)
    if duration == 'longterm':
        if binsize == 60:
            x1 = 1 / len(a) * (day1 - 1) + (1 / len(a))
            x2 = x1 + 1 / len(a) * 11 + 0.005
        else:
            x1 = 1 / len(a) * day1 + (1 / len(a)) + 0.005
            x2 = x1 + 1 / len(a) * 11
        lbo = ax.get_ylim()[0] - shaderange
        ubo = lbo + shaderange
        ax.axhspan(ymin=lbo, ymax=ubo, xmin=0, xmax=x1, color='lightyellow')
        ax.axhspan(ymin=lbo, ymax=ubo, xmin=x1, xmax=x2, color='lightgrey')
        ax.axhspan(ymin=lbo, ymax=ubo, xmin=x2, xmax=1, color='lightyellow')
    else:
        pass
    savename = os.path.join(outdir, (measure + '_' + str(binsize) + 'min_binning_' + group1 + '_' + group2 + '.tiff'))
    fig.savefig(savename, format='tiff', bbox_inches='tight', dpi=600)
    savename = os.path.join(outdir, (measure + '_' + str(binsize) + 'min_binning_' + group1 + '_' + group2 + '.pdf'))
    fig.savefig(savename, format='pdf', bbox_inches='tight', dpi=600)
    savename = os.path.join(outdir, (measure + '_' + str(binsize) + 'min_binning_' + group1 + '_' + group2 + '.svg'))
    fig.savefig(savename, format='svg', bbox_inches='tight', dpi=600)
    plt.close('all')


def plot_sliding_window(dataframe, name, group1, group2, outdir, day1, length,
                        duration, n):
    df = dataframe.copy()

    mean_both = df.groupby(['genotype', 'time']).mean().reset_index()
    sem_both = df.groupby(['genotype', 'time']).sem().reset_index()
    gbm = mean_both.groupby('genotype')
    gbm.groups
    gr_mean = gbm.get_group(group1)
    index = np.arange(len(gr_mean.index))
    gr_mean = gr_mean.set_index(index)
    het_mean = gbm.get_group(group2).set_index(index)

    gbs = sem_both.groupby('genotype')
    gbs.groups
    gr_sem = gbs.get_group(group1).set_index(index)
    het_sem = gbs.get_group(group2).set_index(index)

    total = pd.concat([gr_mean['distance'], gr_sem['distance'],
                       het_mean['distance'], het_sem['distance']], axis=1)
    total.columns = ['gr mean', 'gr sem', 'het mean', 'het sem']
    cut = -int((len(index) / 25))
    new = total.drop(total.index[cut:])
    x_axis = np.arange(len(new.index))

    initial_range = list(it.chain(np.arange((21 - day1), 24), np.arange(0, (23 - day1))))
    labels = ['12:00']
    for j in initial_range:
        labels.append(str(j) + ':00')

    fig, ax = plt.subplots(figsize=(20, 15), dpi=600)
    plt.suptitle('Progression of distance moved over 24 hours',
                 fontsize=48)
    ax.set_xlabel('Time of Day', fontsize=48)
    ax.set_ylabel('Distance moved (mm)', fontsize=48)
    plt.margins(0.01, 0.05)
    ax.plot(x_axis, new['gr mean'], color='b', linewidth=3, label=group1)
    ax.plot(x_axis, new['het mean'], color='g', linewidth=3, label=group2)
    ptitle = 'n = ' + str(n)  # determines legend title and writes the number of animals
    lg = ax.legend(loc=1, fontsize='x-large', bbox_to_anchor=(0.97, 0.97),
                   borderaxespad=0, title=ptitle, frameon=True)
    plt.setp(lg.get_title(), fontsize='xx-large')
    ax.xaxis.set_major_locator(ticker.MultipleLocator(6))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(10))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(2))
    a = ax.get_xticks().tolist()
    a = labels
    ax.set_xticklabels(a)
    setup(ax, ticksize=18)
    if duration == 'longterm':
        x1 = 1 / len(a) * day1 + (1 / len(a)) + 0.005
        x2 = x1 + 1 / len(a) * 11 + 0.005
        lbo = ax.get_ylim()[0] - 15
        ubo = ax.get_ylim()[0] - 5
        ax.axhspan(ymin=lbo, ymax=ubo, xmin=0, xmax=x1, color='lightyellow')
        ax.axhspan(ymin=lbo, ymax=ubo, xmin=x1, xmax=x2, color='lightgrey')
        ax.axhspan(ymin=lbo, ymax=ubo, xmin=x2, xmax=1, color='lightyellow')
    else:
        pass
    ymin1 = new['gr mean'] - new['gr sem']
    ymax1 = new['gr mean'] + new['gr sem']
    ymin2 = new['het mean'] - new['het sem']
    ymax2 = new['het mean'] + new['het sem']
    plt.fill_between(x_axis, ymin1, ymax1, color='lightblue')
    plt.fill_between(x_axis, ymin2, ymax2, color='lightgreen')

    savename = os.path.join(outdir, (name + '_sliding_window_' + group1 + '_' + group2 + name + '.tiff'))
    fig.savefig(savename, format='tiff', bbox_inches='tight', dpi=600)
    savename = os.path.join(outdir, (name + '_sliding_window_' + group1 + '_' + group2 + name + '.pdf'))
    fig.savefig(savename, format='pdf', bbox_inches='tight', dpi=600)
    savename = os.path.join(outdir, (name + '_sliding_window_' + group1 + '_' + group2 + name + '.svg'))
    fig.savefig(savename, format='svg', bbox_inches='tight', dpi=600)
    plt.close('all')
