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
import behavioural_analysis_formulas as bhf
from behavioural_analysis_script import *

__author__ = 'Elena Maria Daniela Hindinger'

def strains(dataframe, n, length, outdir, name):
    df = dataframe.copy()
    group_mean = (df.groupby(['genotype', 'animal']).sum()) / 10 / 100
    group_mean = group_mean.reset_index()
    fig, ax = plt.subplots(1, 1, figsize=(15, 10))
    ax = sns.barplot(x='genotype', y='distance', data=group_mean, units='animal', ax=ax)
    ax.set_xlabel('Strain', fontsize=18)
    ax.set_ylabel('Distance (m)', fontsize=18)
    plt.suptitle('Strain differences in average total distance travelled over ' + str(length) + ' hours (n=' + str(n) + ')', fontsize='xx-large')
    savename = outdir + name + '_quantified.pdf'
    fig.savefig(savename, format='pdf', bbox_inches='tight')
    plt.close('all')


def dose_response(dataframe, n, length, duration, name, outdir, conc):
    df = bhf.rearrange(dataframe.copy())
    df_m = df / 10 / 100
    fig, ax = plt.subplots(1, 1, figsize=(15, 10))
    ptitle = 'n = ' + str(n)  # determines legend title and writes the number of animals
    x = np.arange(4)
    ax.margins(0.1, 0)
    ax.errorbar(x, df_m['gr mean'], yerr=df_m['gr sem'], ecolor='b', color='b', label='gr')
    ax.errorbar(x, df_m['het mean'], yerr=df_m['het sem'], ecolor='g', color='g', label='het')
    plt.xticks(x, conc)
    ax.set_xlabel('Concentration (uM)', fontsize=18)
    ax.set_ylabel('Distance (m)', fontsize=18)
    plt.suptitle('Dose-response profile for average total distance travelled over ' + str(length) + ' hours', fontsize='xx-large')
    lg = ax.legend(loc=5, fontsize='x-large', bbox_to_anchor=(1.15, 0.5),
                    borderaxespad=0, title=ptitle, frameon=True)
    plt.setp(lg.get_title(), fontsize='xx-large')
    setup(ax)
    savename = outdir + name + '_dose_response_profile.pdf'
    fig.savefig(savename, format='pdf', bbox_inches='tight')
    plt.close('all')


def quantitative(dataframe, outdir, name, n, group_names):
    df = dataframe.copy()
    df = (df.groupby(['genotype', 'animal', 'daytime']).sum()) / 10 / 100
    df = df.reset_index()
    df.columns = ['genotype', 'animal', 'Time of Day', 'minute', 'Distance (m)', 'hours']
    sns.set_style('white')
    fig, ax = plt.subplots(1, 1, figsize=(15, 10))
    plt.suptitle('Average total distance travelled during day and night (n=' + str(n) + ')', fontsize='xx-large')
    ax.set_ylabel('Distance (m)', fontsize=18)
    ax = sns.factorplot(x='Time of Day', y='Distance (m)', hue='genotype',
                        data=df, units='animal', kind='bar', ax=ax)
    savename = outdir + name + '_quantified.pdf'
    fig.savefig(savename, format='pdf', bbox_inches='tight')
    plt.close('all')


def setup(ax):
    ''' This function determines parameters for the traceplot. '''
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(which='major', width=1.00, length=5, pad=5)
    ax.tick_params(which='minor', width=0.75, length=2.5)
    ax.tick_params(which='both', right='off')
    ax.yaxis.set_ticks_position('left')
    ax.tick_params(which='major', width=1.00, length=5, pad=5)
    ax.tick_params(which='minor', width=0.75, length=2.5)
    ax.xaxis.labelpad = 30
    ax.yaxis.labelpad = 30


def trace_plot(dataframe, outdir, name, conc, n=48, length=25, duration=None, day1=None):
    df = dataframe.copy(deep=True)
    df['axis'] = (df.frame - df.frame.iat[0]) / np.timedelta64(60, 's')  # gets time series axis into right format by subtracting the date
    ptitle = 'n = ' + str(n)  # determines legend title and writes the number of animals
    sns.set_style('white')
    f, (ax1) = plt.subplots(1, 1, figsize=(15, 10))
    ax1 = sns.tsplot(data=df, time='axis', value='distance', unit='animal',
                     condition='genotype', ax=ax1)  # plots time series
    ax1.set_xlabel('Minutes', fontsize=18)
    ax1.set_ylabel('Distance (mm)', fontsize=18)
    plt.suptitle('Average total distance travelled per minute over ' + str(length) + ' hours', fontsize='xx-large')
    lg = ax1.legend(loc=5, fontsize='x-large', bbox_to_anchor=(1.25, 0.5),
                    borderaxespad=0, title=ptitle, frameon=True)
    plt.setp(lg.get_title(), fontsize='xx-large')
    setup(ax1)

#==============================================================================
#     mi, ma = ax1.get_ylim()
#     lbo = ma - 10
#     ubo = ma - 5
#==============================================================================
    if duration == 'longterm':
        x1 = 1 / length * day1
        x2 = 1 / length * (day1 + 10)
        lbo = 0
        ubo = 10
        ax1.axhspan(ymin=lbo, ymax=ubo, xmin=0, xmax=x1, color='lightyellow')
        ax1.axhspan(ymin=lbo, ymax=ubo, xmin=x1, xmax=x2, color='lightgrey')
        ax1.axhspan(ymin=lbo, ymax=ubo, xmin=x2, xmax=1, color='lightyellow')
    else:
        pass

    ax1.xaxis.set_major_locator(ticker.MultipleLocator(60))
    ax1.xaxis.set_minor_locator(ticker.MultipleLocator(10))
    ax1.yaxis.set_major_locator(ticker.MultipleLocator(10))
    ax1.yaxis.set_minor_locator(ticker.MultipleLocator(2))
    savename = outdir + duration + '_' + name + '_' + conc + '_in_minutes_trace.pdf'
    f.savefig(savename, format='pdf', bbox_inches='tight')
    plt.close('all')
