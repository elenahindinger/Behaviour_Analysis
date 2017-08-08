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
from behavioural_analysis_script import *

__author__ = 'Elena Maria Daniela Hindinger'


def readin(filepath):
    ''' This function reads input data if in txt format. '''
    datalist = []
    with open(filepath, 'rb') as csvfile:
        csvreader = csv.reader((x.replace('\0', '') for x in csvfile),
                               delimiter=';')
        for row in it.islice(csvreader, 35, None):
            datalist.append(row)
    return datalist


def setup(ax, ticksize=24):
    ''' This function determines parameters for the following graphs, such as
    determining labelsize, tick positions, padding between ticks and labels etc. '''
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(labelsize=ticksize)
    ax.tick_params(which='major', width=1.00, length=5, pad=5)
    ax.tick_params(which='minor', width=0.75, length=2.5)
    ax.tick_params(which='both', right='off')
    ax.yaxis.set_ticks_position('left')
    ax.tick_params(which='major', width=1.00, length=5, pad=5)
    ax.tick_params(which='minor', width=0.75, length=2.5)
    ax.xaxis.labelpad = 30
    ax.yaxis.labelpad = 30


def lookup_dict(csv_file, fish, output):
    look_up_table = pd.read_csv(csv_file)
    look_up_table = look_up_table.set_index('fish')
    parent_dictionary = look_up_table.to_dict(orient='dict')
    if output == 'condition':
        dictionary = parent_dictionary['condition']
        for k, v in dictionary.iteritems():
            if k == fish:
                return v
    elif output == 'well':
        dictionary = parent_dictionary['wells']
        for k, v in dictionary.iteritems():
            if k == fish:
                return v
    else:
        assert False


def sort(dataframe):
    df = dataframe.copy()
    sort = [3, 4, 2, 1, 7, 8, 6, 5]
    df['sort'] = sort
    sorted_df = df.sort_values('sort', axis=0, ascending=True)
    return sorted_df


def rearrange(dataframe):
    df = dataframe.copy()
    group_mean = df.groupby(['genotype', 'animal']).sum().reset_index().groupby('genotype').mean()
    group_sem = df.groupby(['genotype', 'animal']).sum().reset_index().groupby('genotype').sem()
    sorted_mean = sort(group_mean)
    sorted_sem = sort(group_sem)
    rearranged = pd.DataFrame()
    means = sorted_mean['distance'].tolist()
    sems = sorted_sem['distance'].tolist()
    rearranged['gr mean'] = means[:4]
    rearranged['gr sem'] = sems[:4]
    rearranged['het mean'] = means[4:]
    rearranged['het sem'] = sems[4:]
    return rearranged


def downsampler(dataframe, binsize, animal_number, condition):
    df = dataframe.copy(deep=True)
    df['frame'] = pd.date_range('1/1/2017 00:00:00', periods=df.shape[0],
                                freq='1T')
    df = df.set_index('frame').resample(binsize, label='right').sum().reset_index()
    df['animal'] = animal_number
    df['genotype'] = condition
    return df


def slider(dataframe, animal_number, condition, increments, binsize):
    dff = dataframe.copy(deep=True)
    df = pd.concat([dff['minute'], dff['distance']], axis=1)
    x = 0
    y = binsize
    result = pd.DataFrame()
    window = np.arange(0, 1500, increments)
    for x in window:
        first_mean = df.iloc[x:y, :].mean().to_frame().T
        result = pd.concat([result, first_mean])
        y += binsize
    result['time'] = np.arange(1, 1501, increments)
    result['animal'] = animal_number
    result['genotype'] = condition
    return result


def time_spent_moving(dataframe, counter, genotype, length):
    df = dataframe.copy(deep=True)
    gb = df.groupby(['genotype', 'minute'])
    gb.groups
    time = []
    for name, df in gb:
        time.append(round(((df != 0).distance.sum(axis=0)), 2))
    time_spent_moving = pd.DataFrame()
    time_spent_moving['minute'] = np.arange(1, length * 60 + 1)
    time_spent_moving['animal'] = counter
    time_spent_moving['genotype'] = genotype
    time_spent_moving['time'] = time
    return time_spent_moving
