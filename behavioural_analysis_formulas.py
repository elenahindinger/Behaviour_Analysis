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
        csvreader = csv.reader((x.replace('\0', '') for x in csvfile), delimiter=';')
        for row in it.islice(csvreader, 35, None):
            datalist.append(row)
    return datalist


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
        normal = pd.concat([gr, het]) # concatenate all tls
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


def tracking_error(list_of_groups, names, out_directory):
    total = []
    non_zeros = []
    # name_index = 0
    for df in list_of_groups:
        non_zeros.append(df.distance.astype(bool).sum(axis=0))
        total.append(len(df.index))
#==============================================================================
#         group_total = []
#         group_non_zeros = []
#         gb = df.groupby('animal')
#         gb.groups
#         all_animals = pd.DataFrame()
#         for animal in gb:
#             temp = gb.get_group(animal)
#             group_total.append(len(temp.index))
#             group_non_zeros.append(temp.distance.astype(bool).sum(axis=0))
#             current_animal = pd.DataFrame()
#             current_animal['total'] = total
#             current_animal['nonzeros'] = non_zeros
#             current_animal['zeros'] = current_animal.total - current_animal.nonzeros
#             current_animal['tracking error'] = (current_animal.zeros / current_animal.total) * 100
#             all_animals = pd.concat([all_animals, current_animal])
#         all_animals['condition'] = names[name_index]
#         name_index += 1
#==============================================================================
    condition = pd.DataFrame()
    condition['group'] = names
    condition['total'] = total
    condition['nonzeros'] = non_zeros
    condition['zeros'] = condition.total - condition.nonzeros
    condition['tracking error'] = (condition.zeros / condition.total) * 100
    condition.to_csv(out_directory + 'tracking_error_per_group.csv')
    # all_animals.to_csv(out_directory + 'tracking_error_per_animal.csv')
    return condition
