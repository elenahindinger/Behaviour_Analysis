__author__ = 'Elena Maria Daniela Hindinger'

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
import behavioural_analysis_plots as bhp
from behavioural_analysis_script import *


def analysis(parent_folder, look_up_table, name_of_experiment, names_of_groups,
             concentrations, number_of_groups, n_per_group,
             length_of_experiment, type_of_experiment, day1):
    input_data_folder = os.path.join(parent_folder, 'data')
    input_data_filelist = nat.natsorted(os.listdir(input_data_folder))
    results_directory = os.path.join(parent_folder, 'results/')
    if not os.path.exists(results_directory):
        os.makedirs(results_directory)
    night = 10
    day2 = int(25 - (day1 + night))

    # builds seconds and minute arrays for later groupby
    x = (length_of_experiment * 60) + 1
    y = length_of_experiment + 1
    minutes = np.vstack([np.arange(1, x) for i in range(1800)]).T.flatten()
    hours = np.vstack([np.arange(1, y) for i in range(60)]).T.flatten()
    d1 = ['daytime'] * day1 * 60
    n = ['night'] * 10 * 60
    d2 = ['daytime'] * day2 * 60

    daytime = list(it.chain(d1, n, d2))
    # ==============================================================================
    #     all_seconds = pd.DataFrame()
    # ==============================================================================

    dfAll = pd.DataFrame()
    #  all_for_zeros = pd.DataFrame()
    counter = 1
    for i in input_data_filelist:
        print str(i) + ' is fish number ' + str(counter)
        individual_filepath = os.path.join(input_data_folder, i)
        if i.endswith('.txt'):
            my_data = bhf.readin(individual_filepath)
            df = pd.DataFrame(my_data)
            df = df.replace(['-', None], 0.0).drop(df.columns[[0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12]], axis=1)
            df.columns = ['distance']
            df['distance'] = pd.to_numeric(df['distance'], errors='coerce')
        elif i.endswith('.csv'):
            df = pd.read_csv(individual_filepath, sep=',', names=['frame', 'distance'])
        else:
            print 'Sorry, this filetype is not supported.'
            assert False
        excess_rows = (length_of_experiment * 60 * 60 * 30) - len(df.index)
        dff = df.drop(df.index[excess_rows:])
        dff['minute'] = minutes
        dff['genotype'] = bhf.lookup_dict(csv_file=look_up_table, fish=counter, output='condition')
        # ==============================================================================
        #         sec_temp = df.groupby(['genotype', 'minute', 'second']).sum()
        #         sec_temp['animal'] = a
        #         sec_temp = sec_temp.reset_index()
        #         series = pd.date_range('1/1/2017 00:00:02', periods=sec_temp.shape[0], freq='1S')
        #         sec_temp['frame'] = series
        #         all_seconds = pd.concat([all_seconds, sec_temp])
        # ==============================================================================
        dff['animal'] = counter
        # ==============================================================================
        #     all_for_zeros = pd.concat([all_for_zeros, df])
        # ==============================================================================
        min_temp = dff.groupby(['genotype', 'minute']).sum()
        min_temp['animal'] = counter
        min_temp = min_temp.reset_index()
        series_m = pd.date_range('1/1/2017 00:02:00', periods=min_temp.shape[0], freq='1T')
        min_temp['frame'] = series_m
        min_temp['hours'] = hours
        if type_of_experiment == 'longterm':
            min_temp['daytime'] = daytime
        else:
            pass
        dfAll = pd.concat([dfAll, min_temp])
        counter += 1

    minute_df = dfAll.copy(deep=True)

    if number_of_groups == 2:
        bhp.trace_plot(dataframe=minute_df, outdir=results_directory, name=name_of_experiment, conc='all', n=n_per_group,
                   length=length_of_experiment, duration=type_of_experiment, day1=day1)
        bhp.quantitative(dataframe=minute_df, outdir=results_directory, name=name_of_experiment, n=n_per_group,
                     group_names=names_of_groups)
    elif number_of_groups == 4:
        normal = bhf.split_by_group(minute_df, names_of_groups, group_number=number_of_groups, output='normal')
        strains = bhf.split_by_group(minute_df, names_of_groups, group_number=number_of_groups, output='strains')
        hetstrains = bhf.split_by_group(minute_df, names_of_groups, group_number=number_of_groups,
                                    output='het with strains')
        bhp.trace_plot(dataframe=minute_df, outdir=results_directory, name=name_of_experiment, conc='all', n=n_per_group,
                   length=length_of_experiment, duration=type_of_experiment, day1=day1)
        bhp.trace_plot(dataframe=normal, outdir=results_directory, name=name_of_experiment, conc='normal', n=n_per_group,
                   length=length_of_experiment, duration=type_of_experiment, day1=day1)
        bhp.trace_plot(dataframe=strains, outdir=results_directory, name=name_of_experiment, conc='strains', n=n_per_group,
                   length=length_of_experiment, duration=type_of_experiment, day1=day1)
        bhp.trace_plot(dataframe=hetstrains, outdir=results_directory, name=name_of_experiment, conc='alleged_wildtype',
                   n=n_per_group, length=length_of_experiment, duration=type_of_experiment, day1=day1)

        strains(dataframe=minute_df, n=n_per_group, length=length_of_experiment, outdir=results_directory,
                name=name_of_experiment)
    elif number_of_groups == 8:
        dmso = bhf.split_by_group(minute_df, names_of_groups, group_number=number_of_groups, output='dmso')
        conc1 = bhf.split_by_group(minute_df, names_of_groups, group_number=number_of_groups, output='con1')
        conc2 = bhf.split_by_group(minute_df, names_of_groups, group_number=number_of_groups, output='con2')
        conc3 = bhf.split_by_group(minute_df, names_of_groups, group_number=number_of_groups, output='con3')
        bhp.trace_plot(dataframe=minute_df, outdir=results_directory, name=name_of_experiment, conc='all', n=n_per_group,
                   length=length_of_experiment, duration=type_of_experiment, day1=day1)
        bhp.trace_plot(dataframe=dmso, outdir=results_directory, name=name_of_experiment, conc=concentrations[0],
                   n=n_per_group, length=length_of_experiment, duration=type_of_experiment, day1=day1)
        bhp.trace_plot(dataframe=conc1, outdir=results_directory, name=name_of_experiment, conc=concentrations[1],
                   n=n_per_group, length=length_of_experiment, duration=type_of_experiment, day1=day1)
        bhp.trace_plot(dataframe=conc2, outdir=results_directory, name=name_of_experiment, conc=concentrations[2],
                   n=n_per_group, length=length_of_experiment, duration=type_of_experiment, day1=day1)
        bhp.trace_plot(dataframe=conc3, outdir=results_directory, name=name_of_experiment, conc=concentrations[3],
                   n=n_per_group, length=length_of_experiment, duration=type_of_experiment, day1=day1)
        if type_of_experiment == 'acute':
            bhp.dose_response(minute_df, n_per_group, length_of_experiment, type_of_experiment, name_of_experiment,
                          results_directory, concentrations)
        elif type_of_experiment == 'longterm':
            bhp.quantitative(dataframe=minute_df, outdir=results_directory, name=name_of_experiment, n=n_per_group,
                         group_names=names_of_groups)
        else:
            pass
    else:
        print 'Sorry, number of groups is not supported.'
        assert False

    # ==============================================================================
    # savename = results_directory + 'dfAll_minutes.csv'
    # dfAll_min.to_csv(savename)
    #
    # ==============================================================================

    print 'Building finished, drawing graphs'

    # This section groups by pigment and produces 2 separate plots

    # ==============================================================================
    # tl = split_by_group(minute_df, output='tl')
    # tln = split_by_group(minute_df, output='tln')
    # trace_plot(tl, results_directory, name='tl', conc='DMSO', n=16) # plots tl
    # trace_plot(tln, results_directory, name='tln', conc='DMSO', n=16) # plots tln
    # ==============================================================================

    # ==============================================================================
    # # calculate means, stds, t-test for any dataframe df
    #
    # gb = df.groupby('animal')
    # gb.groups
    #
    # gr_sums = []
    # het_sums = []
    #
    # for group, y in gb:
    #     temp = gb.get_group(group)
    #     if temp.iloc[0, 0] == 'gr pigmented':
    #         gr_sums.append(temp['distance'][479:1079].sum())
    #     elif temp.iloc[0, 0] == 'het pigmented':
    #         het_sums.append(temp['distance'][479:1079].sum())
    #     else:
    #         pass
    #
    # gr_means = np.mean(gr_sums)
    # het_means = np.mean(het_sums)
    # gr_std = np.std(gr_sums)
    # het_std = np.std(het_sums)
    # ttest = stats.ttest_ind_from_stats(gr_means, gr_std, 48, het_means, het_std, 48)
    # ==============================================================================

    # ==============================================================================
    # # This section calculates tracking errors from DanioVision data
    # print 'Calculating tracking errors'
    # tracking = all_for_zeros.groupby('genotype')
    # tracking.groups
    #
    # gr_dmso = tracking.get_group(names_of_groups[0])
    # gr_conc1 = tracking.get_group(names_of_groups[1])
    # gr_conc2 = tracking.get_group(names_of_groups[2])
    # gr_conc3 = tracking.get_group(names_of_groups[3])
    # het_dmso = tracking.get_group(names_of_groups[4])
    # het_conc1 = tracking.get_group(names_of_groups[5])
    # het_conc2 = tracking.get_group(names_of_groups[6])
    # het_conc3 = tracking.get_group(names_of_groups[7])
    #
    # list_of_groups = [gr_dmso, gr_conc1, gr_conc2, gr_conc3, het_dmso, het_conc1,
    #                   het_conc2, het_conc3]
    #
    # print 'Tracking error per group see below'
    # tracking_error(list_of_groups, names_of_groups, results_directory)
    # ==============================================================================

    print 'All done.'
