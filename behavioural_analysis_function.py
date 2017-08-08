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

__author__ = 'Elena Maria Daniela Hindinger'


def analysis(parent_folder, look_up_table, name_of_experiment, names_of_groups,
             concentrations, number_of_groups, n_per_group,
             length_of_experiment, type_of_experiment, day1):
    input_data_folder = os.path.join(parent_folder, 'data')
    input_data_filelist = nat.natsorted(os.listdir(input_data_folder))
    results_directory = os.path.join(parent_folder, 'results/')
    if not os.path.exists(results_directory):
        os.makedirs(results_directory)
    # ==============================================================================
    #     os.makedirs(os.path.join(results_directory, 'individual fish'))
    # ==============================================================================
    night = 10
    day2 = length_of_experiment - (day1 + night)

    # builds seconds and minute arrays for later groupby
    x = (length_of_experiment * 60) + 1
    y = length_of_experiment + 1
    minutes = np.vstack([np.arange(1, x) for i in range(1800)]).T.flatten()
    hours = np.vstack([np.arange(1, y) for i in range(60)]).T.flatten()
    d1 = ['day1'] * int(day1 * 60)
    n = ['night'] * 10 * 60
    d2 = ['day2'] * int(day2 * 60)

    daytime = list(it.chain(d1, n, d2))
    dfAll_1int = pd.DataFrame()
    sliding_window_6010 = pd.DataFrame()  # 60 min average slid by 10 min increments
    sliding_window_305 = pd.DataFrame()  # 30 min average slid by 5 min increments
    sliding_window_3010 = pd.DataFrame()  # 30 min average slid by 10 min increments
    dfAll_5int = pd.DataFrame()
    dfAll_10int = pd.DataFrame()
    dfAll_1hint = pd.DataFrame()

    time_1hint = pd.DataFrame()
    dfAll_time = pd.DataFrame()
    #  all_for_zeros = pd.DataFrame()
    counter = 1
    for i in input_data_filelist:
        print str(i) + ' is fish number ' + str(counter)
        individual_filepath = os.path.join(input_data_folder, i)
        if i.endswith('.txt'):
            my_data = readin(individual_filepath)
            df = pd.DataFrame(my_data)
            df = df.replace(['-', None], 0.0).drop(df.columns[[0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12]], axis=1)
            df.columns = ['distance']
            df['distance'] = pd.to_numeric(df['distance'], errors='coerce')
        elif i.endswith('.csv'):
            if name_of_experiment == 'longterm_fluoxetine' or name_of_experiment == 'longterm_ketamine':
                df = pd.read_csv(individual_filepath, sep=',', names=['frame', 'distance'])
            else:
                df = pd.read_csv(individual_filepath, sep=';', names=['frame', 'distance'])
        else:
            print 'Sorry, this filetype is not supported.'
            assert False
        excess_rows = (length_of_experiment * 60 * 60 * 30) - len(df.index)
        dff = df.drop(df.index[excess_rows:])
        genotype = lookup_dict(csv_file=look_up_table, fish=counter,
                               output='condition')
        dff['minute'] = minutes
        dff['genotype'] = genotype

        ''' This section is on time spent moving '''
        time = time_spent_moving(dff, counter, genotype, length_of_experiment)
        time_temp = time.copy(deep=True)
        time_temp['daytime'] = daytime
        dfAll_time = pd.concat([dfAll_time, time_temp])
        time_int1h = downsampler(dataframe=time, binsize='1H',
                                 animal_number=counter, condition=genotype)
        total = time_int1h['time'].tolist()
        percentage = []
        for i in total:
            percentage.append(i / 108000 * 100)
        time_int1h['percentage'] = percentage
        time_1hint = pd.concat([time_1hint, time_int1h])

        ''' This section deals with distance moved '''

        min_temp = dff.groupby(['genotype', 'minute']).sum().reset_index()

        int1 = min_temp.copy(deep=True)
        int1['animal'] = counter
        int1['frame'] = pd.date_range('1/1/2017 00:02:00',
                                      periods=int1.shape[0], freq='1T')
        int1['hours'] = hours
        if type_of_experiment == 'longterm':
            int1['daytime'] = daytime
        else:
            pass
        # ==============================================================================
        #         fig = plt.figure()
        #         plt.plot(int1['distance'])
        #         plt.savefig(os.path.join(results_directory, 'individual fish', 'fish_%s' % str(counter) + '.pdf'), format='pdf')
        #         plt.close('all')
        #
        # ==============================================================================
        dfAll_1int = pd.concat([dfAll_1int, int1])

        int5 = downsampler(dataframe=min_temp, binsize='5T',
                           animal_number=counter, condition=genotype)
        dfAll_5int = pd.concat([dfAll_5int, int5])

        int10 = downsampler(dataframe=min_temp, binsize='10T',
                            animal_number=counter, condition=genotype)
        dfAll_10int = pd.concat([dfAll_10int, int10])

        int1h = downsampler(dataframe=min_temp, binsize='1H',
                            animal_number=counter, condition=genotype)
        dfAll_1hint = pd.concat([dfAll_1hint, int1h])

        slider_temp = slider(dataframe=min_temp, animal_number=counter,
                             condition=genotype, increments=10, binsize=60)
        sliding_window_6010 = pd.concat([sliding_window_6010, slider_temp])

        slider_temp_305 = slider(dataframe=min_temp, animal_number=counter,
                                 condition=genotype, increments=5, binsize=30)
        sliding_window_305 = pd.concat([sliding_window_305, slider_temp_305])

        slider_temp_3010 = slider(dataframe=min_temp, animal_number=counter,
                                  condition=genotype, increments=10, binsize=30)
        sliding_window_3010 = pd.concat([sliding_window_3010, slider_temp_3010])

        counter += 1

    minute_df = dfAll_1int.copy(deep=True)
    print 'Building finished, drawing graphs'

    if number_of_groups == 2:
        trace_plot(dataframe=minute_df, outdir=results_directory,
                   name=name_of_experiment, conc='all', n=n_per_group,
                   length=length_of_experiment, duration=type_of_experiment,
                   day1=day1)
        binning_plot(measure='distance', dataframe=dfAll_5int,
                     group1=names_of_groups[0],
                     group2=names_of_groups[1], n=n_per_group,
                     duration=type_of_experiment, day1=day1,
                     outdir=results_directory, binsize=5,
                     length=length_of_experiment)
        binning_plot(measure='distance', dataframe=dfAll_10int,
                     group1=names_of_groups[0],
                     group2=names_of_groups[1], n=n_per_group,
                     duration=type_of_experiment, day1=day1,
                     outdir=results_directory, binsize=10,
                     length=length_of_experiment)
        binning_plot(measure='distance', dataframe=dfAll_1hint,
                     group1=names_of_groups[0],
                     group2=names_of_groups[1], n=n_per_group,
                     duration=type_of_experiment, day1=day1,
                     outdir=results_directory, binsize=60,
                     length=length_of_experiment)
        binning_plot(measure='time', dataframe=time_1hint,
                     group1=names_of_groups[0],
                     group2=names_of_groups[1], n=n_per_group,
                     duration=type_of_experiment, day1=day1,
                     outdir=results_directory, binsize=60,
                     length=length_of_experiment)
        plot_sliding_window(dataframe=sliding_window_6010, name='60-10',
                            group1=names_of_groups[0],
                            group2=names_of_groups[1], n=n_per_group,
                            duration=type_of_experiment, day1=day1,
                            outdir=results_directory,
                            length=length_of_experiment)
        plot_sliding_window(dataframe=sliding_window_305, name='30-5',
                            group1=names_of_groups[0],
                            group2=names_of_groups[1], n=n_per_group,
                            duration=type_of_experiment, day1=day1,
                            outdir=results_directory,
                            length=length_of_experiment)
        plot_sliding_window(dataframe=sliding_window_3010, name='30-10',
                            group1=names_of_groups[0],
                            group2=names_of_groups[1], n=n_per_group,
                            duration=type_of_experiment, day1=day1,
                            outdir=results_directory,
                            length=length_of_experiment)
        quantitative_3(measure='distance', dataframe=minute_df,
                       outdir=results_directory,
                       name=name_of_experiment, n=n_per_group,
                       group_names=names_of_groups,
                       number_of_groups=number_of_groups, day1=day1)
        quantitative_3(measure='time', dataframe=dfAll_time,
                       outdir=results_directory,
                       name=name_of_experiment, n=n_per_group,
                       group_names=names_of_groups,
                       number_of_groups=number_of_groups, day1=day1)
    elif number_of_groups == 4:
        normal = split_by_group(minute_df, names_of_groups,
                                group_number=number_of_groups, output='normal')
        strains = split_by_group(minute_df, names_of_groups,
                                 group_number=number_of_groups,
                                 output='strains')
        hetstrains = split_by_group(minute_df, names_of_groups,
                                    group_number=number_of_groups,
                                    output='het with strains')
        trace_plot(dataframe=minute_df, outdir=results_directory,
                   name=name_of_experiment, conc='all', n=n_per_group,
                   length=length_of_experiment, duration=type_of_experiment,
                   day1=day1)
        trace_plot(dataframe=normal, outdir=results_directory,
                   name=name_of_experiment, conc='normal', n=n_per_group,
                   length=length_of_experiment, duration=type_of_experiment,
                   day1=day1)
        trace_plot(dataframe=strains, outdir=results_directory,
                   name=name_of_experiment, conc='strains', n=n_per_group,
                   length=length_of_experiment, duration=type_of_experiment,
                   day1=day1)
        trace_plot(dataframe=hetstrains, outdir=results_directory,
                   name=name_of_experiment, conc='alleged_wildtype',
                   n=n_per_group, length=length_of_experiment,
                   duration=type_of_experiment, day1=day1)

        strains(dataframe=minute_df, n=n_per_group,
                length=length_of_experiment, outdir=results_directory,
                name=name_of_experiment)
    elif number_of_groups == 8:
        dmso = split_by_group(minute_df, names_of_groups,
                              group_number=number_of_groups, output='dmso')
        conc1 = split_by_group(minute_df, names_of_groups,
                               group_number=number_of_groups, output='con1')
        conc2 = split_by_group(minute_df, names_of_groups,
                               group_number=number_of_groups, output='con2')
        conc3 = split_by_group(minute_df, names_of_groups,
                               group_number=number_of_groups, output='con3')
        trace_plot(dataframe=minute_df, outdir=results_directory,
                   name=name_of_experiment, conc='all', n=n_per_group,
                   length=length_of_experiment, duration=type_of_experiment,
                   day1=day1)
        trace_plot(dataframe=dmso, outdir=results_directory,
                   name=name_of_experiment, conc=concentrations[0],
                   n=n_per_group, length=length_of_experiment,
                   duration=type_of_experiment, day1=day1)
        trace_plot(dataframe=conc1, outdir=results_directory,
                   name=name_of_experiment, conc=concentrations[1],
                   n=n_per_group, length=length_of_experiment,
                   duration=type_of_experiment, day1=day1)
        trace_plot(dataframe=conc2, outdir=results_directory,
                   name=name_of_experiment, conc=concentrations[2],
                   n=n_per_group, length=length_of_experiment,
                   duration=type_of_experiment, day1=day1)
        trace_plot(dataframe=conc3, outdir=results_directory,
                   name=name_of_experiment, conc=concentrations[3],
                   n=n_per_group, length=length_of_experiment,
                   duration=type_of_experiment, day1=day1)
        if type_of_experiment == 'acute':
            dose_response(minute_df, n_per_group, length_of_experiment,
                          type_of_experiment, name_of_experiment,
                          results_directory, concentrations)
        elif type_of_experiment == 'longterm':
            quantitative_3(measure='distance', dataframe=minute_df,
                           outdir=results_directory,
                           name=name_of_experiment, n=n_per_group,
                           group_names=names_of_groups,
                           number_of_groups=number_of_groups, day1=day1)
            quantitative_3(measure='time', dataframe=dfAll_time,
                           outdir=results_directory,
                           name=name_of_experiment, n=n_per_group,
                           group_names=names_of_groups,
                           number_of_groups=number_of_groups, day1=day1)
            binning_plot(measure='distance', dataframe=dfAll_5int,
                         group1=names_of_groups[0],
                         group2=names_of_groups[4], n=n_per_group,
                         duration=type_of_experiment, day1=day1,
                         outdir=results_directory, binsize=5,
                         length=length_of_experiment)
            binning_plot(measure='distance', dataframe=dfAll_10int,
                         group1=names_of_groups[0],
                         group2=names_of_groups[4], n=n_per_group,
                         duration=type_of_experiment, day1=day1,
                         outdir=results_directory, binsize=10,
                         length=length_of_experiment)
            binning_plot(measure='distance', dataframe=dfAll_1hint,
                         group1=names_of_groups[0],
                         group2=names_of_groups[4], n=n_per_group,
                         duration=type_of_experiment, day1=day1,
                         outdir=results_directory, binsize=60,
                         length=length_of_experiment)
            binning_plot(measure='time', dataframe=time_1hint,
                         group1=names_of_groups[0],
                         group2=names_of_groups[4], n=n_per_group,
                         duration=type_of_experiment, day1=day1,
                         outdir=results_directory, binsize=60,
                         length=length_of_experiment)
            plot_sliding_window(dataframe=sliding_window_6010, name='60-10',
                                group1=names_of_groups[0],
                                group2=names_of_groups[4], n=n_per_group,
                                duration=type_of_experiment, day1=day1,
                                outdir=results_directory,
                                length=length_of_experiment)
            plot_sliding_window(dataframe=sliding_window_305, name='30-5',
                                group1=names_of_groups[0],
                                group2=names_of_groups[4], n=n_per_group,
                                duration=type_of_experiment, day1=day1,
                                outdir=results_directory,
                                length=length_of_experiment)
            plot_sliding_window(dataframe=sliding_window_3010, name='30-10',
                                group1=names_of_groups[0],
                                group2=names_of_groups[4], n=n_per_group,
                                duration=type_of_experiment, day1=day1,
                                outdir=results_directory,
                                length=length_of_experiment)
        else:
            pass
    else:
        print 'Sorry, number of groups is not supported.'
        assert False