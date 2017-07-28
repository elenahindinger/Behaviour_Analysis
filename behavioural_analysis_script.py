from __future__ import division
import os
import behavioural_analysis_function as ba

__author__ = 'Elena Maria Daniela Hindinger'


# MODULAR
parent_folder = r'E:\24h trial analysis\pilot experiments\trial1 pigmented\NEW ANALYSIS'
look_up_table = os.path.join(parent_folder, 'randomization list', '17-04-17_DanioVision_matched_list_gr_het_longterm.csv')
name_of_experiment = 'longterm_pilot'
names_of_groups = ['gr', 'het']
concentrations = ['dmso']
number_of_groups = 2
n_per_group = int(96 / number_of_groups)
type_of_experiment = 'longterm'  # longterm or acute
length_of_experiment = 24  # in hours
day1 = 4

ba.analysis(parent_folder, look_up_table, name_of_experiment, names_of_groups,
            concentrations, number_of_groups, n_per_group, length_of_experiment,
            type_of_experiment, day1)