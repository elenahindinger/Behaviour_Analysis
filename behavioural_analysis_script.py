__author__ = 'Elena Maria Daniela Hindinger'

from __future__ import division
import os
import behavioural_analysis_function as ba



# MODULAR
parent_folder = r'I:\Elena H\DanioVision Custom Tracker Output\longterm_fluoxetine'
look_up_table = os.path.join(parent_folder, 'randomization list', '13-06-17_DanioVision_matched_list_fluoxetine_long_term.csv')
name_of_experiment = 'longterm_fluoxetine'
names_of_groups = ['gr dmso', 'gr 4.6 uM', 'gr 10 uM', 'gr 20 uM', 'het dmso',
                   'het 4.6 uM', 'het 10 uM', 'het 20 uM']
concentrations = ['dmso', '4.6 uM', '10 uM', '20 uM']
number_of_groups = 8
n_per_group = int(96 / number_of_groups)
type_of_experiment = 'longterm'  # longterm or acute
length_of_experiment = 25  # in hours
day1 = 8

ba.analysis(parent_folder, look_up_table, name_of_experiment, names_of_groups,
            concentrations, number_of_groups, n_per_group, length_of_experiment,
            type_of_experiment, day1)