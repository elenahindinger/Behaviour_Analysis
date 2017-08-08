from __future__ import division
import os
import pandas as pd
import numpy as np
import natsort as nat
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from scipy import stats
import seaborn as sns
import behavioural_analysis_function as ba

__author__ = 'Elena Maria Daniela Hindinger'

# MODULAR
parent_folder = r'I:\Elena H\BEHAVIOUR\Spontaneous Locomotion Experiments\DanioVision Custom Tracker Output\longterm_H2B\DanioVision tracked'
look_up_table = os.path.join(parent_folder, 'randomization list', '28-06-17_DanioVision_matched_list_H2B_gr_het_longterm.csv')
name_of_experiment = 'longterm_H2B'
#==============================================================================
# names_of_groups = ['gr dmso', 'gr 20 uM', 'gr 50 uM', 'gr 100 uM',
#                    'het dmso', 'het 20 uM', 'het 50 uM', 'het 100 uM']
# concentrations = ['dmso', '20 uM', '50 uM', '100 uM']
#==============================================================================
names_of_groups = ['gr H2B', 'het H2B']
concentrations = ['dmso']
number_of_groups = 2
n_per_group = int(96 / number_of_groups)
type_of_experiment = 'longterm'  # longterm or acute
length_of_experiment = 24  # in hours
day1 = 7

ba.analysis(parent_folder, look_up_table, name_of_experiment, names_of_groups,
            concentrations, number_of_groups, n_per_group, length_of_experiment,
            type_of_experiment, day1)
