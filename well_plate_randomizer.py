__author__ = 'Elena Maria Daniela Hindinger'

import random
import itertools as it
import pandas as pd
import numpy as np

''' ENTER PARAMETERS HERE '''

plate_size = 96
group_number = 8  # ENTER NUMBER OF GROUPS
names_of_groups = ['gr dmso', 'gr 4.6 uM', 'gr 10 uM', 'gr 20 uM',
                   'het dmso', 'het 4.6 uM', 'het 10 uM', 'het 20 uM']  # ENTER GROUP NAMES
savename = r'P:/randomization lists/25-07-17_DanioVision_matched_list_fluoxetine_longterm_repeat.csv'  # ENTER SAVE NAME AND LOCATION


''' FROM HERE ON DO NOT NEED TO MODIFY '''

all_groups = []
for group in names_of_groups:
    temp = [group] * (plate_size/group_number)
    all_groups = list(it.chain(all_groups, temp))

wells_96 = ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10', 'A11', 'A12',
            'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9', 'B10', 'B11', 'B12',
            'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10', 'C11', 'C12',
            'D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8', 'D9', 'D10', 'D11', 'D12',
            'E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E8', 'E9', 'E10', 'E11', 'E12',
            'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12',
            'G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7', 'G8', 'G9', 'G10', 'G11', 'G12',
            'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12']

wells_48 = ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'A8',
            'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8',
            'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8',
            'D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8',
            'E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E8',
            'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8']
random.shuffle(all_groups)

df = pd.DataFrame()
df['condition'] = all_groups
if plate_size == 96:
    df['wells'] = wells_96
elif plate_size == 48:
    df['wells'] = wells_48
else:
    print 'Sorry, this plate size is not supported.'
    assert False
df['fish'] = np.arange(1, plate_size+1)
df.to_csv(savename)