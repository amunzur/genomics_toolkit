import pandas as pd
import numpy as np
import os
import re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

path_mutations="/groups/wyattgrp/users/amunzur/lu_chip/results/variant_calling/CHIP_SSCS2_curated.csv"
muts=pd.read_csv(path_mutations)
muts=muts[muts["Diagnosis"]=="Prince"]

class mutations_object_single_group:
    def __init__(self, muts_df, bar_color, group1_ntotal, group2_ntotal):
        self.muts_df = muts_df
        self.bar_color = bar_color
        self.group1_ntotal = group1_ntotal
        self.group2_ntotal = group2_ntotal
        self.group1_name = "LuPSMA"
        self.group2_name = "Cabazitaxel"

class mutations_object_two_arms:
    def __init__(self, muts_df, grouping_col_name, color_dict, group1_ntotal, group2_ntotal):
        self.muts_df = muts_df
        self.grouping_col_name = grouping_col_name
        self.color_dict = color_dict
        self.group1_ntotal = group1_ntotal
        self.group2_ntotal = group2_ntotal
        self.group1_name = "LuPSMA"
        self.group2_name = "Cabazitaxel"
