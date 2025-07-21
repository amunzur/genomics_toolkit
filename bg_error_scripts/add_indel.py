#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 15:17:05 2024

@author: jbacon
"""

import pandas as pd

error=pd.read_csv("/groups/wyattgrp/users/jbacon/err_rates/urine_panel/error_rate/error_rates.tsv", sep='\t')

error.rename(columns={'chrom': 'CHROM'}, inplace=True)
error.rename(columns={'pos': 'POS'}, inplace=True)
error.rename(columns={'ref': 'REF'}, inplace=True)

error.to_csv("/groups/wyattgrp/users/jbacon/err_rates/urine_panel/error_rate/error_rate.fixed.bed", index=False, sep='\t')