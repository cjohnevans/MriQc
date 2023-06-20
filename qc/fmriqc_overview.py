# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 07:51:37 2023

@author: sapje1
"""
data_path = 'C:\\Users\\sapje1\\data\\mriqc_7t'
data_path = 'C:\\Users\\sapje1\\data\\mriqc_3te'
code_path = 'C:\\Users\\sapje1\\code\\python_mrdatamethods\\qc'

import os
import sys
import pandas as pd
sys.path.append(code_path)
import mriqc

qc_over = mriqc.FmriQcOverview(data_path)
#qc_over.oview_qc.info()
print(qc_over.oview_qc['scanner'].value_counts())

