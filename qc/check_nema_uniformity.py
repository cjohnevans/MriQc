# -*- coding: utf-8 -*-
"""
check_nema_uniformity_t2

Created on Sat Sep 16 11:12:32 2023

Temporary testing script for adding uniformity and T2 checks into mriqc

Method:
    get niis from 3TE and 3TW folders (one experiment in each)
    ? don't use echo1.  signal is less than echo 2. partial echo??
    get snr from echo2
    get uniformity from echo2
    get T2 curve from echoes 2 .. 32

@author: sapje1
"""

import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
sys.path.append('C:\\Users\\sapje1\\code\\python_mrdatamethods\\qc')
import mriqc

#
data_path = 'C:\\Users\\sapje1\\data\\2023\\2309_nema_snr_newproto'

# snr (code already in place)
f1 = os.path.join(data_path, 'QA3TE', 'QA3TE_2_NEMA_SNR_Uniformity_T2_1_e2.nii')
f2 = os.path.join(data_path, 'QA3TE', 'QA3TE_3_NEMA_SNR_Uniformity_T2_2_e2.nii')
print(f1)
print(f2)

qc = mriqc.BasicQc()
qc.snr_nema(f1,f2)

# uniformity
qc.uniformity_nema(f1)
