#!/home/sapje1/miniconda2/envs/mri/bin/python3

''' 
 autoqc_status.py
 
 Check the status of existing analyses

 CJE Jan 2024
'''

import sys, os
sys.path.append('/home/sapje1/code/MriQc')

qc_types = ['fmriqc', 'spikehead', 'spikebody']
data_path = '/cubric/collab/108_QA/'
scannerlist = ['QA7T', 'QA3TM', 'QA3TW','QA3TE']

# Alternative method, which works for spikehead, spikebody too:
# check for latest directory in the proc folder

for scanner in scannerlist:
    print('-----------------------------------------')
    for qc in qc_types:
        dd = os.path.join(data_path, scanner, qc, 'proc')
        qc_dates = []
        proc_files = os.listdir(dd)
        if len(proc_files) > 0: 
            for f in proc_files:
                fsplt = f.split('-')[0:2] # take first two (date, time)
                if len(fsplt[0]) > 8:  # it's XA yyyy_mm_dd convention: reduce to yy_mm_dd
                    fsplt0new = fsplt[0][2:]
                else:
                    fsplt0new = fsplt[0]
        #        print(fsplt[0] + '   >>>    ' + fsplt0new)
                qc_dates.append(fsplt0new + '-' + fsplt[1])
            qc_dates.sort()
            latest_date = qc_dates[-1]
        else:
            latest_date = '!!no qc data'
        if qc == 'fmriqc':
            pad = '   '
        else:
            pad = ''
        print(scanner + ' last ' + qc + pad +' ' + latest_date)
        

