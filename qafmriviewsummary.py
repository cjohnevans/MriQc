# qafmriviewsummary.py
# purpose: visualise summary data from RoutineQA

import csv
import matplotlib.pyplot as plt
import numpy as np
from datetime import date

filepath = '/home/sapje1/data_sapje1/QC/RoutineQC/summary/7TQC_epi.txt'

class QA:
    measdate = []
    sfnr = []
    fileid = []
    scanner = ''
    sfnr = []
    snr = []
    nmeas = []

    def loadQA(self, filepath):
        with open( filepath) as csvfile:
            qaread = csv.DictReader(csvfile)
            for row in qaread:
                self.scanner = row['Scanner']
                self.sfnr.append(float(row['SFNR_roi_nodrift']))
                self.snr.append(float(row['SNR']))
                print(("20" + row['Date'][0:2]) , \
                                     (row['Date'][3:5]) , \
                                     (row['Date'][-2:]) )
                self.measdate.append( date(int("20" + row['Date'][:2]), \
                                     int(row['Date'][3:5]), \
                                     int(row['Date'][-2:]) ) )
        self.nmeas = len(self.sfnr)

plt.close('all')
qa7t = QA()
#qa7t.loadQA('/home/sapje1/data_sapje1/QC/RoutineQC/summary/7TQC_epi.txt')
qa7t.loadQA('./7TQC_epi.txt')
print (qa7t.measdate)

fig1 = plt.figure(1)
ax1 = fig1.add_subplot(1,1,1)
ax1.plot(qa7t.measdate, qa7t.sfnr, 'x')
plt.xlabel('Date')
plt.ylabel('SFNR')
plt.title(qa7t.scanner)

fig2 = plt.figure(2)
ax2 = fig2.add_subplot(1,1,1)
ax2.plot(qa7t.measdate, qa7t.snr, 'o')
plt.xlabel('Date')
plt.ylabel('SNR')
plt.title(qa7t.scanner)

plt.show()



