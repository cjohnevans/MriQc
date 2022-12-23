import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt
import os

# get last log file from directory
logdir = '/cubric/scanners/mri/7t/magnet_logs/'
logdirlist = os.listdir(logdir)
logdirlist.sort()
logname = logdirlist[-1]
print('Latest log: ' + logdirlist[-1])
logfile = os.path.join(logdir,logname)

def strtime_to_datetime(t):
    '''
    take a list of strings in format [dd, mm, yy, hh, mm, ss] and convert to datetime object
    '''
    t_dt = dt.datetime(int(t[2]),int(t[1]),int(t[0]),int(t[3]),int(t[4]),int(t[5]))
    return t_dt

log7t = pd.read_csv(logfile, delimiter='\t')
log7t = log7t[log7t['Pressure (mB)']!=0]

meas_time = log7t['Date & time'].str.split('/|:| ', regex=True, expand=False)
# convert to datetime object
log7t['datetime'] = meas_time.apply(strtime_to_datetime)
cryo_T = log7t[['Service Shield Link Top Temp 4ABb (K)',\
                      'Patient Shield Link Top Temp 3ABb (K)',\
                      'Service Endplate Temp 2AB (K)',\
                      'Patient Endplate Temp 1AB (K)',\
                      'Service Outer Tube Bottom Temp 2ABb (K)',\
                      'Service Coldhead 1st Stage Temp 6AB (K)',\
                      'Patient Coldhead 1st Stage Temp 5AB (K)' \
                     ]]

plt.plot(log7t['datetime'].values, log7t['Pressure (mB)'].values)
plt.xticks(rotation=20)
plt.title('7T He can pressure (mB)')
plt.savefig(os.path.join(logdir,'cryostat_pressure.png'))
#plt.show()
plt.close()

plt.plot(log7t['datetime'],log7t['Service Shield Link Top Temp 4ABb (K)'], \
        log7t['datetime'],log7t['Service Endplate Temp 2AB (K)'], \
        log7t['datetime'],log7t['Patient Endplate Temp 1AB (K)'], \
        log7t['datetime'],log7t['Service Outer Tube Bottom Temp 2ABb (K)'], \
        log7t['datetime'],log7t['Patient Coldhead 1st Stage Temp 5AB (K)'])
plt.xticks(rotation=20)
#plt.show()
plt.title('Cryostat Temperatures (K)')
plt.savefig(os.path.join(logdir,'cryostat_T.png'))

print(log7t[['Date & time', 'Pressure (mB)','Service Shield Link Top Temp 4ABb (K)' ]].tail(6))
