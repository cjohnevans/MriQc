import datetime as dt
import time
import os 

# needs to be running on the Aglient PC to transfer the magnet data file across to the transfer 
# directory.  Parsing and plots happen on the linux side.
# Runs via the batch file magmon7T_aglient_run.bat
# CJE April 2023

# read data from magmon file
magmon_infile = 'C:\\transfer\\022771data.txt'
output_path = 'Y:\\magnet_logs\\'

def read_write_magmon(now):
    f_in = open(magmon_infile, 'r')
    magmon = f_in.readlines()
    f_in.close()
    
    now_time = str(now.date()) + '_' + str(now.time())
    now_time = now_time.replace(':','')[0:17]
    print('Latest readings at ' + now_time + ':')
    print(magmon[-1])

    magmon_outfile = 'magmon7T_' + now_time + '.txt'
    magmon_outpath = os.path.join(output_path, magmon_outfile)
    f_out = open(magmon_outpath, 'w')
    f_out.write(magmon[0])
    for line in magmon[-101:]:
        f_out.write(line)
    f_out.close()
    print('Written to ' + magmon_outpath)
    
    # sleep for 60s here, so that data will only be written once while minute=10 is true
    time.sleep(60)

# loop and check the time here.  If 10 past, run the read_write function
while True:
    now = dt.datetime.now()
    print((now.hour, now.minute))
    # only check at 10 mins past the hour
    if now.minute == 10:
        print('running')
        read_write_magmon(now)
    # sleep for 20s.  Not longer to make sure we don't miss 10-past.
    time.sleep(20)
