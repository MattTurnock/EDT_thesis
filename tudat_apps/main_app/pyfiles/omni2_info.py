import numpy as np
import os
import datetime
from astropy import time


# Set various project directories
pyfiles_dir = os.path.dirname(os.path.realpath(__file__))
main_app_dir = os.path.abspath(os.path.join(pyfiles_dir, os.pardir))
tudat_apps_dir = os.path.abspath(os.path.join(main_app_dir, os.pardir))

# Set omni2 directory and datafile(s)
omni2_dir = os.path.abspath(os.path.join(tudat_apps_dir, "omni2_data"))
omni2_magfield_datafile = os.path.join(omni2_dir, "omni2_monthly_magfield_1964_2020.lst")

# Load omni2 data, and remove non-data rows
omni2DataRaw = np.genfromtxt(omni2_magfield_datafile)
omni2DataTrimmed = np.delete(omni2DataRaw, np.where(omni2DataRaw==999.9), axis=0)

# Convert year DOY in omni2 data to datetime objects
omni2DataDatetimes = np.zeros(np.size(omni2DataTrimmed))
for i in range(len(omni2DataTrimmed)):
    year = omni2DataTrimmed[i, 0]
    DOY = omni2DataTrimmed[i, 1]
    hour = omni2DataTrimmed[i, 2]
    timestring = "%s:%s:%s:%s:%s" (year, DOY, hour, 0, 0)

    print(time.Time(timestring, format="yday", scale="utc"))



omni2DataDatetimes[:0] = [datetime.datetime]

print(omni2DataTrimmed)

