
#################### README PLEASE ##########################################################################################
# Use this file if I want to do some fancy stuff from big omni2 files - is taking too long atm so generating individually
#################### README PLEASE ##########################################################################################



import numpy as np
import os
from matplotlib import pyplot as plt
import numpy.polynomial.polynomial as poly
import scipy.optimize



# Set various project directories
pyfiles_dir = os.path.dirname(os.path.realpath(__file__))
main_app_dir = os.path.abspath(os.path.join(pyfiles_dir, os.pardir))
tudat_apps_dir = os.path.abspath(os.path.join(main_app_dir, os.pardir))



raw2all = True
##################### Open full raw omni2 file, trim out nan values, and change year day hour to decimal year #########



# Set omni2 directory and raw datafile
omni2_dir = os.path.abspath(os.path.join(tudat_apps_dir, "omni2_data"))
omni2_datafile = os.path.join(omni2_dir, "omni2_all_data_raw.lst")


if raw2all:

    #################### Load omni2 data, and change time columns to decimal year ################################
    print("===== Loading omni2 data ======")
    # Load data
    omni2DataRaw = np.genfromtxt(omni2_datafile)


    print("===== Converting to decimal year ======")
    # Create model Newtimes data array
    omni2DataNewtimes = np.zeros((np.size(omni2DataRaw[:,0]), np.size(omni2DataRaw[0,2:])))
    omni2DataNewtimes[:,1:] = omni2DataRaw[:, 3:] # Columns after time the same as usual

    # Change each time column to be a decimal year
    for i in range(len(omni2DataRaw)):
        year = omni2DataRaw[i, 0]
        DOY_years = omni2DataRaw[i, 1] / 365.25
        hour_years = (omni2DataRaw[i, 2]/24)/365.25
        yearDecimal = year + DOY_years + hour_years

        omni2DataNewtimes[i, 0] = yearDecimal

    print("===== Saving newtimes omni2 data ======")

    # Save newTimes array
    np.savetxt(omni2_dir + "/omni2_all_data.txt", omni2DataNewtimes, delimiter=" ", fmt="%1.3f")


    print("===== Creating .fmt for newtimes omni2 data ======")

    NominalAllItems = ["YEAR",
                "DOY",
                "Hour",
                "Scalar B, nT",
                "Vector B Magnitude,nT",
                "Lat. Angle of B (GSE)",
                "Long. Angle of B (GSE)",
                "BX, nT (GSE, GSM)",
                "BY, nT (GSE)",
                "BZ, nT (GSE)",
                "BY, nT (GSM)",
                "BZ, nT (GSM)",
                "RMS_magnitude, nT",
                "RMS_field_vector, nT",
                "RMS_BX_GSE, nT",
                "RMS_BY_GSE, nT",
                "RMS_BZ_GSE, nT",
                "SW Plasma Temperature, K",
                "SW Proton Density, N/cm^3",
                "SW Plasma Speed, km/s",
                "SW Plasma flow long. angle",
                "SW Plasma flow lat. angle",
                "sigma-V, km/s",
                "Flow pressure",
                "E elecrtric field",
                "R (Sunspot No.)",
                "Dst-index, nT"]

    CustomItems = ["Decimal Year"]

    newtimesFmtList = [CustomItems[0]] + NominalAllItems[3:]
    newtimesFmtListIndices = np.array((range(len(newtimesFmtList))))

    newtimesFmtArray = np.empty( (len(newtimesFmtList), 2), dtype=object )
    newtimesFmtArray[:,0] = newtimesFmtListIndices
    newtimesFmtArray[:,1] = newtimesFmtList

    np.savetxt(omni2_dir + "/omni2_all_data.fmt", newtimesFmtArray, delimiter=" ", fmt="%s")





print("===== Creating magfield data files ======")

dataFileToUse = omni2_dir + "/omni2_all_data.txt"
indicesToPull = [0, 1, 2]
trimming = True
fileToSave = "omni2_mag_data.txt"
trimValue = 999.9


omni2Data = np.genfromtxt(dataFileToUse)
omni2DataSubset = omni2Data[:, indicesToPull]


if trimming:
    np.delete(omni2DataSubset, np.where(omni2DataSubset==trimValue), axis=0)

np.savetxt(omni2_dir + "/" + fileToSave, omni2DataSubset, delimiter=" ", fmt="%1.1f")
