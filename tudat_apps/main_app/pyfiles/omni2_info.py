import numpy as np
import os
from matplotlib import pyplot as plt
import numpy.polynomial.polynomial as poly
import scipy.optimize
pi=np.pi


def fit_sin(tt, yy):
    '''Fit sin to the input time sequence, and return fitting parameters "amp", "omega", "phase", "offset", "freq", "period" and "fitfunc"'''
    tt = np.array(tt)
    yy = np.array(yy)
    ff = np.fft.fftfreq(len(tt), (tt[1]-tt[0]))   # assume uniform spacing
    Fyy = abs(np.fft.fft(yy))
    guess_freq = abs(ff[np.argmax(Fyy[1:])+1])   # excluding the zero frequency "peak", which is related to offset
    guess_amp = np.std(yy) * 2.**0.5
    guess_offset = np.mean(yy)
    guess = np.array([guess_amp, 2.*np.pi*guess_freq, 0., guess_offset])

    def sinfunc(t, A, w, p, c):  return A * np.sin(w*t + p) + c
    popt, pcov = scipy.optimize.curve_fit(sinfunc, tt, yy, p0=guess)
    A, w, p, c = popt
    f = w/(2.*np.pi)
    fitfunc = lambda t: A * np.sin(w*t + p) + c
    return {"amp": A, "omega": w, "phase": p, "offset": c, "freq": f, "period": 1./f, "fitfunc": fitfunc, "maxcov": np.max(pcov), "rawres": (guess,popt,pcov)}





# Set various project directories
pyfiles_dir = os.path.dirname(os.path.realpath(__file__))
plotfiles_dir = os.path.join(pyfiles_dir, "pyplots")
main_app_dir = os.path.abspath(os.path.join(pyfiles_dir, os.pardir))
tudat_apps_dir = os.path.abspath(os.path.join(main_app_dir, os.pardir))

# Set omni2 directory and datafile(s)
omni2_dir = os.path.abspath(os.path.join(tudat_apps_dir, "omni2_data"))
omni2YearlyFilepath = os.path.join(omni2_dir, "omni2_yearly_magfield_1964_2020.lst")
omni2MonthlyFilepath = os.path.join(omni2_dir, "omni2_monthly_magfield_1964_2020.lst")
omni2DailyFilepath = os.path.join(omni2_dir, "omni2_daily_magfield_1964_2020.lst")
omni2HourlyFilepath = os.path.join(omni2_dir, "omni2_hourly_magfield_1964_2020.lst")

def loadOmniDataFile(filePath, trimValue=None, convertTimes=True, startyear=1964, endyear=2020):

    # Load omni2 data, and remove non-data rows
    omni2DataRaw = np.genfromtxt(filePath)

    if trimValue is not None:
        omni2Data = np.delete(omni2DataRaw, np.where(omni2DataRaw==trimValue), axis=0)
    else:
        omni2Data = omni2DataRaw


    if convertTimes:
        # Convert year DOY in omni2 data to decimal year and add to new numpy array
        omni2DataFinal = np.zeros((np.size(omni2Data[:,0]), np.size(omni2Data[0,2:])))
        omni2DataFinal[:,1:] = omni2Data[:, 3:]

        for i in range(len(omni2Data)):
            year = omni2Data[i, 0]
            DOY_years = omni2Data[i, 1] / 365.25
            hour_years = (omni2Data[i, 2]/24)/365.25
            yearDecimal = year + DOY_years + hour_years

            omni2DataFinal[i, 0] = yearDecimal

    # remove unwanted data years
    omni2DataFinal = np.delete(omni2DataFinal, np.where( np.bitwise_or( (omni2DataFinal[:,0]>endyear), (omni2DataFinal[:,0]<startyear) ) )[0], 0)
    return omni2DataFinal

################ DO things to do with magfield strength directly ##################

startyear = 1964
endyear=2020

omni2DataYearly = loadOmniDataFile(omni2YearlyFilepath, trimValue=999.9, convertTimes=True, startyear=startyear, endyear=endyear)
omni2DataMonthly = loadOmniDataFile(omni2MonthlyFilepath, trimValue=999.9, convertTimes=True, startyear=startyear, endyear=endyear)
# omni2DataHourly = loadOmniDataFile(omni2HourlyFilepath, trimValue=999.9, convertTimes=True, startyear=startyear, endyear=endyear)

# Get best fit for magfield strength (using general sine, and double sine)
times = omni2DataYearly[:,0]
longtimes = np.linspace(omni2DataYearly[0,0], 2060, 1000)
lontimeOnes = np.ones(len(longtimes))
timesOnes = np.ones(len(times))
magfieldStrengths = omni2DataYearly[:,1]

magfieldStrengthsMonthly = omni2DataMonthly[:,1]
timesMonthly = omni2DataMonthly[:,0]

# magfieldStrengthsHourly = omni2DataHourly[:,1]
# timesHourly = omni2DataHourly[:,0]

def twosines(x, a1, a2, b1, b2, c1, c2, d):
    y = a1*np.sin(b1*x + c1) + a2*np.sin(b2*x + c2) + d
    return y


# bounds = ([0, 0, 0, 0, 0, 0, 0], [2, 1, 5, 2, 10, 1, 10])

b1LowerBound = (2*pi)/150
b1UpperBound = (2*pi)/80

b2LowerBound = (2*pi)/12
b2UpperBound = (2*pi)/10

bounds = ([-100, -100, b1LowerBound, b2LowerBound, -100, -100, -100], [100, 100, b1UpperBound, b2UpperBound, 100, 100, 100])
# bounds = (-100,100)
# popt, pcov = scipy.optimize.curve_fit(twosines, times, magfieldStrengths, bounds=bounds)
popt, pcov = scipy.optimize.curve_fit(twosines, timesMonthly, magfieldStrengthsMonthly, bounds=bounds)
# popt, pcov = scipy.optimize.curve_fit(twosines, timesHourly, magfieldStrengthsHourly, bounds=bounds)
a1, a2, b1, b2, c1, c2, d = popt
print("Found popt: ", popt)
print(a1, a2, b1, b2, c1, c2, d)
testYear = 2001.3
print("B0 for year %s: %s" %(testYear, twosines(testYear, a1, a2, b1, b2, c1, c2, d)))
# print(fittedData)

a1 = 1
b1 = 0.1
c1 = 4
# a2 = 1
# b2 = 10
# c2 = 0
a2 = 1.5
b2 = 0.57
c2 = -2.9
d = 6.2

#
# pguess = [a1, b1, c1, a2, b2, c2, d]
# print("pguess: ", pguess)

print("\n=========\nPopt being used manually (each parameter): ")
print("a1: %s" %a1)
print("a2: %s" %a2)
print("b1: %s" %b1)
print("b2: %s" %b2)
print("T1: %s" %(2*pi/b1))
print("T2: %s" %(2*pi/b2))
print("c1: %s" %c1)
print("c2: %s" %c2)
print("d: %s" %d)
print("==============\n")

# print("period1: ", 2*pi/b1)
# print("period2: ", 2*pi/b2)

figsize = (8, 5)

plt.figure(figsize=figsize)
# plt.title("$B_0$ over time, with sinusoidal fit function")
plt.xlabel("Year")
plt.ylabel("$B_0$ (nT)")

legend1 = []
# Plot magfield strength
plt.plot(longtimes, twosines(longtimes, a1, a2, b1, b2, c1, c2, d))
# print("LONGTIMES: %s" %longtimes)
legend1.append("Least-squares approximation")
# plt.plot(times, sinFitFunction(times))
plt.plot(omni2DataYearly[:,0], omni2DataYearly[:,1])
legend1.append("OMNI2 data (yearly average)")

# plt.plot(omni2DataMonthly[:,0], omni2DataMonthly[:,1])
# legend1.append(("OMNI2 data (monthly average)"))

plt.legend(legend1)

plt.savefig(os.path.join(plotfiles_dir, "B0_fitFunction_yearly.pdf") )
plt.savefig(os.path.join(plotfiles_dir, "B0_fitFunction_yearly.png") )


plt.figure(figsize=figsize)
# plt.title("$B_0$ over time, with sinusoidal fit function")
plt.xlabel("Year")
plt.ylabel("$B_0$ (nT)")

legend2 = []
# Plot magfield strength
plt.plot(longtimes, twosines(longtimes, a1, a2, b1, b2, c1, c2, d))
legend2.append("Sinusoidal fit-function")
# plt.plot(times, sinFitFunction(times))
# plt.plot(omni2DataYearly[:,0], omni2DataYearly[:,1])
# legend1.append("OMNI2 data (yearly average)")

plt.plot(omni2DataMonthly[:,0], omni2DataMonthly[:,1])
legend2.append(("OMNI2 data (monthly average)"))

plt.legend(legend2)

plt.savefig(os.path.join(plotfiles_dir, "B0_fitFunction_monthly.png") )



################### THings to do with magfield direction, ie phi0 #########################

# Calculate phi0 from solar wind speed (parker)

phisYearly = []
for i in range(len(omni2DataYearly[:,0])):
    vsw = omni2DataYearly[i, 8]
    phi_rad = np.arctan2(405, vsw)
    phi_deg = np.rad2deg(phi_rad)
    phisYearly.append(phi_deg)

phisMonthly = []
for i in range(len(omni2DataMonthly[:,0])):
    vsw = omni2DataMonthly[i, 8]
    phi_rad = np.arctan2(405, vsw)
    phi_deg = np.rad2deg(phi_rad)
    phisMonthly.append(phi_deg)

print(phisYearly)

meanSWYearly = np.mean(omni2DataYearly[:, 8])
meanPhiDegYearly = np.rad2deg( np.arctan2(405, meanSWYearly) )
print("Yearly mean phi ", meanPhiDegYearly)

meanSWYearly = np.mean(omni2DataYearly[:, 8])
meanPhiDegYearly = np.rad2deg( np.arctan2(405, meanSWYearly) )
print("Monthly mean phi ", meanPhiDegYearly)


b1LowerBoundPhi = -100#(2*pi)/150
b1UpperBoundPhi = 100#(2*pi)/80

b2LowerBoundPhi = (2*pi)/12
b2UpperBoundPhi = (2*pi)/10

# boundsPhi = ([-100, -100, b1LowerBoundPhi, b2LowerBoundPhi, -100, -100, -100], [100, 100, b1UpperBoundPhi, b2UpperBoundPhi, 100, 100, 100])
boundsPhi = ([-100, -100, -100000, b2LowerBoundPhi, -100, -100, -100], [100, 100, 100000, b2UpperBoundPhi, 100, 100, 100])

# popt, pcov = scipy.optimize.curve_fit(twosines, times, magfieldStrengths, bounds=bounds)
poptPhi, pcovPhi = scipy.optimize.curve_fit(twosines, timesMonthly, phisMonthly, bounds=boundsPhi)
# popt, pcov = scipy.optimize.curve_fit(twosines, timesHourly, magfieldStrengthsHourly, bounds=bounds)
a1Phi, a2Phi, b1Phi, b2Phi, c1Phi, c2Phi, dPhi = poptPhi
print(" poptPhi: ", poptPhi)

print("period1 phi: ", 2*pi/b1Phi)
print("period2 phi: ", 2*pi/b2Phi)




legend3=[]
plt.figure(figsize=figsize)
# plt.title("$\phi_0$ over time, derived from solar windspeed")
plt.xlabel("Year")
plt.ylabel("$\phi_0$ (deg)")

# plt.plot(omni2DataMonthly[:,0], omni2DataMonthly[:, 4])
# plt.plot(omni2DataYearly[:,0], omni2DataYearly[:, 8])

# plt.plot(omni2DataYearly[:,0], omni2DataYearly[:,4])
# legend3.append("Long angle of B")

plt.plot(omni2DataYearly[:,0], phisYearly)
legend3.append("$\phi_0$ (Yearly average)")

plt.plot(longtimes, twosines(longtimes, a1Phi, a2Phi, b1Phi, b2Phi, c1Phi, c2Phi, dPhi))
legend3.append("Sinusoidal fit-function")

# plt.plot(omni2DataMonthly[:,0], phisMonthly)
# legend3.append("Phi (Parker) monthly")

plt.plot(times, timesOnes*meanPhiDegYearly)
legend3.append("$\phi_0$ (Global average)")

plt.legend(legend3)

plt.savefig(os.path.join(plotfiles_dir, "phi0_fitFunction_yearly.pdf") )


legend4=[]
plt.figure(figsize=figsize)
# plt.title("$\phi_0$ over time, derived from solar windspeed")
plt.xlabel("Year")
plt.ylabel("$\phi_0$ (deg)")

# plt.plot(omni2DataMonthly[:,0], omni2DataMonthly[:, 4])
# plt.plot(omni2DataYearly[:,0], omni2DataYearly[:, 8])

# plt.plot(omni2DataYearly[:,0], omni2DataYearly[:,4])
# legend3.append("Long angle of B")

plt.plot(omni2DataMonthly[:,0], phisMonthly)
legend4.append("$\phi_0$ (Monthly average)")

plt.plot(longtimes, twosines(longtimes, a1Phi, a2Phi, b1Phi, b2Phi, c1Phi, c2Phi, dPhi))
legend4.append("Sinusoidal fit-function")

# plt.plot(omni2DataMonthly[:,0], phisMonthly)
# legend3.append("Phi (Parker) monthly")

plt.plot(times, timesOnes*meanPhiDegYearly)
legend4.append("$\phi_0$ (Global average)")

plt.legend(legend4)

plt.savefig(os.path.join(plotfiles_dir, "phi0_fitFunction_monthly.png") )


plotting = False

if plotting:
    plt.show()
