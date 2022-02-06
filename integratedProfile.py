### loading from comparesingleDAy
import numpy as np
import crawlerPlots
import crawlerTools
from matplotlib import pyplot as plt
from scipy import signal
import os, sys
import datetime as DT
from testbedutils import sblib as sb
sys.path.append('/home/spike/repos')
from getdatatestbed import getDataFRF
################################
GPSfname = "/home/spike/data/20211019/20211019_174456.103_telemetry.gssbin_GPS_STAT_2.csv"  # 0.69
offset = 4.6482 + 0.0843 + 0.32 - 0.0254
lineWindow = 10  # defines how much wiggle room to identify a crawler line from above y definitions.
searchRadius = 15  # definition of comparison window between crawler and survey

########################
start = DT.datetime.strptime(os.path.dirname(GPSfname).split('/')[-1], "%Y%m%d")
if start == DT.datetime(2021, 10, 20, 0, 0):
    # time: 2021-10-20, Hs: 0.74
    fname = "/data/FRF/geomorphology/elevationTransects/survey/FRF_geomorphology_elevationTransects_survey_20211021.nc"
    go.epochd1 = (go.d1 - DT.timedelta(days=2)).timestamp()
    go.end = DT.datetime(2021, 10, 22)
    go.epochd2 = go.end.timestamp()
    lineNumbers = [ 79,  128, 172, 224, 266, 311, 357, 407, 457] # complicated: -65, -16, 25,
    lineAngles = [78, 260]
    lineWindow = 8  # defines how much wiggle room to identify a crawler line from above y definitions.
    searchRadius = 15

elif start == DT.datetime(2021, 10, 18):   #manually derived values
    # time: 2021-10-18, Hs: 1.22
    fname = "/data/FRF/geomorphology/elevationTransects/survey" \
            "/FRF_geomorphology_elevationTransects_survey_20211019.nc"
    idxSurvey = 14 # 20211019.nc'

    # go.end = start + DT.timedelta(days=3)
    # go.epochd2 = go.end.timestamp()
    lineAngles = [75, 255]
    lineNumbers = [592, 637]
    # complicated line numbers:  592, 685
    go.end = start + DT.timedelta(days=5)
    go.epochd2 = go.end.timestamp()
    lineWindow = 8  # defines how much wiggle room to identify a crawler line from above y definitions.

elif start == DT.datetime(2021, 10, 19):
    fname = "/data/FRF/geomorphology/elevationTransects/survey/FRF_geomorphology_elevationTransects_survey_20211019.nc"
    lineAngles = [71, 252]
    lineNumbers = [597, 731, 778, 870]
    # complicated line numbers: 640, 688, 824, 915, 959, 1000

elif start == DT.datetime(2021, 10, 5): # date in ['20211005']:
    # time: 2021-10-05, Hs: 0.65
    fname = '/data/FRF/geomorphology/elevationTransects/survey/FRF_geomorphology_elevationTransects_survey_20211006.nc'
    lineNumbers = [639, 453, 413]
    lineAngles = [72, 249]
    idxSurvey = 10
    go.end = start + DT.timedelta(days=3)
    go.epochd2 = go.end.timestamp()

elif start == DT.datetime(2021, 10, 25):  #date in ['20211025']:
    # time: 2021-10-25, Hs: 0.76
    fname = '/data/FRF/geomorphology/elevationTransects/survey/FRF_geomorphology_elevationTransects_survey_20211024.nc'
    lineAngles = [75, 256]
    lineNumbers = [957, 866, 823, 775, 731, 686]     # complicated line number: 917
    idxSurvey = 16
    searchRadius = 15

elif start == DT.datetime(2021, 9, 28): # in ['20210928']:
    #time: 2021-09-28, Hs: 0.77
    fname = '/home/spike/data/FRF/geomorphology/elevationTransects/survey' \
            f'/FRF_geomorphology_elevationTransects_survey_20211006.nc'
    lineNumbers = [558, 639]
    lineAngles = [81, 255]
    idxSurvey = 8     # 20210928.nc
    go.end = start + DT.timedelta(days=10)
    go.epochd2 = go.end.timestamp()

else:
    bathy, topo = None, None


data = crawlerTools.loadAndMergePriorityFiles(GPSfname, verbose=False, combineDays=True)
data = crawlerTools.cleanDF(data)
data = crawlerTools.convert2FRF(data)
data = crawlerTools.rotateTranslatePoints(data, offset)
data = crawlerTools.identifyCrawlerProfileLines(data, angleWindow=25, lineLengthThreshold=40,
        consecutivePointThresh=50, fname=os.path.join(os.path.dirname(GPSfname),  ''.join(os.path.basename(
            GPSfname).split('.')[0])+f"IdentifyProfileLines.png"), lineNumbers=lineNumbers, lineAngles=lineAngles,
                                                lineWindow=lineWindow)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
idx = data['profileNumber'] == 731
time = data['time'][idx]
pitch = np.deg2rad(data['attitude_pitch_deg'][idx])
yaw = np.deg2rad(data['attitude_pitch_deg'][idx])
speed = data['NAVSoln_speed_over_ground'][idx]
xFRF = data['xFRF'][idx]
yFRF = data['yFRF'][idx]
pressure = data['NAVSoln_depth'][idx]
elevationNAVD88 = data['elevation_NAVD88_m'][idx]
print("using crawler as trugh (not survey)")
timeStep = np.diff(time)
distOverGround, delev, dforward, dlateral = [], [], [], []
bias = 0.04411434 #
for i in range(timeStep.shape[0]):
    distOverGround.append(timeStep[i]/np.timedelta64(1, 's') * speed.iloc[i])
    delev.append(np.sin(pitch.iloc[i]-bias) * distOverGround[i])
    dforward.append(np.cos(pitch.iloc[i]) * distOverGround[i])
    dlateral.append(np.sin(yaw.iloc[i]) * distOverGround[i])

# Integrate the individual increments
elevationC = np.cumsum(delev, dtype=float)
elevationC = np.insert(elevationC, 0, 0) # pad first cell w/ zero
crossShoreC = np.cumsum(dforward, dtype=float)
crossShoreC = np.insert(crossShoreC, 0, 0) # pad first cell w/ zero
# laterl = np.cumsum(dlateral, dtype=float)

# get WL to correct depths
go = getDataFRF.getObs(time.iloc[0].to_pydatetime(), time.iloc[-1].to_pydatetime())
WL = go.getWL()
depthTruth = elevationNAVD88 - np.mean(WL['WL']) # correct elevations to depths (with WL)

# create FRF to local coordinate where start of [idx] is zero
truthX = xFRF- xFRF.min()
if elevationC[-1] - elevationC[0] > 0: # profile goes up (coming from offshore)
    truthX = truthX[::-1]

plt.figure();
for win in [50]:
    # plt.plot(sb.running_mean(crossShoreC, win), -sb.running_mean(pressure, win), label=r"$\bar{P}$".format(win))
    runningP = signal.filtfilt(np.ones(win)/win, 1,  pressure)
    plt.plot(crossShoreC, -runningP, label=r'$\bar{P}$')
# for order in [3, 5, 7]:
#     for freq in [10, 20, 50, 1000]:
#         a,b = signal.butter(order, freq, 'low', analog=True )
#         lowpassP = signal.filtfilt(b, a, pressure)
#         plt.plot(crossShoreC, -lowpassP, '.-', ms =1, label=f'lowPass freq:{freq}, order:{order}')
plt.plot(truthX, depthTruth, '.', label='Crawler GPS')
plt.plot(crossShoreC, speed, linestyle='-', linewidth=1, label='speed')
plt.plot(crossShoreC, elevationC-pressure.iloc[0], '.', label='Integrated Soln')
plt.legend()
plt.xlabel('Vehicle Local Coordinate system [m]')
plt.ylabel('Depth [m]')