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
idx = data['profileNumber'] == 778


# get WL to correct depths
# go = getDataFRF.getObs(time.iloc[0].to_pydatetime(), time.iloc[-1].to_pydatetime())
# WL = go.getWL()
def developProfFromPressureAndIMUprofile(data, WL, plotting=False, **kwargs):
    """Function develops profile from pressure signal and IMU pitch values from already subseted python dataframe.
    Uses a running average 2-way (filt-filt)
    Args:
        data: already subset data frame
        WL: waterlevel to adjust the elevations from the crawler
        plotting: True/false or filename.  (Default=False)
    Keyword Args:
        "integratedProfile": use IMU to derive a profile (try to identify bias and remove that as well)
    Returns:

    """
    runAvgWin = [10] #, 40, 60]
    integratedProfile = kwargs.get('integratedProfile', False)
    # first unpack data
    time = data['time']
    pitch = np.deg2rad(data['attitude_pitch_deg'])
    yaw = np.deg2rad(data['attitude_pitch_deg'])
    speed = data['NAVSoln_speed_over_ground']
    xFRF = data['xFRF']
    # yFRF = data['yFRF']
    pressure = data['NAVSoln_depth']
    elevationNAVD88 = data['elevation_NAVD88_m']
    # print("using crawler as truth (not survey)")
    timeStep = np.diff(time)
    biasWindow = [-0.25, 0.25] # bias window
    #next identify the bias for this section
    # counts, bins, _ = plt.hist(data['attitude_pitch_deg'], bins=50)
    hist, bin_edges = np.histogram(data['attitude_pitch_deg'], bins=500, density=False)
    binCenters = (bin_edges + np.diff(bin_edges).mean())[:-1]
    mask = (binCenters < biasWindow[1]) & (binCenters > biasWindow[0])
    if mask.any():
        bias = binCenters[mask][hist[mask].argmax()]
        print(f'      Bias Removed {bias:.3f}')
    else:
        bias=0
    distOverGround, delev, dforward, dlateral = [], [], [], []
    for i in range(timeStep.shape[0]):
        distOverGround.append(timeStep[i]/np.timedelta64(1, 's') * speed.iloc[i])
        delev.append(np.sin(pitch.iloc[i]-bias) * distOverGround[i])
        dforward.append(np.cos(pitch.iloc[i]) * distOverGround[i])
        dlateral.append(np.sin(yaw.iloc[i]) * distOverGround[i])
    
    # Integrate the individual increments
    elevationC = np.cumsum(delev, dtype=float)  # elevation
    elevationC = np.insert(elevationC, 0, 0) # pad first cell w/ zero
    crossShoreC = np.cumsum(dforward, dtype=float)  # cross-shore coordinate
    crossShoreC = np.insert(crossShoreC, 0, 0) # pad first cell w/ zero
    
    depthTruth = elevationNAVD88 - np.mean(WL['WL']) # correct elevations to depths (with WL)
    
    # move all of the data to FRF and NAVD88 coordinate systems #  to local coordinate where start of [idx] is zero
    truthX = xFRF
    # trueShoreLine = xFRF[np.abs(elevationNAVD88).argmin()]
    if elevationC[-1] - elevationC[0] > 0: # profile goes up (coming from offshore)
        runningP, shorelineP, crossShoreP= {}, {}, {}
        # crossShoreC = crossShoreC + trueShoreLine
        crossShoreP = xFRF.max() - crossShoreC.max() + crossShoreC  # adjust cross-shore coordinate for pressure to xFRF
        for win in runAvgWin:
            # subtract first point (in air) to "zero"/calibrate the signal on this profile the sig
            p_avg = signal.filtfilt(np.ones(win)/win, 1,  pressure[::-1]) - pressure.iloc[-1]
            runningP[win] = p_avg - np.mean(WL['WL'])  # adjust for waterlevel
        # develop offset to adjust IMU derived profiles to pressure  loc at offshore location
        pressureOffset = -p_avg.max()
        # if integrated profile is too high, # re-compute with double bias
        if integratedProfile == True:
            for ii in range(10):
                if (elevationC.max() + pressureOffset > elevationNAVD88.max()) :
                    bias = bias / 2
                    delev = []
                    for i in range(timeStep.shape[0]):
                        delev.append(np.sin(pitch.iloc[i]-bias) * distOverGround[i])
                    elevationC = np.cumsum(delev, dtype=float)  # elevation # Integrate the individual increments
                    elevationC = np.insert(elevationC, 0, 0) # pad first cell w/ zero
                
            elevationC = elevationC[::-1]  # flip shoreward orientation (for traveling towards shore lines)
        speed = speed[::-1]
    else: # profile is derived going offshore
        # first id pressure
        runningP, shorelineP, crossShoreP= {}, {}, {}
        # crossShoreC = crossShoreC + trueShoreLine
        crossShoreP = xFRF.max() - crossShoreC.max() + crossShoreC  # adjust cross-shore coordinate for pressure to xFRF
        for win in runAvgWin:
            # plt.plot(sb.running_mean(crossShoreC, win), -sb.running_mean(pressure, win), label=r"$\bar{P}$".format(win))
            p_avg = signal.filtfilt(np.ones(win)/win, 1,  pressure) - pressure.iloc[0] # subtract first point to "zero"
            # the sig
            runningP[win] = p_avg - np.mean(WL['WL'])  # adjust for Waterlevel
        pressureOffset = -elevationC[-1] - p_avg.max()  # adjusting the offshore integrated solution to colocated to
        # pressure
        if integratedProfile == True:
            for ii in range(10):
                if (elevationC.max() + pressureOffset > elevationNAVD88.max()) :
                    bias = bias / 2
                    delev = []
                    for i in range(timeStep.shape[0]):
                        delev.append(np.sin(pitch.iloc[i]-bias) * distOverGround[i])
                    elevationC = np.cumsum(delev, dtype=float)  # elevation # Integrate the individual increments
                    elevationC = np.insert(elevationC, 0, 0) # pad first cell w/ zero
    data[f'elevation_pressure'] = -runningP[win]
    data[f'xFRF_pressure'] = crossShoreP
    if plotting is not False:
        fs = 12 # fontsize
        fig = plt.figure();
        ax = plt.subplot2grid((5,1), (0,0), rowspan=4)
        ax.set_title(f"Derived Profiles for {data['time'].iloc[0].strftime('%Y-%m-%d')} profile:"
                f" {data['profileNumber'].mean()}", fontsize=fs)
        for win in runAvgWin:
            # plt.plot(sb.running_mean(crossShoreC, win), -sb.running_mean(pressure, win), label=r"$\bar{P}$".format(win))
            ax.plot(crossShoreP, -runningP[win], label=r"$\bar{P}_{%s}$" %(str(win))) #_{%s}$'.format(win))
        # for order in [3, 5, 7]:
        #     for freq in [10, 20, 50, 1000]:
        #         a,b = signal.butter(order, freq, 'low', analog=True )
        #         lowpassP = signal.filtfilt(b, a, pressure)
        #         plt.plot(crossShoreC, -lowpassP, '.-', ms =1, label=f'lowPass freq:{freq}, order:{order}')
        ax.plot(truthX, depthTruth, '.', label='Crawler GPS')
        if integratedProfile == True:
            ax.plot(crossShoreP, elevationC + pressureOffset, '.', label=r'IMU drift removed: ${%f}^{\circ}$' %(
                    bias))
        xlims = ax.get_xlim()
        plt.plot(xlims, [WL['WL'], WL['WL']], 'b--', label='Water Level' )
        
        ax.legend(fontsize=fs)
        plt.ylabel('Elevation NAVD88 [m]', fontsize=fs)
        
        ax2 = plt.subplot2grid((5,1), (4,0), sharex=ax)
        ax2.plot(crossShoreP, speed, linestyle='-', linewidth=1, label='speed')
        ax2.set_ylabel('speed [m/s]', fontsize=fs)
        ax2.set_xlabel('xFRF [m]', fontsize=fs)
        ax2.tick_params(axis='both', which='both', labelsize=fs, size=5)
        ax.tick_params(axis='y', which='both', labelsize=fs)
        ax.tick_params(
                axis='x',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                bottom=False,      # ticks along the bottom edge are off
                top=False,         # ticks along the top edge are off
                labelbottom=False)
        plt.tight_layout(h_pad=0)

        if plotting is not True:
            plt.savefig(plotting)
            plt.close()
    
    return data