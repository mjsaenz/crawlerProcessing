""" This code analyses a single profile repetedly for the DUNEX time period.  Specficially will outouut a tuple with
statistics/residuals for multiple comparisons  and plots comparing the single profiles. """
import math
import sys, os
import matplotlib.pyplot as plt
import numpy as np
import crawlerPlots
from testbedutils import geoprocess as gp
sys.path.append('/home/spike/repos')
from getdatatestbed import getDataFRF
import crawlerTools
import datetime as DT
import pandas as pd
import glob
import pickle
###############################
dateString = '20211020' #'20211020'
getdata = False
# mast [86+97 inch] + antenna centroid [8.42 cm] + deck to floor [32 cm] - tread height [1 in]
offset = 4.6482 + 0.0843 + 0.32 - 0.0254
yMin = 0 # used for defining which data to "keep"
yMax = 1000  # then subset
yRange = 10  # in meters distance in alongshore to consider points "valid" for comparion
savePath = "plots/DUNEXcomparisons"
searchRadius = 15  # definition of comparison window between crawler and survey
lineWindow = 10  # defines how much wiggle room to identify a crawler line from above y definitions.

# flist = glob.glob(f'/data/{dateString}/*GPS_STAT_2.csv')
flist = ['/data/20211019/20211019_175717.075_telemetry.gssbin_GPS_STAT_2.csv']
flist = ["/home/spike/data/20211020/20211020_153003.335_telemetry.gssbin_GPS_STAT_2.csv"]
# flist = ["/home/spike/data/20211020/20211020_153003.335_telemetry.gssbin_GPS_STAT_2.csv",
#         "/home/spike/data/20211019/20211019_174456.103_telemetry.gssbin_GPS_STAT_2.csv",
#         "/home/spike/data/20211018/20211018_152949.708_telemetry.gssbin_GPS_STAT_2.csv",
#         "/home/spike/data/20210928/20210928_185238.280_telemetry.gssbin_GPS_STAT_2.csv",
flist = ["/home/spike/data/20211005/20211005_192624.153_telemetry.gssbin_GPS_STAT_2.csv",
         "/home/spike/data/20211025/20211025_160206.428_telemetry.gssbin_GPS_STAT_2.csv"]

GPSfname = "/home/spike/data/20211019/20211019_174456.103_telemetry.gssbin_GPS_STAT_2.csv"  # 0.69
GPSfname = "/home/spike/data/20210928/20210928_185238.280_telemetry.gssbin_GPS_STAT_2.csv"  # 0.77
GPSfname = "/home/spike/data/20211005/20211005_192624.153_telemetry.gssbin_GPS_STAT_2.csv"  # 0.65
GPSfname = "/home/spike/data/20211025/20211025_160206.428_telemetry.gssbin_GPS_STAT_2.csv"  # 0.76
GPSfname =  "/home/spike/data/20211018/20211018_152949.708_telemetry.gssbin_GPS_STAT_2.csv" # 1.22
GPSfname = "/home/spike/data/20211020/20211020_153003.335_telemetry.gssbin_GPS_STAT_2.csv"  # 0.74

for GPSfname in flist:
    
    ########################################################################
    # figure out start/end times gather background data
    start = DT.datetime.strptime(os.path.dirname(GPSfname).split('/')[-1], "%Y%m%d")
    print(f"\n\n Working on date {start}")
    end = start + DT.timedelta(days=1)
    go = getDataFRF.getObs(start, end)
    wave = go.getWaveData()
    print(f"time: {wave['time'][0].strftime('%Y-%m-%d')}, Hs: {wave['Hs'][0]:.2f}, Tp: {1/wave['peakf'][0]:.2f}")
    
    if getdata is True:
        topo = go.getLidarDEM()  # get topo data
        bathy = go.getBathyTransectFromNC(method=0) #get bathy data
        bathy = pd.DataFrame.from_dict(bathy)
    
        # elif os.path.isfile(f'data/{start.strftime("%Y_%m_%d")}_bathytopo.pickle'):
    #     bathy, topo = pickle.load(open(f'data/{start.strftime("%Y_%m_%d")}_bathytopo.pickle', 'rb'))
    #
    elif start == DT.datetime(2021, 10, 20, 0, 0):
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
        go.d1 = start - DT.timedelta(days=2)
        go.epochd1 = go.d1.timestamp()
        # go.end = start + DT.timedelta(days=3)
        # go.epochd2 = go.end.timestamp()
        lineAngles = [75, 255]
        lineNumbers = [592, 637]
        # complicated line numbers:  592, 685
        go.end = start + DT.timedelta(days=5)
        go.epochd2 = go.end.timestamp()
        lineWindow = 8  # defines how much wiggle room to identify a crawler line from above y definitions.

    elif start == DT.datetime(2021, 10, 19):
        go.epochd2 = go.d2.timestamp()
        fname = "/data/FRF/geomorphology/elevationTransects/survey/FRF_geomorphology_elevationTransects_survey_20211019.nc"
        go.end = start + DT.timedelta(days=3)
        go.epochd2 = go.end.timestamp()
        idxSurvey = 14  # 20211019.nc'
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
        go.end = start + DT.timedelta(days=3)
        go.epochd2 = go.end.timestamp()
        go.d1 = start - DT.timedelta(days=2)
        go.epochd1 = go.d1.timestamp()
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

    bathy = go.getBathyTransectFromNC(fname=fname, forceReturnAll=True)
    bathy = pd.DataFrame.from_dict(bathy)
    topo = None
    
    ###############################################
    ## first load file and correct elipsoid values
    data = crawlerTools.loadAndMergePriorityFiles(GPSfname, verbose=False, combineDays=True)
    #plt.figure(); plt.hist(data['attitude_roll_deg'], bins=200); plt.title(f'roll {start}'); plt.show()
    plt.figure(); plt.hist(data['attitude_pitch_deg'], bins=200); plt.title(f'pitch {start}'); plt.show()
    plt.show()
    # rawCrawler = data.copy()
    # rawCrawler = crawlerTools.convert2FRF(rawCrawler)
    data = crawlerTools.convert2FRF(data)
    # data = crawlerTools.loadAndMergeFiles(GPSfname, verbose=False)  # no control of output variable names
    if data is None:
        raise NotImplementedError('Need files to load!')
    data = crawlerTools.cleanDF(data)
    
    # now rotate translate for orientation
    data_og = crawlerTools.TranslateOnly_Wrong(data.copy(), offset)
    data = crawlerTools.rotateTranslatePoints(data, offset)
    if data is None:
        print(f'skipping {fname} no good data\n--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=')
        raise EnvironmentError("no Good Crawler Data")
    
    # quick comparison plots
    crawlerPlots.bathyEnvalopeComparison(GPSfname, data, bathy)
    fname = os.path.join(os.path.dirname(GPSfname), ''.join(os.path.basename(GPSfname).split('.')[0])+
                         "_withLocalObs_XY.png")
    crawlerPlots.bathyPlanViewComparison(fname, data, bathy, topo)
    
    ## identify profile linesdata.hist('yFRF', bins=50)
    data = crawlerTools.identifyCrawlerProfileLines(data, angleWindow=25, lineLengthThreshold=40,
          consecutivePointThresh=50, fname=os.path.join(os.path.dirname(GPSfname),  ''.join(os.path.basename(
                 GPSfname).split('.')[0])+f"IdentifyProfileLines.png"), lineNumbers=lineNumbers, lineAngles=lineAngles,
                                                    lineWindow=lineWindow)
    #
    # data_og = crawlerTools.identifyCrawlerProfileLines(data_og, angleWindow=25, fname=os.path.join(os.path.dirname(GPSfname),
    #                                 ''.join(os.path.basename(GPSfname).split('.')[0])+f"IdentifyProfileLines_OG.png"))
    for profile in sorted(bathy['profileNumber'].unique()):
        # if profile <= 1: continue  #handles panda's error AttributeError: 'UnaryOp' object has no attribute 'evaluate'
        # subSetLogic = f'(yFRF <= {profile + yRange}) & (yFRF >={profile - yRange}) & (profiles==True)'
        subB = bathy.query(f'(profileNumber == {profile})')
        subC = crawlerTools.searchPointsInRadius(subB, data, radius=searchRadius, removeDuplicates=False,
                                                 searchOnlyLinePoints=True)
        # subC = data.query(f'(profileNumber <= {profile + yRange}) & (profileNumber >= {profile - yRange})')
        # subC_og = crawlerTools.searchPointsInRadius(subB, data_og, radius=searchRadius, removeDuplicates=False)
        #data.query(f'('profileNumber <= {profile + yRange}) & (profileNumber >= {profile - yRange})')
        if (subC is not None and subB is not None) and subC.__len__() > 1  and (type(subC) is pd.core.frame.DataFrame
                                  and not subC.empty) and (type(subB) is pd.core.frame.DataFrame and not subB.empty):
            fOut = os.path.join(savePath, f"singleProfile_{subC['time'].iloc[0].strftime('%Y_%m_%d')}"
                                          f"_{profile.astype(int):04}")
            print(f'    makingFile comparison  {fOut}')
            stats, newX, surveyInterp, crawlInterp, pitch, roll = crawlerPlots.profileCompare(subC=subC, subB=subB,
                                                            fname=fOut, plotRaws=True)# subC_og=subC_og)
            logStats.append((GPSfname.split('/')[4], profile, stats, newX, surveyInterp, crawlInterp, pitch, roll))
        else:
            print(f'No crawler data for survey {profile}')
    
    
print('Do something with LogStats')
pickle.dump(logStats, open("logStats.pkl", 'wb'))
for i in range(len(logStats)):
    print(f" date {logStats[i][0]}, profileNumber {logStats[i][1]}")
    
    # ### interpolate and compare
    # from scipy import interpolate
    # dxy=1
    # xmin = np.ceil(max(data['xFRF'].min(), bathy['xFRF'].min()))
    # xmax = np.floor(min(data['xFRF'].max(), bathy['xFRF'].max()))
    # ymin = np.ceil(max(data['yFRF'].min(), bathy['yFRF'].min()))
    # ymax = np.floor(min(data['yFRF'].max(), bathy['yFRF'].max()))
    #
    # yPoints = np.linspace(xmin, xmax, int(xmax-xmin), dxy)
    # xPoints = np.linspace(ymin, ymax, int(ymax-ymin), dxy)
    # yy, xx = np.meshgrid(xPoints, yPoints)
    # grid_c = interpolate.griddata((data['xFRF'], data['yFRF']), data['elevation_NAVD88_m'], xi=(xx, yy),
    #                               method='spline')
    # grid_b = interpolate.griddata((bathy['xFRF'], bathy['yFRF']), bathy['elevation'], xi=(xx,yy))
    #
    #
    # plt.figure()
    # ax1 = plt.subplot2grid((1,3), (0,0))
    # a1 = plt.pcolormesh(xx,yy, grid_b, vmin=data['elevation_NAVD88_m'], vmax=data['elevation_NAVD88_m'],
    #                                             cmap='ocean')
    # ax1.plot(bathy['xFRF'], bathy['yFRF'], 'k.', ms=1)
    # plt.colorbar(a1, ax=ax1)
    # plt.ylabel('yFRF')
    # plt.ylim([ymin, ymax])
    # plt.xlim([xmin, xmax])
    #
    # ax2 = plt.subplot2grid((1,3), (0,1), sharey = ax1, sharex=ax1)
    # a2 = ax2.pcolormesh(xx,yy, grid_c, cmap='ocean')
    # plt.colorbar(a2, ax=ax2)
    # ax2.plot(data['xFRF'], data['yFRF'], 'k.', ms=1)
    #
    # ax3 = plt.subplot2grid((1,3), (0,2), sharex=ax1, sharey=ax1)
    # a3 = ax3.pcolormesh(xx, yy, grid_c-grid_b, cmap='RdBu')
    # plt.colorbar(a3, ax=ax3)
    #