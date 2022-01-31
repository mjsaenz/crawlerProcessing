""" This code is just a general development script for doing a single day/profile comparison"""
import math
import sys, os
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
dateString = '20211019' #'20211020'
getdata = False
# mast [86+97 inch] + antenna centroid [8.42 cm] + deck to floor [32 cm] - tread height [1 in]
offset = 4.6482 + 0.0843 + 0.32 - 0.0254
yMin = 0 # used for defining which data to "keep"
yMax = 1000  # then subset
yRange = 10  # in meters distance in alongshore to consider points "valid" for comparion
savePath = "plots/DUNEXcomparisons"

# flist = glob.glob(f'/data/{dateString}/*GPS_STAT_2.csv')
flist = ['/data/20211019/20211019_175717.075_telemetry.gssbin_GPS_STAT_2.csv']
for GPSfname in flist:
    print(f"\n\nWorking on {dateString}\n\n")
    data = crawlerTools.loadAndMergePriorityFiles(GPSfname, verbose=False, combineDays=True)
    # data = crawlerTools.loadAndMergeFiles(GPSfname, verbose=False)  # no control of output variable names

    ########################################################################
    # figure out start/end times gather background data
    start = DT.datetime.strptime(os.path.dirname(GPSfname).split('/')[-1], "%Y%m%d")
    end = start + DT.timedelta(days=1)
    go = getDataFRF.getObs(start, end)
    if getdata is True:
        topo = go.getLidarDEM()  # get topo data
        bathy = go.getBathyTransectFromNC(method=0) #get bathy data
        bathy = pd.DataFrame.from_dict(bathy)
    elif os.path.isfile(f'data/{start.strftime("%Y_%m_%d")}_bathytopo.pickle'):
        bathy, topo = pickle.load(open(f'data/{start.strftime("%Y_%m_%d")}_bathytopo.pickle', 'rb'))
    else:
        bathy, topo = None, None
    ###############################################
    ## first load file and correct elipsoid values
    if data is None:
        raise NotImplementedError('Need files to load!')
    data = crawlerTools.cleanDF(data)
    
    # add FRF coords to data
    coords = gp.FRFcoord(data.longitude.to_numpy(), data.latitude.to_numpy())
    data['xFRF'] = coords['xFRF']
    data['yFRF'] = coords['yFRF']
    
    # now rotate translate for orientation
    data_og = crawlerTools.TranslateOnly_Wrong(data.copy(), offset)
    data = crawlerTools.rotateTranslatePoints(data, offset)
    if data is None:
        print(f'skipping {fname} no good data\n--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=')
        continue
    
    # quick comparison plots
    crawlerPlots.bathyEnvalopeComparison(GPSfname, data_og, bathy)
    fname = os.path.join(os.path.dirname(GPSfname), ''.join(os.path.basename(GPSfname).split('.')[0])+
                         "_withLocalObs_XY.png")
    crawlerPlots.bathyPlanViewComparison(fname, data_og, bathy, topo)
    
    ## identify profile lines
    # data = crawlerTools.identifyCrawlerProfileLines(data, angleWindow=25, lineLengthThreshold=40,
    #      consecutivePointThresh=50, fname=os.path.join(os.path.dirname(GPSfname),  ''.join(os.path.basename(
    #             GPSfname).split('.')[0])+f"IdentifyProfileLines.png"))
    #
    # data_og = crawlerTools.identifyCrawlerProfileLines(data_og, angleWindow=25, fname=os.path.join(os.path.dirname(GPSfname),
    #                                 ''.join(os.path.basename(GPSfname).split('.')[0])+f"IdentifyProfileLines_OG.png"))
    # ## angular Window
    bathy = pickle.load(open('/home/spike/repos/crawlerProcessing/LarcData.pickle', 'rb'))
    for profile in bathy['profileNumber'].unique():
        if profile <= 1: continue  #handles panda's error AttributeError: 'UnaryOp' object has no attribute 'evaluate'
        subSetLogic = f'(yFRF <= {profile + yRange}) & (yFRF >={profile - yRange}) & (profiles==True)'
        subB = bathy.query(f'(profileNumber == {profile})')
        crawlerTools.searchPointsInRadius(subB, data)
        subC = data.query(f'(profileNumber <= {profile + yRange}) & (profileNumber >= {profile - yRange})')
        subC_og = data.query(f'(profileNumber <= {profile + yRange}) & (profileNumber >= {profile - yRange})')
        if not subC.empty and not subB.empty:
            fOut = os.path.join(savePath, f"singleProfile_{subC['time'].iloc[0].strftime('%Y_%m_%d')}"
                                          f"_{profile.astype(int)}")
            print(f'makingFile {fOut}')
            stats = crawlerPlots.profileCompare(subC=subC, subB=subB, fname=fOut, subC_og=subC_og)
        else:
            print(f'No crawler data for {profile}')