""" This code is just a general development script for doing a single day/profile comparison"""
import sys, os

import crawlerPlots
from testbedutils import geoprocess as gp

sys.path.append('/home/spike/repos')
from getdatatestbed import getDataFRF
import crawlerTools
import datetime as DT
import pandas as pd

###############################
# GPSfname = "data/20211019/20211019_181020.816_telemetry.gssbin_GPS_STAT_2.csv"
GPSfname = "data/20210928/20210928_185238.280_telemetry.gssbin_GPS_STAT_2.csv"
#GPSfname = "data/20211005/20211005_191258.345_telemetry.gssbin_GPS_STAT_2.csv"
# mast [86+97 inch] + antenna centroid [8.42 cm] + deck to floor [32 cm] - tread height [1 in]
offset = 4.6482 + 0.0843 + 0.32 - 0.0254
yMin = 0 # used for defining which data to "keep"
yMax = 1000  # then subset
yRange = 10  # in meters distance in alongshore to consider points "valid" for comparion
savePath = "plots/DUNEXcomparisons"
getdata = False
########################################################################
# figure out start/end times gather background data
start = DT.datetime.strptime(GPSfname.split('/')[1], "%Y%m%d")
end = start + DT.timedelta(days=1)
go = getDataFRF.getObs(start, end)
if getdata is True:
    topo = go.getLidarDEM()  # get topo data
    bathy = go.getBathyTransectFromNC(method=0) #get bathy data
    bathy = pd.DataFrame.from_dict(bathy)
else:
    fname = "data/2021_09_28_bathytopo.pickle"
    import pickle
    bathy, topo = pickle.load(open(fname, 'rb'))
###############################################
## first load file and correct elipsoid values
print(f'working on {start}')
data = crawlerTools.loadAndMergeFiles(GPSfname, verbose=False)
data = crawlerTools.cleanDF(data)

# add FRF coords to data
coords = gp.FRFcoord(data.longitude.to_numpy(), data.latitude.to_numpy())
data['xFRF'] = coords['xFRF']
data['yFRF'] = coords['yFRF']

# now rotate translate for orientation
data_og = crawlerTools.TranslateOnly_Wrong(data.copy(), offset)
data = crawlerTools.rotateTranslatePoints(data, offset)

# quick comparison plots
crawlerPlots.bathyEnvalopeComparison(GPSfname, data_og, bathy)
fname = os.path.join(os.path.dirname(GPSfname), ''.join(os.path.basename(GPSfname).split('.')[0])+
                     "_withLocalObs_XY.png")
crawlerPlots.bathyPlanViewComparison(fname, data_og, bathy, topo)

## identify profile lines
data = crawlerTools.identifyCrawlerProfileLines(data, angleWindow=25)
data_og = crawlerTools.identifyCrawlerProfileLines(data_og, angleWindow=25)
## angular Window
for profile in bathy['profileNumber'].unique():
    subSetLogic = f'(yFRF <= {profile + yRange}) & (yFRF >={profile - yRange}) & (profiles==True)'
    subB = bathy.query(f'(profileNumber == {profile})')
    subC = data.query(f'(profileNumber <= {profile + yRange}) & (profileNumber >= {profile - yRange})')
    subC_og = data.query(f'(profileNumber <= {profile + yRange}) & (profileNumber >= {profile - yRange})')
    if not subC.empty and not subB.empty:
        fOut = os.path.join(savePath, f"singleProfile_{subC['time'][0].strftime('%Y_%m_%d')}_{profile.astype(int)}")
        print(f'makingFile {fOut}')
        stats = crawlerPlots.profileCompare(subC=subC, subB=subB, fname=fOut, subC_og=subC_og)