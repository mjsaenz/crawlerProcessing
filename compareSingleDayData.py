""" This code is just a general development script for doing a single day/profile comparison"""
import sys, os
import crawlerPlots
from testbedutils import geoprocess as gp
from testbedutils import sblib as sb
sys.path.append('/home/spike/repos')
from getdatatestbed import getDataFRF
import crawlerTools
import datetime as DT
from matplotlib import pyplot as plt
import glob
import pandas as pd
import numpy as np

###############################
# IMUfname = "data/20211019/20211019_181020.816_telemetry.gssbin_OPENINS_IMU_STAT.csv"
# NAVfname = "data/20211019/20211019_181020.816_telemetry.gssbin_OPENINS_NAV_SOLUTION.csv"
# GPSfname = "data/20211019/20211019_181020.816_telemetry.gssbin_GPS_STAT_2.csv"
GPSfname = "data/20210928/20210928_175458.243_telemetry.gssbin_GPS_STAT_2.csv"
# NAVfname = "data/20210928/20210928_175458.243_telemetry.gssbin_OPENINS_NAV_SOLUTION.csv"
# IMUfname = "data/20210928/20210928_175458.243_telemetry.gssbin_OPENINS_IMU_STAT.csv"
# mast [86+97 inch] + antenna centroid [8.42 cm] + deck to floor [32 cm] - tread height [1 in]
offset = 4.6482 + 0.0843 + 0.32 - 0.0254
yMin = 956  # used for defining which data to "keep"
yMax = 962  # then subset
yRange = 5  # in meters distance in alongshore to consider points "valid" for comparion
savePath = "plots/DUNEXcomparisons"
########################################################################
# figure out start/end times gather background data
start = DT.datetime(2021, 11, 22)
end = DT.datetime(2021, 11, 22)
# go = getDataFRF.getObs(start, end)
# topo = go.getLidarDEM()  # get topo data
# bathy = go.getBathyTransectFromNC(method=0) #get bathy data
# bathy = pd.DataFrame.from_dict(bathy)
###############################################
## first load file and correct elipsoid values
print(f'working on {start}')
data = crawlerTools.loadAndMergeFiles(GPSfname)
# ldata = crawlerTools.loadCorrectEllipsoid(GPSfname, geoidFile='data/g2012bu8.bin', plot=False)
# NAVdf = pd.read_csv(NAVfname, header=4)
# IMUdf =  pd.read_csv(IMUfname, header=4)oadCorrectEllipsoid(GPSfname, geoidFile='data/g2012bu8.bin', plot=False)
# NAVdf = pd.read_csv(NAVfname, header=4)
# IMUdf =  pd.read_csv(IMUfname, header=4)
# NAVdf = pd.read_csv(NAVfname, header=4)

def interpDataFrames(timeStamp2Interp, df): #  = IMUdf , timeStamp2Interp = data['UNIX_timestamp']
    out = pd.DataFrame()
    for key in df:
        if key != 'UNIX_timestamp':
            try:
                out[key] = np.interp(timeStamp2Interp, df.UNIX_timestamp.astype(float), df[key])
            except TypeError:  # typically strings that are all the same
                if (df[key][0] == df[key]).all():
                    out[key] = df[key][:len(timeStamp2Interp)]
            finally:
                out['UNIX_timestamp'] = timeStamp2Interp
    return out

# NAVdf = interpDataFrames(data.UNIX_timestamp, NAVdf)
# IMUdf = interpDataFrames(data.UNIX_timestamp, IMUdf)
## now merge data
# data = pd.merge(data, IMUdf, how='left', on="UNIX_timestamp")
# data = pd.merge(data, NAVdf, how='left', on="UNIX_timestamp")
data = crawlerTools.cleanDF(data)

# add FRF coords to data
coords = gp.FRFcoord(data.longitude.to_numpy(), data.latitude.to_numpy())
data['xFRF'] = coords['xFRF']
data['yFRF'] = coords['yFRF']
data = crawlerTools.rotateTranslateAntenna2Ground(data, offset)
print("don't forget to account for pitch/roll")
crawlerPlots.bathyEnvalopeComparison(GPSfname, data, bathy)


########### find subset to focus on ########3
for profile in np.unique(bathy.profileNumber):
    subSetLogic = f'(yFRF <= {profile + yRange}) & (yFRF >={profile - yRange})'

    subB = bathy.query(subSetLogic)
    subC = data.query(subSetLogic)

    #####################################################
    if not subC.empty:
        profileComparison = crawlerPlots.singleProfileComparison(savePath, subB, subC)