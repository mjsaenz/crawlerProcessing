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
IMUfname = "data/20211019/20211019_181020.816_telemetry.gssbin_OPENINS_IMU_STAT.csv"
NAVfname = "data/20211019/20211019_181020.816_telemetry.gssbin_OPENINS_NAV_SOLUTION.csv"
GPSfname = "data/20211019/20211019_181020.816_telemetry.gssbin_GPS_STAT_2.csv"
# mast [86+97 inch] + antenna centroid [8.42 cm] + deck to floor [32 cm] - tread height [1 in]
offset = 4.6482 + 0.0843 + 0.32 - 0.0254
yMin = 956  # used for defining which data to "keep"
yMax = 962  # then subset
savePath = "plots/DUNEXcomparisons"
########################################################################
# figure out start/end times gather background data
start = DT.datetime(2021, 10, 19)
end = DT.datetime(2021, 10, 20)
go = getDataFRF.getObs(start, end)
topo = go.getLidarDEM()  # get topo data
bathy = go.getBathyTransectFromNC(method=0) #get bathy data
bathy = pd.DataFrame.from_dict(bathy)
###############################################
## first load file and correct elipsoid values
print(f'working on {start}')
data = crawlerTools.loadCorrectEllipsoid(GPSfname, geoidFile='data/g2012bu8.bin', plot=False)
data = pd.merge(data, pd.read_csv(IMUfname, header=4), how='left', on="UNIX_timestamp")
data = pd.merge(data, pd.read_csv(NAVfname, header=4), how='left', on="UNIX_timestamp")
data = crawlerTools.cleanDF(data)

# add FRF coords to data
coords = gp.FRFcoord(data.longitude.to_numpy(), data.latitude.to_numpy())
data['xFRF'] = coords['xFRF']
data['yFRF'] = coords['yFRF']
data = crawlerTools.rotateTranslateAntenna2Ground(data, offset)
print("don't forget to account for pitch/roll")
crawlerPlots.bathyEnvalopeComparison(GPSfname, data, bathy)


########### find subset to focus on ########3
subSetLogic = f'(yFRF <= {yMax}) & (yFRF >={yMin})'
subB = bathy.query()
# sb.reduceDict(bathy, np.argwhere((bathy['profileNumber'] >= yMin) & (bathy['profileNumber'] <=Max)).squeeze(), exemptList=[])
subC = data.query()

#####################################################
