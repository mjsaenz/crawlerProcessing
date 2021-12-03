""" This code is just a general development script for playing around with these data"""
import sys, os

import crawlerPlots
from crawlerPlots import cbar, fname1

sys.path.append('/home/spike/repos')
from testbedutils import geoprocess as gp
from getdatatestbed import getDataFRF
import crawlerTools
import datetime as DT
from matplotlib import pyplot as plt
import glob
import numpy as np

########################################################################
flist = glob.glob("/data/robots/crawlerArchive/2021/*/*GPS_STAT*.csv")
for fname in flist:
    print(f'working on {fname}')
    offset = 0
    # yMin = 450  # used for defining which data to "keep"
    # yMax = 465  # then subset
    ########################################################################
    ## first load file and correct elipsoid values
    data = crawlerTools.loadCorrectEllipsoid(fname, geoidFile='data/g2012bu8.bin', plot=False)
    if data is not None:
        data = crawlerTools.cleanDF(data)
        print("then add rotation and translation of elevations due to mast")
        coords = gp.FRFcoord(data.longitude.to_numpy(), data.latitude.to_numpy())
        data['xFRF'] = coords['xFRF']
        data['yFRF'] = coords['yFRF']
        
        #### subset data to one profile
        # data = data[(data['yFRF'] > yMin) & (data['yFRF'] < yMax)]
        
        go = getDataFRF.getObs(data.time.iloc[0].to_pydatetime() - DT.timedelta(days=1),
                               data['time'].iloc[-1].to_pydatetime() + DT.timedelta(days=3))
        # get topo data
        topo = go.getLidarDEM()
        
        #get bathy data
        bathy = go.getBathyTransectFromNC(method=0)
        # idxBathy = (bathy['profileNumber'] > yMin) & (bathy['profileNumber'] < yMax)
        # bathy = sb.reduceDict(bathy, idxBathy, exemptList = [])
        #
        
        
        crawlerPlots.bathyEnvalopeComparison(fname, data, bathy)
 
    else:
        print(f'--- ERROR: no lat/lon data in file {os.path.basename(fname)}')
        ##################### lo
