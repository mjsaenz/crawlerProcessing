""" This code is meant to be used to plot available background data in topo/bathy and what was collected by the
crawler.  It is meant as a first level investigation to show what is available for further investigation."""
import os
import sys
import crawlerPlots
sys.path.append('/home/spike/repos')
from testbedutils import geoprocess as gp
from getdatatestbed import getDataFRF
import crawlerTools
import datetime as DT
import glob
import numpy as np
from testbedutils import sblib as sb

########################################################################
flist = glob.glob("/data/robots/crawlerArchive/2021/*/*GPS_STAT_2*.csv")
yMin = 0
yMax = 1200
offset = 4.6482 + 0.0843 + 0.32 - 0.0254
for fname in sorted(flist):
    print(f'working on {fname}')
    
    ## first load file and correct elipsoid values
    data = crawlerTools.loadCorrectEllipsoid(fname, geoidFile='data/g2012bu8.bin', plot=False)
    data = crawlerTools.cleanDF(data)
    if data is not None and np.size(data) > 0:
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
        
        # get bathy data
        bathy = go.getBathyTransectFromNC(method=0)
        
        # idxBathy = (bathy['profileNumber'] > yMin) & (bathy['profileNumber'] < yMax)
        # bathy = sb.reduceDict(bathy, idxBathy, exemptList=[])
        #
        # bathy=None
        # topo=None
        data = crawlerTools.TranslateOnly_Wrong(data, offset)
        crawlerPlots.bathyEnvalopeComparison(fname, data, bathy)
        
        fname = fname.split('.')[0] + 'withLocalObs_XY.png'
        crawlerPlots.bathyPlanViewComparison(fname, data, bathy, topo)
    
    
    else:
        print(f'--- ERROR: no appropriate data in file {os.path.basename(fname)}')
        ##################### lo