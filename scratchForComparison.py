""" This code is just a general development script for playing around with these data"""
import sys
sys.path.append('/home/spike/repos')
from testbedutils import geoprocess as gp
from getdatatestbed import getDataFRF
import crawlerTools
import datetime as DT
from matplotlib import pyplot as plt
import glob
from testbedutils import sblib as sb
import numpy as np



########################################################################
flist = glob.glob('/data/robots/CrawlerFiles/dataCSV/2021/20211019/*GPS_STAT_2.csv')
fname = "/data/robots/CrawlerFiles/GPSOffsetData/20210702_140416.002_telemetry.gssbin_GPS_STAT_2.csv" #flist[2]
offset = 2
yMin = 725  # used for defining which data to "keep"
yMax = 740  # then subset
########################################################################
## first load file and correct elipsoid values
data = crawlerTools.correctEllipsoid(fname, geoidFile='data/g2012bu8.bin')
data = crawlerTools.cleanDF(data)
print("then add rotation and translation of elevations due to mast")
coords = gp.FRFcoord(data.longitude.to_numpy(), data.latitude.to_numpy())
data['xFRF'] = coords['xFRF']
data['yFRF'] = coords['yFRF']

#### subset data to one profile
data = data[(data['yFRF'] > yMin) & (data['yFRF'] < yMax)]

go = getDataFRF.getObs(data.time.iloc[0].to_pydatetime() - DT.timedelta(days=1),
                       data['time'].iloc[-1].to_pydatetime() + DT.timedelta(days=3))
# get topo data
topo = go.getLidarDEM()

#get bathy data
bathy = go.getBathyTransectFromNC(forceReturnAll=True)
idxBathy = (bathy['profileNumber'] > yMin) & (bathy['profileNumber'] < yMax)
bathy = sb.reduceDict(bathy, idxBathy, exemptList = [])

plt.figure()
plt.scatter(bathy['xFRF'], bathy['yFRF'], c=bathy['elevation'], vmin=-2, vmax=2, label='survey')
plt.pcolormesh(topo['xFRF'], topo['yFRF'], np.mean(topo['elevation'], axis=0), vmin=-2, vmax=2, label='topo')
# cmap = plt.scatter(data.xFRF, data.yFRF, c=data.elevation_NAVD88_m-offset, marker='x', vmin=-2, vmax=2, label='crawler')
cmap = plt.scatter(data.xFRF, data.yFRF, c=data['time'], marker='x', label='crawler')
cbar = plt.colorbar(cmap)
cbar.set_label('elevation NAVD88')
plt.xlabel('xFRF')
plt.ylabel('yFRF')
plt.title(f'crawler comparison for {data.time.iloc[0].to_pydatetime().strftime("%Y-%m-%dT%H:%M:%SZ")}')
plt.legend()
plt.xlim([30, 300])
plt.ylim([720, 745])


plt.figure()
plt.plot(data.xFRF, data.elevation_NAVD88_m, 'x', label='crawler')
plt.plot(bathy['xFRF'], bathy['elevation'], '.', label='survey')
plt.xlim([0, 250])
plt.legend()

