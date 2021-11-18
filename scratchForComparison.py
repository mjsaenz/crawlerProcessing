""" This code is just a general development script for playing around with these data"""
import sys
sys.path.append('/home/spike/repos')
from testbedutils import sblib as sb
from getdatatestbed import getDataFRF
import crawlerTools
## first load file and correct elipsoid values
data = crawlerTools.correctEllipsoid("data/20211020_155001.380_telemetry.gssbin_GPS_STAT_2.csv",
                                     geoidFile='data/g2012bu8.bin')
data = crawlerTools.cleanDF(data)
print("then add rotation and translation of elevations due to mast")

go = getDataFRF.getObs(start, end)