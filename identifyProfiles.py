""" This code is just a general development script for doing a single day/profile comparison"""
import math
import sys, os
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
from scipy.cluster import vq
from matplotlib import pyplot as plt
import heapq
###############################
dateString = '20211020' #'20211020'
getdata = False
# mast [86+97 inch] + antenna centroid [8.42 cm] + deck to floor [32 cm] - tread height [1 in]
offset = 4.6482 + 0.0843 + 0.32 - 0.0254
yMin = 0 # used for defining which data to "keep"
yMax = 1000  # then subset
yRange = 10  # in meters distance in alongshore to consider points "valid" for comparion
savePath = "plots/DUNEXcomparisons"
searchRadius = 15
# flist = glob.glob(f'/data/{dateString}/*GPS_STAT_2.csv')
flist = ['/data/20211019/20211019_175717.075_telemetry.gssbin_GPS_STAT_2.csv']
flist = ["/home/spike/data/20211020/20211020_153003.335_telemetry.gssbin_GPS_STAT_2.csv"]
GPSfname = "/home/spike/data/20211019/20211019_174456.103_telemetry.gssbin_GPS_STAT_2.csv"
GPSfname = "/home/spike/data/20211020/20211020_153003.335_telemetry.gssbin_GPS_STAT_2.csv"
print(f"\n\nWorking on {dateString}\n\n")

# ##############################################
# ## first load file and correct elipsoid values
# data = crawlerTools.loadAndMergePriorityFiles(GPSfname, verbose=False, combineDays=True)
# data = crawlerTools.cleanDF(data)
#
# from sklearn import cluster
# weightedFactor = [.01, .8, .1, .1] #  @Dylan Anderson @Adam Collins  right here
# for i in [2, 4, 15]:
#     # cluster count
#     clusterCount = i
#     clusterIn = np.array([data['xFRF'].to_numpy(),data['yFRF'].to_numpy(),
#                           np.sin(data['attitude_heading_deg'].to_numpy()),
#                           np.cos(data['attitude_heading_deg'].to_numpy())],
#                          dtype=float).T
#     # data['UNIX_timestamp'].to_numpy()], dtype=float).T
#     whitened = vq.whiten(clusterIn)
#     whitened = whitened * weightedFactor
#     kma = cluster.KMeans(n_clusters=clusterCount, n_init=2000).fit(clusterIn)
#     dba = cluster.DBSCAN(eps=1, min_samples=10000).fit(clusterIn)
#     _, groupSize = np.unique(kma.labels_, return_counts=True)
#     _, groupSize = np.unique(dba.labels_, return_counts=True)
#     d_groups = {}
#     for k in range(clusterCount):
#         d_groups[f'{k}'] = np.where(kma.labels_ == k)[0]
#
#     plt.figure()
#     for key in d_groups.keys():
#         plt.plot(data['xFRF'].iloc[d_groups[key]], data['yFRF'].iloc[d_groups[key]], '.', label=key)
#     plt.legend()

bathyFiles= sorted(glob.glob('/home/spike/data/FRF/geomorphology/elevationTransects/survey/*.nc'))
# flist = ["/home/spike/data/20211020/20211020_153003.335_telemetry.gssbin_GPS_STAT_2.csv",
#         "/home/spike/data/20211019/20211019_174456.103_telemetry.gssbin_GPS_STAT_2.csv",
#         "/home/spike/data/20211018/20211018_152949.708_telemetry.gssbin_GPS_STAT_2.csv",
#         "/home/spike/data/20210928/20210928_185238.280_telemetry.gssbin_GPS_STAT_2.csv",
flist = ["/home/spike/data/20211005/20211005_192624.153_telemetry.gssbin_GPS_STAT_2.csv",
         "/home/spike/data/20211025/20211025_160206.428_telemetry.gssbin_GPS_STAT_2.csv"]

for fname in flist:
    crwl = crawlerTools.loadAndMergePriorityFiles(fname, verbose=False)
    crwl = crawlerTools.cleanDF(crwl)
    crwl = crawlerTools.convert2FRF(crwl)
    date = fname.split('/')[-2]
    print(f'RUNNING DAtTE {date}')
    # identify which survey file we're interested in (manually)
    if date in ['20211020']:  #manually derived values
        idxSurvey = 15  # 2021-10-19
        lineNumbers = [-60, -18, 24, 76,  130, 172, 225, 268, 309, 357, 404, 458]
        lineAngles = [78, 260]
    elif date in ['20211019']:
        idxSurvey = 14  # 20211019.nc'
        lineAngles = [71, 252]
        lineNumbers = [597, 731, 778, 870]
        # complicated line numbers: 640, 688, 824, 915, 959, 1000
    elif date in ['20211018']:
        idxSurvey = 14 # 20211019.nc'
        lineAngles = [75, 255]
        lineNumbers = [637]
        # complicated line numbers:  592, 685
    elif date in ['20211005']:
        lineNumbers = [639, 453, 413]
        lineAngles = [72, 249]
        idxSurvey = 10
    elif date in ['20211025']:
        lineAngles = [75, 256]
        lineNumbers = [957, 866, 823, 775, 731, 686]
        # complicated line number: 917
        idxSurvey = 16
    elif date in ['20210928']:
        lineNumbers=[558, 639]
        lineAngles =[81, 255]
        idxSurvey = 8  # 20210928.nc
        
    else:
        idxSurvey = np.argmin(abs(np.array([int(f.split('_')[-1].split('.nc')[0]) - int(date) for f in bathyFiles])))
    # load survey data
    go = getDataFRF.getObs(DT.datetime.strptime(date, "%Y%m%d") - DT.timedelta(days=6),
                           DT.datetime.strptime(date, "%Y%m%d") + DT.timedelta(days=6))
    
    bathy = go.getBathyTransectFromNC(fname=bathyFiles[idxSurvey], forceReturnAll=True)
    countsHeading, binsHeading, _  = plt.hist(crwl['attitude_heading_deg'], bins=50)
    val1, val2 = heapq.nlargest(2, countsHeading)
    countsY, binsY, _ = plt.hist(crwl['yFRF'], bins=50)

    plotFnme = fname.split('_')[0] + "_Comparison_in_XY_withBathy.png"
    crawlerPlots.bathyPlanViewComparison(plotFnme, crwl, bathy=bathy, topo=None, lineNumbers=lineNumbers, plotShow=False)
    

