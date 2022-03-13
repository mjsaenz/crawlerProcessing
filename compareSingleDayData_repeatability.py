""" This code analyses a single profile repeatedly.  this is for the paper """
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
from matplotlib.colors import LogNorm
###############################
dateString = '20220113' #'20211020'
getdata = True
# mast [86+97 inch] + antenna centroid [8.42 cm] + deck to floor [32 cm] - tread height [1 in]
offset = 4.6482 + 0.0843 + 0.32 - 0.0254
yMin = 0 # used for defining which data to "keep"
yMax = 1000  # then subset
yRange = 10  # in meters distance in alongshore to consider points "valid" for comparion
savePath = "plots/DUNEXcomparisons"
searchRadius = 15  # definition of comparison window between crawler and survey
lineWindow = 10  # defines how much wiggle room to identify a crawler line from above y definitions.

GPSfname = "/data/20220113/20220113_151649.080_telemetry.gssbin_GPS_STAT_2.csv"

########################################################################
# figure out start/end times gather background data
start = DT.datetime(2022, 1, 13)
print(f"\n\n Working on date {start}")
end = start + DT.timedelta(days=1)
go = getDataFRF.getObs(start, end)
wave = go.getWaveData()
print(f"time: {wave['time'][0].strftime('%Y-%m-%d')}, Hs: {wave['Hs'][0]:.2f}, Tp: {1/wave['peakf'][0]:.2f}, "
      f"Dm: {wave['waveDm'][0]:.02f}")

# topo = go.getLidarDEM()  # get topo data
bathy = go.getBathyTransectFromNC(method=0) #get bathy data
bathy = pd.DataFrame.from_dict(bathy)
fname = None
topo = None

###############################################
## first load file and correct elipsoid values
# data = crawlerTools.loadAndMergePriorityFiles(GPSfname, verbose=False, combineDays=True)
# #plt.figure(); plt.hist(data['attitude_roll_deg'], bins=200); plt.title(f'roll {start}'); plt.show()
# plt.figure(); plt.hist(data['attitude_pitch_deg'], bins=200); plt.title(f'pitch {start}'); plt.show()
# plt.show()
# # rawCrawler = data.copy()
# # rawCrawler = crawlerTools.convert2FRF(rawCrawler)
# data = crawlerTools.convert2FRF(data)
# # data = crawlerTools.loadAndMergeFiles(GPSfname, verbose=False)  # no control of output variable names
# if data is None:
#     raise NotImplementedError('Need files to load!')
# data = crawlerTools.cleanDF(data)
#
# # now rotate translate for orientation
# # data_og = crawlerTools.TranslateOnly_Wrong(data.copy(), offset)
# data = crawlerTools.rotateTranslatePoints(data, offset)
# if data is None:
#     print(f'skipping {fname} no good data\n--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=')
#     raise EnvironmentError("no Good Crawler Data")
# fnameSave = 'data/repeatabilityDataCrawler.pickle'
# pickle.dump(data,open(fnameSave, 'wb'))
data = pd.read_csv("repeatabilitydemo.csv")  # load matthew's output csv with unique profile numbers
data['time'] = [DT.datetime.strptime(data['time'][i], "%Y-%m-%d %H:%M:%S") for i in range(len(data['time']))]
# quick comparison plots
crawlerPlots.bathyEnvalopeComparison(GPSfname, data, bathy)
fname = os.path.join(os.path.dirname(GPSfname), ''.join(os.path.basename(GPSfname).split('.')[0])+
                     "_withLocalObs_XY.png")
crawlerPlots.bathyPlanViewComparison(fname, data, bathy, topo, plotShow=True)

# only find the profile that is near the crawler data
idxRepeatability = (bathy['profileNumber'] < 375) & (bathy['profileNumber']> 355)
bathy = bathy[idxRepeatability]
## remove the part of the cralwer data where the vehicle approaches first way point
idxCrawl = data['time'] > DT.datetime(2022, 1, 13, 15, 18, 30)
data = data[idxCrawl]
# plot cross-shore envelope
##########
crawlerMS, surveyMS = 1, 3
########################################################
plt.figure()
ax1 = plt.subplot2grid((2,3), (1,0), colspan=2)
ax1.plot(bathy['xFRF'], bathy['yFRF'], ms=surveyMS, marker='.', linestyle='', color='k')
# ax1.hist2d(data['xFRF'], data['yFRF'], bins=[np.arange(91, 310, 1), np.arange(360, 375, 1)], cmap='YlOrBr',
#            norm=LogNorm())
ax1.plot(data['xFRF'], data['yFRF'], ms=crawlerMS, marker='.', linestyle='', color='r')
ax1.set_ylim([355, 375])
ax2 = plt.subplot2grid((2,3), (0,0), colspan=2, sharex=ax1)
ax2.plot(bathy['xFRF'], bathy['elevation'], ms=surveyMS, marker='.', linestyle='', color='k')
# ax2.plot(data['xFRF'], data['elevation_NAVD88_m'], ms=crawlerMS, marker='.', linestyle='', color='r')
_, _, _, m = ax2.hist2d(data['xFRF'], data['elevation_NAVD88_m'],bins=[np.arange(91, 310, 1), np.arange(-4, 2, 0.1) ],
           cmap='YlOrBr', norm=LogNorm())
cbar = plt.colorbar(mappable=m)
ax2.set_xlim([75, 350])
ax2.set_ylim([-4, 2])
ax3=plt.subplot2grid((2,3), (0,2))

# interpolate the crawler and bathy lines to a common x step
from scipy import interpolate
xFRFbase = np.arange(50,310, 1)
f = interpolate.interp1d(bathy['xFRF'], bathy['elevation'], kind='linear',
                         bounds_error=False, fill_value=np.nan)
surveyE = f(xFRFbase)
########################################
# ax1 = plt.subplot2grid((3,1), (0,0))
fig, axs = plt.subplots(3, 1, constrained_layout=True)
ax1, ax2, ax3 = axs[0], axs[1], axs[2]
ax1.plot(bathy['xFRF'], bathy['yFRF'], ms=surveyMS, marker='.', linestyle='', color='k')
# ax2 = plt.subplot2grid((3,1), (1,0), sharex=ax1)
ax2.plot(xFRFbase, surveyE, 'k', linewidth=0.25, label='Survey')
_, _, _, m = ax2.hist2d(data['xFRF'], data['elevation_NAVD88_m'],bins=[np.arange(91, 310, 1), np.arange(-4, 2, 0.1) ],
                        cmap='YlOrBr', norm=LogNorm())
ax1.tick_params(labelbottom=False, bottom=False)
ax2.tick_params(bottom=False, labelbottom=False)
ax2.legend()
cbar = fig.colorbar(mappable=m, ax=ax2, location='bottom', use_gridspec=True)
cbar.set_label('counts')
# ax3 = plt.subplot2grid((3,1), (2,0), sharex=ax1)
ax3.plot(xFRFbase, np.zeros_like(xFRFbase), 'k:')
eleSave, surSave = [], []
for profNum in np.unique(data['id']):
    if profNum != -1:
        idx = np.where(data['id'] == profNum)[0]
        ax1.plot(data['xFRF'].iloc[idx], data['yFRF'].iloc[idx], ms=crawlerMS)

        xFRFnew = np.arange(data['xFRF'].iloc[idx].min(), data['xFRF'].iloc[idx].max(), 1)
        f = interpolate.interp1d(data['xFRF'].iloc[idx], data['elevation_NAVD88_m'].iloc[idx], kind='linear',
                                 bounds_error=False, fill_value=np.nan)
        elev = f(xFRFbase)
        residualZ = elev - surveyE
        ax3.plot(xFRFbase, residualZ, label=f'interp_{profNum}')
        eleSave.append(elev)
        surSave.append(surveyE)

ax1.set_ylabel('yFRF [m]')
ax1.set_xlim([90, 325])
ax2.set_ylabel('elevation [m]')
ax2.set_xlim([90, 325])
ax3.set_xlim([90, 325])

ax3.set_ylabel('residual')
ax3.set_xlabel('xFRF [m]')
from testbedutils import sblib as sb
stats = sb.statsBryant(surSave, eleSave)
plt.savefig('plots/PaperInfoPlots/repeatability.tif')
############################################################




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