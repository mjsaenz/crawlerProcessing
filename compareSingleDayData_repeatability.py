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

# GPSfname = "/data/20220113/20220113_151649.080_telemetry.gssbin_GPS_STAT_2.csv"
#
# ########################################################################
# # figure out start/end times gather background data
#
# GPSfname = "/home/spike/data/20211019/20211019_174456.103_telemetry.gssbin_GPS_STAT_2.csv"  # 0.69
# GPSfname = "/home/spike/data/20210928/20210928_185238.280_telemetry.gssbin_GPS_STAT_2.csv"  # 0.77
# GPSfname = "/home/spike/data/20211005/20211005_192624.153_telemetry.gssbin_GPS_STAT_2.csv"  # 0.65
# GPSfname = "/home/spike/data/20211025/20211025_160206.428_telemetry.gssbin_GPS_STAT_2.csv"  # 0.76
# GPSfname =  "/home/spike/data/20211018/20211018_152949.708_telemetry.gssbin_GPS_STAT_2.csv" # 1.22
# GPSfname = "/home/spike/data/20211020/20211020_153003.335_telemetry.gssbin_GPS_STAT_2.csv"  # 0.74
#
#
# GPSfname = "/data/20220316/20220316_133050.186_telemetry.gssbin_GPS_STAT_2.csv"
GPSfname = "/data/20220322_backup/20220322_192544.302_telemetry.gssbin_GPS_STAT_2.csv"
###############################################
## first load file and correct elipsoid values
data = crawlerTools.loadAndMergePriorityFiles(GPSfname, verbose=False, combineDays=True)
data = crawlerTools.convert2FRF(data)
data = crawlerTools.cleanDF(data)
data, offset = crawlerTools.identfyMastOffset(data)

start = data['time'].iloc[0].to_pydatetime()
print(f"\n\n Working on date {start}")
end = data['time'].iloc[-1].to_pydatetime()  #DT.timedelta(days=1)
go = getDataFRF.getObs(start, end)
wave = go.getWaveData('8m-array')
print(f"time: {wave['time'][0].strftime('%Y-%m-%d')}, Hs: {wave['Hs'][0]:.2f}, Tp: {1/wave['peakf'][0]:.2f}, "
      f"Dm: {wave['waveDm'][0]:.02f}")

# topo = go.getLidarDEM()  # get topo data
# surveyFname = "/data/FRF/geomorphology/elevationTransects/survey/FRF_geomorphology_elevationTransects_survey_20220314.nc"
# bathy = go.getBathyTransectFromNC(fname=surveyFname, forceReturnAll=False) #get bathy data
# bathy = pd.DataFrame.from_dict(bathy)
fname = None
topo = None

#
# # now rotate translate for orientation
data = crawlerTools.rotateTranslatePoints(data, offset)
if data is None:
    print(f'skipping {fname} no good data\n--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=')
    raise EnvironmentError("no Good Crawler Data")
# fnameSave = 'data/repeatabilityTest2DataCrawler.pickle'
# pickle.dump(data, open(fnameSave, 'wb'))
# data = crawlerTools.transectSelection(data)
fnameSave = 'data/repeatabilityTest2DataCrawler_withNegative70.pickle'
# pickle.dump(data, open(fnameSave, 'wb'))
data=pickle.load(open(fnameSave,'rb'))
# do some logic here to see if files are exisiting of profile numbers are good
# fnameLoad = 'data/RepeatabilityTest2Transects.xlsx'  # test 2
# fnameLoad = 'data/iso3.xlsx' # tst 3
# if os.path.exists(fnameLoad):
#     data = pd.read_excel(fnameLoad)
    # data = pd.read_excel('data/RepeatabilityTest2Transects.xlsx')

    # this is for repeatability test #1
    #data = pd.read_csv("repeatabilitydemo.csv")  # load matthew's output csv with unique profile numbers
    # try:
    #     data['time'] = [DT.datetime.strptime(data['time'][i], "%Y-%m-%d %H:%M:%S") for i in range(len(data['time']))]
    # except TypeError:
    #     pass
# quick comparison plots
# crawlerPlots.bathyEnvalopeComparison(GPSfname, data, bathy)
# fname = os.path.join(os.path.dirname(GPSfname), ''.join(os.path.basename(GPSfname).split('.')[0])+
#                      "_withLocalObs_XY.png")
# crawlerPlots.bathyPlanViewComparison(fname, data, bathy, topo, plotShow=True)
go2 = getDataFRF.getObs(data['time'].iloc[0].to_pydatetime(), data['time'].iloc[-1].to_pydatetime())
wl = go2.getWL()
# only find the profile that is near the crawler data
# idxRepeatability = (bathy['profileNumber'] < 450) & (bathy['profileNumber']> 400)
# bathy = bathy[idxRepeatability]

# plot cross-shore envelope
##########
crawlerMS, surveyMS = 1, 3
########################################################
# interpolate the crawler and bathy lines to a common x step
from scipy import interpolate
xFRFbase = np.arange(50,310, 1)
# f = interpolate.interp1d(bathy['xFRF'], bathy['elevation'], kind='linear',
#                          bounds_error=False, fill_value=np.nan)
# surveyE = f(xFRFbase)
plt.style.use('fivethirtyeight')
plt.style.use('seaborn')

########################################
# ax1 = plt.subplot2grid((3,1), (0,0))

# ax2.plot(bathy['xFRF'], bathy['yFRF'], ms=surveyMS, marker='.', linestyle='', color='k')
# ax2 = plt.subplot2grid((3,1), (1,0), sharex=ax1)
# ax1.plot(xFRFbase, surveyE, 'k', linewidth=0.25, label='Survey')

eleSave, surSave = [], []
for tt, profNum in enumerate(np.unique(data['profileNumber'].dropna())):
    if (profNum != -1) & (not np.isnan(profNum)):
        idx = np.where(data['profileNumber'] == profNum)[0]
        # if profNum < 410: # vechile is traveling seaward
        #     ax2.plot(data['xFRF'].iloc[idx], data['yFRF'].iloc[idx], ms=crawlerMS, color=colors1[tt],
        #              linewidth=1.5)
        # else:
        #     ax2.plot(data['xFRF'].iloc[idx], data['yFRF'].iloc[idx], ms=crawlerMS, color=colors2[tt%10], linewidth=1.5)
        #
        xFRFnew = np.arange(data['xFRF'].iloc[idx].min(), data['xFRF'].iloc[idx].max(), 1)
        f = interpolate.interp1d(data['xFRF'].iloc[idx], data['elevation_NAVD88_m'].iloc[idx], kind='slinear',
                                 bounds_error=False, fill_value=np.nan)
        elev = f(xFRFbase)
        eleSave.append(elev)
        # surSave.append(medProfile) #surveyE)
# make the median profile
medProfile = np.nanmedian(np.array(eleSave), axis=0)
a_special = 0.25  # constants from
a_exclusive = 0.15
b = 0.0075
d = -medProfile
TVU_special = np.sqrt(a_special**2 + (b*d)**2)
TVU_exclusive = np.sqrt(a_exclusive**2 + (b*d)**2)

####### make plot
fig, axs = plt.subplots(4, 1, constrained_layout=True, figsize=[12,8])
ax1, ax2, ax3, ax4 = axs[0], axs[1], axs[2], axs[3]
ax1.plot(xFRFbase, np.ones_like(xFRFbase) * wl['predictedWL'].mean(), 'b', linestyle=':',linewidth=2, label='Water '
                                                                                                            'Level')
_, _, _, m = ax1.hist2d(data['xFRF'], data['elevation_NAVD88_m'],bins=[np.arange(91, 310, 1), np.arange(-4, 2, 0.1) ],
                        cmap='YlOrBr', norm=LogNorm())
ax1.tick_params(labelbottom=False, bottom=False)
ax2.tick_params(bottom=False, labelbottom=False)
ax3.tick_params(bottom=False, labelbottom=False)
ax1.legend()
cbar = fig.colorbar(mappable=m, ax=ax1, location='top', use_gridspec=True)
cbar.set_label('counts')
# ax3 = plt.subplot2grid((3,1), (2,0), sharex=ax1)
ax3.plot(xFRFbase, np.zeros_like(xFRFbase), 'k:', label='median')
colors1 = plt.cm.winter_r(np.linspace(0, 1, int(len(data['profileNumber'].dropna().unique())/int(2))))
colors2 = plt.cm.autumn_r(np.linspace(0, 1, int(len(data['profileNumber'].dropna().unique())/int(2))))
ax3.legend()
eleSave, surSave = [], []
for tt, profNum in enumerate(np.unique(data['profileNumber'].dropna())):
    if (profNum != -1) & (not np.isnan(profNum)):
        idx = np.where(data['profileNumber'] == profNum)[0]
        if profNum < 410: # vechile is traveling seaward
            ax2.plot(data['xFRF'].iloc[idx], data['yFRF'].iloc[idx], ms=crawlerMS, color=colors1[tt],
                     linewidth=1.5)
        else:
            ax2.plot(data['xFRF'].iloc[idx], data['yFRF'].iloc[idx], ms=crawlerMS, color=colors2[tt%10], linewidth=1.5)

        xFRFnew = np.arange(data['xFRF'].iloc[idx].min(), data['xFRF'].iloc[idx].max(), 1)
        f = interpolate.interp1d(data['xFRF'].iloc[idx], data['elevation_NAVD88_m'].iloc[idx], kind='slinear',
                                 bounds_error=False, fill_value=np.nan)
        elev = f(xFRFbase)
        residualZ = elev - medProfile #surveyE
        if profNum < 410: # vechile is traveling seaward
            ax3.plot(xFRFbase, residualZ, label=f'interp_{profNum}', color=colors1[tt], linewidth=1.5)
        else:
            ax3.plot(xFRFbase, residualZ, label=f'interp_{profNum}', color=colors2[tt%10],
                     linewidth=1.5)

        eleSave.append(elev)
        # surSave.append(medProfile) #surveyE)


ax4.plot(xFRFbase, 2*np.nanstd(np.array(eleSave)[0:20:2], axis=0), color='C1', label='north cluster', linewidth=1.5)
ax4.plot(xFRFbase, 2*np.nanstd(np.array(eleSave)[1:20:2], axis=0), color='b', label='south cluster', linewidth=1.5)
ax4.plot(xFRFbase, TVU_special, label='special-IHO', linewidth=1.5, linestyle='--')
ax4.plot(xFRFbase, TVU_exclusive, label='exclusive-IHO', linewidth=1.5, linestyle='--')
ax4.legend()
# ax4.plot(xFRFbase, np.nanmean(np.array(eleSave) - np.array(surSave), axis=0), label='bias')

ax1.set_ylabel('elevation [m]')
ax1.set_xlim([90, 325])
ax1.set_ylim([-4.5, 2])
ax2.set_ylabel('yFRF [m]')
ax2.set_xlim([90, 325])
ax3.set_xlim([90, 325])
ax4.set_ylabel('uncertainty [m]')
ax4.set_xlim([90, 325])
ax3.set_ylabel('residual [m]')
ax4.set_xlabel('xFRF [m]')
fig.axes[0].text(92.5, 1, 'a)')
fig.axes[1].text(92.5, 415, 'b)')
fig.axes[2].text(92.5, 0.12, 'c)')
fig.axes[3].text(92.5, 0.2, 'd)')
# from testbedutils import sblib as sb
# stats = sb.statsBryant(surSave, eleSave)
plt.savefig('plots/PaperInfoPlots/repeatability.eps', format='eps')
############################################################




#
#
# data = crawlerTools.identifyCrawlerProfileLines(data, angleWindow=25, lineLengthThreshold=40,
#       consecutivePointThresh=50, fname=os.path.join(os.path.dirname(GPSfname),  ''.join(os.path.basename(
#              GPSfname).split('.')[0])+f"IdentifyProfileLines.png"), lineNumbers=lineNumbers, lineAngles=lineAngles,
#                                                 lineWindow=lineWindow)
# #
# # data_og = crawlerTools.identifyCrawlerProfileLines(data_og, angleWindow=25, fname=os.path.join(os.path.dirname(GPSfname),
# #                                 ''.join(os.path.basename(GPSfname).split('.')[0])+f"IdentifyProfileLines_OG.png"))
# for profile in sorted(bathy['bathyprofileNumber'].unique()):
#     # if profile <= 1: continue  #handles panda's error AttributeError: 'UnaryOp' object has no attribute 'evaluate'
#     # subSetLogic = f'(yFRF <= {profile + yRange}) & (yFRF >={profile - yRange}) & (profiles==True)'
#     subB = bathy.query(f'(bathyprofileNumber == {profile})')
#     subC = crawlerTools.searchPointsInRadius(subB, data, radius=searchRadius, removeDuplicates=False,
#                                              searchOnlyLinePoints=True)
#     # subC = data.query(f'(bathyprofileNumber <= {profile + yRange}) & (bathyprofileNumber >= {profile - yRange})')
#     # subC_og = crawlerTools.searchPointsInRadius(subB, data_og, radius=searchRadius, removeDuplicates=False)
#     #data.query(f'('bathyprofileNumber <= {profile + yRange}) & (bathyprofileNumber >= {profile - yRange})')
#     if (subC is not None and subB is not None) and subC.__len__() > 1  and (type(subC) is pd.core.frame.DataFrame
#                               and not subC.empty) and (type(subB) is pd.core.frame.DataFrame and not subB.empty):
#         fOut = os.path.join(savePath, f"singleProfile_{subC['time'].iloc[0].strftime('%Y_%m_%d')}"
#                                       f"_{profile.astype(int):04}")
#         print(f'    makingFile comparison  {fOut}')
#         stats, newX, surveyInterp, crawlInterp, pitch, roll = crawlerPlots.profileCompare(subC=subC, subB=subB,
#                                                         fname=fOut, plotRaws=True)# subC_og=subC_og)
#         logStats.append((GPSfname.split('/')[4], profile, stats, newX, surveyInterp, crawlInterp, pitch, roll))
#     else:
#         print(f'No crawler data for survey {profile}')
#
#
# print('Do something with LogStats')
# pickle.dump(logStats, open("logStats.pkl", 'wb'))
# for i in range(len(logStats)):
# print(f" date {logStats[i][0]}, bathyprofileNumber {logStats[i][1]}")
#
# # ### interpolate and compare
# # from scipy import interpolate
# # dxy=1
# # xmin = np.ceil(max(data['xFRF'].min(), bathy['xFRF'].min()))
# # xmax = np.floor(min(data['xFRF'].max(), bathy['xFRF'].max()))
# # ymin = np.ceil(max(data['yFRF'].min(), bathy['yFRF'].min()))
# # ymax = np.floor(min(data['yFRF'].max(), bathy['yFRF'].max()))
# #
# # yPoints = np.linspace(xmin, xmax, int(xmax-xmin), dxy)
# # xPoints = np.linspace(ymin, ymax, int(ymax-ymin), dxy)
# # yy, xx = np.meshgrid(xPoints, yPoints)
# # grid_c = interpolate.griddata((data['xFRF'], data['yFRF']), data['elevation_NAVD88_m'], xi=(xx, yy),
# #                               method='spline')
# # grid_b = interpolate.griddata((bathy['xFRF'], bathy['yFRF']), bathy['elevation'], xi=(xx,yy))
# #
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