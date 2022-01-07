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
# GPSfname = "data/20211019/20211019_181020.816_telemetry.gssbin_GPS_STAT_2.csv"
# GPSfname = "data/20210928/20210928_185238.280_telemetry.gssbin_GPS_STAT_2.csv"
GPSfname = "data/20211005/20211005_191258.345_telemetry.gssbin_GPS_STAT_2.csv"
# mast [86+97 inch] + antenna centroid [8.42 cm] + deck to floor [32 cm] - tread height [1 in]
offset = 4.6482 + 0.0843 + 0.32 - 0.0254
yMin = 956  # used for defining which data to "keep"
yMax = 962  # then subset
yRange = 10  # in meters distance in alongshore to consider points "valid" for comparion
savePath = "plots/DUNEXcomparisons"
########################################################################
# figure out start/end times gather background data
start = DT.datetime(2021, 9, 28)
end = DT.datetime(2021, 9, 29)
go = getDataFRF.getObs(start, end)
topo = go.getLidarDEM()  # get topo data
bathy = go.getBathyTransectFromNC(method=0) #get bathy data
bathy = pd.DataFrame.from_dict(bathy)
###############################################
## first load file and correct elipsoid values
print(f'working on {start}')
data = crawlerTools.loadAndMergeFiles(GPSfname, verbose=False)
data = crawlerTools.cleanDF(data)

# add FRF coords to data
coords = gp.FRFcoord(data.longitude.to_numpy(), data.latitude.to_numpy())
data['xFRF'] = coords['xFRF']
data['yFRF'] = coords['yFRF']
data = crawlerTools.rotateTranslateAntenna2Ground(data, offset)
print("don't forget to account for pitch/roll")
crawlerPlots.bathyEnvalopeComparison(GPSfname, data, bathy)
fname = os.path.join(os.path.dirname(GPSfname), ''.join(os.path.basename(GPSfname).split('.')[0])+
                     "_withLocalObs_XY.png")
crawlerPlots.bathyPlanViewComparison(fname, data, bathy, topo)


########### find subset to focus on ########3
for profile in np.unique(bathy.profileNumber):
    subSetLogic = f'(yFRF <= {profile + yRange}) & (yFRF >={profile - yRange})'
    
    subB = bathy.query(subSetLogic)
    subC = data.query(subSetLogic)
    print(f"found {subC.shape[0]} in range {subSetLogic}")
    #####################################################
    if not subC.empty:
        profileComparison = crawlerPlots.singleProfileComparison(savePath, subB, subC)

###################################3
pitch = data.attitude_pitch_deg
roll = data.attitude_roll_deg
idx = np.argmin(np.abs(subC['xFRF'] - 117))  # pick single point
x = subC['xFRF'].iloc[idx]
y = subC['yFRF'].iloc[idx]
z = subC['elevation_NAVD88_m'].iloc[idx]
pitch_i = pitch.iloc[idx]
roll_i = roll.iloc[idx]
translate = np.ones_like(roll) * offset

T =  np.matrix(([0, 0, 0, 0],
                [0, 0, 0, 0],
                [0, 0, 0, offset],
                [0, 0, 0, 1]))

Rx = np.matrix([[1, 0, 0, 0],
                [0, np.cos(roll.iloc[idx]), -np.sin(roll.iloc[idx]), 0 ],
                [0, np.sin(roll.iloc[idx]), np.cos(roll.iloc[idx]),  0 ],
                [0, 0,                      0,                       1 ],
                ])
Ry = np.matrix(([np.cos(pitch.iloc[idx]),  0,    np.sin(pitch.iloc[idx]), 0],
                [0,                        1,    0,                       0],
                [-np.sin(pitch.iloc[idx]), 0,    np.cos(pitch.iloc[idx]), 0],
                [ 0                      , 0,    1,                       1]))
input = np.matrix([x, y, z, 1])
Ry @ Ry @ T @ input.T


# this is code from the above plotting routine used to analyze one profile
profile = 558 #  fixing the profile number for comparison
xStart = max(subC['xFRF'].min(), subB['xFRF'].min())
xStop = min(max(subC['xFRF']), max(subB['xFRF']))
profileNumber = np.unique(subB['profileNumber']).squeeze()
date = subB['time'].iloc[0].date().strftime("%Y-%m-%d")
saveFname = os.path.join(savePath, f'SingleProfileCompare_{date}_{profileNumber}.png')
title = f"Comparison on {date} of profile number {profileNumber}"
############3
dx =  0.6 #np.min(np.diff(subB['xFRF'].squeeze()).mean(), np.median(np.diff(subC['xFRF'])) )
newX = np.linspace(xStart, xStop, np.round((xStop-xStart)/dx).astype(int))
crawlInterp = np.interp(newX, subC.sort_values(by='xFRF')['xFRF'], subC.sort_values(by='xFRF')['elevation_NAVD88_m'])
surveyInterp = np.interp(newX, subB.sort_values(by='xFRF')['xFRF'], subB.sort_values(by='xFRF')['elevation'])
crawlInterpY = np.interp(newX, subC.sort_values(by='xFRF')['xFRF'], subC.sort_values(by='xFRF')['yFRF'])
surveyInterpY = np.interp(newX, subB.sort_values(by='xFRF')['xFRF'], subB.sort_values(by='xFRF')['yFRF'])
alongshoreResidual = crawlInterpY - surveyInterpY
totalPitchInterpY = np.interp(newX, subC.sort_values(by='xFRF')['xFRF'], subC.sort_values(by='xFRF').attitude_pitch_deg)
totalRollInterpY = np.interp(newX, subC.sort_values(by='xFRF')['xFRF'], subC.sort_values(by='xFRF').attitude_roll_deg)


plt.figure()
plt.suptitle(title)
ax1 = plt.subplot(211)
ax1.plot(subB['xFRF'], subB['elevation'], '.', label='survey - raw')
ax1.plot(subC['xFRF'], subC['elevation_NAVD88_m'], '.', ms=1, label='crawler - raw')
c = ax1.scatter(newX, crawlInterp, c=totalPitchInterpY, vmin=-8, vmax=8, label='crawler - interp', cmap='RdBu')
ax1.plot(newX, surveyInterp, '.', ms=1, label='survey - interp ')
cbar = plt.colorbar(c, ax=ax1)
cbar.set_label('pitch')
ax1.legend()
ax1.set_xlabel('xFRF [m]')
ax1.set_ylabel('elevation [m]')
ax1.set_xlim([xStart-20, xStop+20])
ax1.set_ylim([subC['elevation_NAVD88_m'].min()-0.5, subC['elevation_NAVD88_m'].max()+0.5])

ax2 = plt.subplot(223)
ax2.plot(subB['xFRF'], subB['yFRF'], '.', label='survey')
ax2.plot(subC['xFRF'], subC['yFRF'], '.', label='crawler')
ax2.set_xlim([xStart-5, xStop+5])
ax2.set_xlabel('xFRF')
ax2.set_ylabel('yFRF')

ax3 = plt.subplot(224)
# c = ax3.scatter(crawlInterp, surveyInterp, c=np.abs(alongshoreResidual), vmin=0, vmax=7, cmap='bone')
c = ax3.scatter(crawlInterp, surveyInterp, c=np.abs(totalPitchInterpY), vmin=0, vmax=7, cmap='bone')
ax3.plot([-3, 2], [-3, 2], 'k--')
stats = sb.statsBryant(surveyInterp, crawlInterp)
ax3.text(-2.75, 0.5, f"RMSE: {stats['RMSEdemeaned']:.2f}[m]\nbias:{stats['bias']:.2f}[m]")
ax3.set_xlabel('elevation survey')
ax3.set_ylabel('elevation crawler')
cbar = plt.colorbar(c)
cbar.set_label('alonshore residual')
ax3.set_ylim([subC['elevation_NAVD88_m'].min()-0.5, subC['elevation_NAVD88_m'].max()+0.5])
ax3.set_xlim([subC['elevation_NAVD88_m'].min()-0.5, subC['elevation_NAVD88_m'].max()+0.5])
