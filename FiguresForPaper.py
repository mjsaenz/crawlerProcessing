""""this scrips assumse a load from """
import pickle
import numpy as np
from matplotlib import pyplot as plt
import os
from testbedutils import sblib as sb
from matplotlib.colors import LogNorm
fpath = "plots/DunexPres"
# load the pickle with my data
logStats = pickle.load(open("logStats.pkl", 'rb'))
# unpack my data
datestring, residualsE, surveyE, newXE, crawlerE, surveyprofileNumber, statsAll, newX, survey, crawler, residuals, \
pitchE, \
rollE, surveyVessel, surveyVesselE, subBAll, subCAll = [], [], [], [], [], [], [],  [], [], [], [], [], [], [], [], \
                                                      [], []
surveyE_purged, crawlerE_purged= [], []
for i in range(len(logStats)):
    datestring.append(logStats[i][0])
    surveyprofileNumber.append(logStats[i][1])
    statsAll.append(logStats[i][2])
    newX.append(logStats[i][3])
    survey.append(logStats[i][4])
    crawler.append(logStats[i][5])
    surveyE.extend(logStats[i][4])
    crawlerE.extend(logStats[i][5])
    print(f" date {datestring[i]}, profileNumber {surveyprofileNumber[i]}")
    residuals.append(crawler[i] - survey[i])
    residualsE.extend(crawler[i] - survey[i])
    newXE.extend(logStats[i][3])
    pitchE.extend(logStats[i][4])
    rollE.extend(logStats[i][5])
    surveyVessel.append(logStats[i][8])
    surveyVesselE.extend(np.tile(logStats[i][8],len(logStats[i][4])))
    subCAll.append(logStats[i][9])
    subBAll.append(logStats[i][10])
    if (logStats[i][0], logStats[i][1]) not in [('20211005', 411), ('20211020', 91), ('20211025', 731)]:
        surveyE_purged.extend(logStats[i][4])
        crawlerE_purged.extend(logStats[i][5])
#calculate overall statistis
statsT = sb.statsBryant(surveyE, crawlerE)
stats_purged = sb.statsBryant(surveyE_purged, crawlerE_purged)
print(f" bias:{statsT['bias']:.3f}m, RMSE: {statsT['RMSE']:.3f}m, RMSE_demeaned: {statsT['RMSEdemeaned']:.3f}m")
print(f" PURGED\nbias:{stats_purged['bias']:.3f}m, RMSE: {stats_purged['RMSE']:.3f}m, RMSE_demeaned: "
      f"{stats_purged['RMSEdemeaned']:.3f}m")

#########################################
#exploreable plots
plt.style.use('seaborn-paper')
fig, axs = plt.subplots(1,2, constrained_layout=True, figsize=[12,6]) #2grid((2,1), (0,0))
ax1, ax2 = axs[0], axs[1]
_,_,_,mca = ax1.hist2d(surveyE, crawlerE, bins=200, cmap='YlOrBr', norm=LogNorm())
a=fig.colorbar(mca, ax=ax1, location='right')
a.set_label('counts');
ax1.text(-2.75, 2.5, 'a)')
# plt.title('depth comparison')
ax1.plot([-3, 3], [-3,3], color='k', linestyle='-', linewidth=1)
ax1.set_xlabel('survey elevation NAVD88 [m]')
ax1.set_ylabel('crawler elevation NAVD88 [m]')

ax2.hist(residualsE, bins=100);
ax2.text(-0.95, 505, 'b)')
ax2.set_xlabel('residuals [m]')
ax2.set_xlim([-1,1])
plt.savefig(os.path.join('plots/PaperInfoPlots', "SkillComparison.eps"), format='eps')
plt.close()
#####################################################
print(f" bias:{statsT['bias']:.3f}m, RMSE: {statsT['RMSE']:.3f}m, RMSE_demeaned: {statsT['RMSEdemeaned']:.3f}m")


################### single pROFILE ################
### good profile example
profile, datestringToCompare = 777, '20211019'
profile, datestringToCompare =  869, '20211019'
profile, datestringToCompare = 686, '20211024'
######## pre processing
idx = np.argwhere((np.array(datestring) == datestringToCompare) & (np.array(surveyprofileNumber) == profile)).squeeze()
subC = subCAll[idx]
subB = subBAll[idx]
stats = statsAll[idx]
xStart = max(subC['xFRF'].min(), subB['xFRF'].min())
xStop = min(max(subC['xFRF']), max(subB['xFRF']))
profileNumber = np.unique(subB['profileNumber']).squeeze()
date = subB['time'].iloc[0].date().strftime("%Y-%m-%d")
saveFname = os.path.join('plots/PaperInfoPlots', f'goodProfile_{date}_{profileNumber}.eps')
# title = f"Comparison on {date} of profile number {profileNumber}"
############3
dx =  1       #np.min(np.diff(subB['xFRF'].squeeze()).mean(), np.median(np.diff(subC['xFRF'])) )
newX = np.linspace(xStart, xStop, np.round((xStop-xStart)/dx).astype(int), endpoint=True)
crawlInterp = np.interp(newX, subC['xFRF'], subC['elevation_NAVD88_m'])
surveyInterp = np.interp(newX, subB['xFRF'], subB['elevation'])
crawlInterpY = np.interp(newX, subC['xFRF'], subC['yFRF'])
surveyInterpY = np.interp(newX, subB['xFRF'], subB['yFRF'])
totalPitchInterpX = np.interp(newX, subC['xFRF'], subC.attitude_pitch_deg)
totalRollInterpY = np.interp(newX, subC['xFRF'], subC.attitude_roll_deg)
alongshoreResidual = crawlInterpY - surveyInterpY

#############
surveyMS = 5
crawlerMS = 12
yWindow = 15
############
plt.figure(figsize=[14,5])
ax1 = plt.subplot2grid((1,4), (0,0), colspan=3)  #plt.subplot(211)
c = ax1.scatter(newX, crawlInterp, c=totalPitchInterpX, s=25, #vmin=-8, vmax=8,
                label='crawler',
                cmap='Spectral')
ax1.plot(newX, surveyInterp, 'k.', ms=surveyMS, label='survey')
cbar = plt.colorbar(c, ax=ax1)
cbar.set_label('pitch [$\degree$]')
ax1.legend()
ax1.set_xlabel('xFRF [m]')
ax1.set_ylabel('elevation [m]')
ax1.set_xlim([90, 225])
ax1.set_ylim([-2.25, 1])
ax1.text(95, 0.75, 'a)')
# ax2 = plt.subplot2grid((2,4), (1,0), colspan=3, sharex=ax1)  #pplt.subplot(223)
# ax2.plot(subB['xFRF'], subB['yFRF'], 'k.',ms=surveyMS, label='survey')
# c2 = ax2.scatter(newX, crawlInterpY, c=totalRollInterpY, s=25,vmin=-5, vmax=5, cmap='Spectral', label='crawler')
# ax2.set_xlim([xStart-5, xStop+5])
# ax2.set_xlabel('xFRF')
# ax2.set_ylabel('yFRF')
# ax2.set_ylim([subB['profileNumber'].mean() - yWindow, subB['profileNumber'].mean() + yWindow])
# ax2.legend()
# cbar = plt.colorbar(c2, ax=ax2)
# cbar.set_label('Roll [$\degree$]')

ax3 = plt.subplot2grid((1,4), (0, 3)) # plt.subplot(224)
c = ax3.scatter(surveyInterp, crawlInterp, c=np.abs(alongshoreResidual), s=10,# vmin=0, vmax=10,
                cmap='inferno', label='Translate/Rotate')

ax3.plot([-3, 2], [-3, 2], 'k--')
ax3.set_xlabel('elevation survey')
ax3.set_ylabel('elevation crawler')
cbar = plt.colorbar(c, ax=ax3)
cbar.set_label('alonshore residual $[m]$')
ax3.set_ylim([-2.25, 0.8])
ax3.set_xlim([-2.25, 0.8])
ax3.text(-2.2, 0.6, 'b)')
# ax4 = plt.subplot2grid((2,4), (1,3))
# ax4.set_axis_off()
# ax4.text(0, 0.5, f"Profile Statistics:\n\nRMSE: {stats['RMSEdemeaned']:.2f}[m]\nbias:{stats['bias']:.2f}[m]",
#          fontsize=12)
stats = sb.statsBryant(surveyInterp,crawlInterp)
print(f"Profile Statistics:\n\nRMSE: {stats['RMSEdemeaned']:.2f}[m]\nbias:{stats['bias']:.2f}[m]")
plt.tight_layout(rect=[0.01, 0.01, 0.99, 0.98])
plt.savefig(saveFname, format='eps')
##########################################################################
##########################################################################
###################example "bad profiles"#################################
##########################################################################
##########################################################################
##########################################################################
#### Bad Data
surveyMS = 5
crawlerMS = 25
profile, datestringToCompare =  91, '20211020'  # 15 cm
profile, datestringToCompare = 731, '20211025'  # 17 cm
profile, datestringToCompare = 411, '20211005'  # 15 cm

idx = np.argwhere((np.array(datestring) == datestringToCompare) & (np.array(surveyprofileNumber) == profile)).squeeze()
subC = subCAll[idx]
subB = subBAll[idx]
stats = statsAll[idx]
xStart = max(subC['xFRF'].min(), subB['xFRF'].min())
xStop = min(subC['xFRF'].max(), max(subB['xFRF']))
profileNumber = np.unique(subB['profileNumber']).squeeze()
date = subB['time'].iloc[0].date().strftime("%Y-%m-%d")
saveFname = os.path.join('plots/PaperInfoPlots', f'badProfileExamples.eps')
from scipy import interpolate
dx =  1       #np.min(np.diff(subB['xFRF'].squeeze()).mean(), np.median(np.diff(subC['xFRF'])) )
newX = np.linspace(xStart, xStop, np.round((xStop-xStart)/dx).astype(int), endpoint=True)
crawlInterp = interpolate.interp1d(subC['xFRF'], subC['elevation_NAVD88_m'])(newX)
surveyInterp = interpolate.interp1d(subB['xFRF'], subB['elevation'])(newX)
crawlInterpY = interpolate.interp1d(subC['xFRF'], subC['yFRF'])(newX)
surveyInterpY = interpolate.interp1d(subB['xFRF'], subB['yFRF'])(newX)
totalPitchInterpX = interpolate.interp1d(subC['xFRF'], subC.attitude_pitch_deg)(newX)
totalRollInterpY = interpolate.interp1d(subC['xFRF'], subC.attitude_roll_deg)(newX)
alongshoreResidual = crawlInterpY - surveyInterpY


### example "bad profiles"
plt.figure(figsize=[14,5])
ax1 = plt.subplot2grid((1,4), (0,0), colspan=3)  #plt.subplot(211)
c = ax1.scatter(newX, crawlInterp, c=totalPitchInterpX, s=crawlerMS, #vmin=-8, vmax=8,
                label='crawler',
                cmap='Spectral')
ax1.plot(newX, surveyInterp, 'k.', ms=surveyMS, label='survey')
cbar = plt.colorbar(c, ax=ax1)
cbar.set_label('pitch [$\degree$]')
ax1.legend()
ax1.set_xlabel('xFRF [m]')
ax1.set_ylabel('elevation [m]')
ax1.set_xlim([105, 265])
ax1.set_ylim([-3, 2.5])
ax1.text(107, 2.25, 'a)')

ax11 = plt.subplot2grid((1,4), (0, 3), sharey=ax1) # plt.subplot(224)
c = ax11.scatter(surveyInterp, crawlInterp, c=np.abs(alongshoreResidual), s=crawlerMS,# vmin=0, vmax=10,
                cmap='inferno', label='Translate/Rotate')

ax11.plot([-3, 3], [-3, 3], 'k--')
ax11.set_xlabel('elevation survey')
ax11.set_ylabel('elevation crawler')
cbar = plt.colorbar(c, ax=ax11)
cbar.set_label('alonshore residual $[m]$')
ax11.set_ylim([-3, 2.5])
ax11.set_xlim(ax11.get_ylim())
ax11.text(-2.8, 2.2, 'b)')
stats = sb.statsBryant(surveyInterp,crawlInterp)
print(f"Profile Statistics:\n\nRMSE: {stats['RMSEdemeaned']:.2f}[m]\nbias:{stats['bias']:.2f}[m]")
plt.tight_layout(rect=[0.01, 0.01, 0.99, 0.98])
plt.savefig(saveFname, format='eps')

## now recompute stats


############### examine if there's signal from crab/larc against crawler #######################

plt.figure()
for ii in range(len(crawlerE)):
    if surveyVesselE[ii] == "CRAB":
        plt.plot(surveyE[ii], crawlerE[ii], 'r.')
    elif surveyVesselE[ii] == 'LARC':
        plt.plot(surveyE[ii], crawlerE[ii], 'b.', ms=3)
    else:
        print('fuckedup!')
plt.ylabel('crawler elevation')
plt.xlabel('survey elevation')
plt.title('CRAB in red, LARC in blue')


# ax2.plot(subB['xFRF'], subB['yFRF'], '.', label='survey')
# ax2.plot(subC['xFRF'], subC['yFRF'], '.', label='crawler')
# ax2.set_xlim([xStart-5, xStop+5])
# ax2.set_xlabel('xFRF')
# ax2.set_ylabel('yFRF')
#
# ax3 = plt.subplot(224)
# # c = ax3.scatter(crawlInterp, surveyInterp, c=np.abs(alongshoreResidual), vmin=0, vmax=7, cmap='bone')
# c = ax3.scatter(crawlInterp, surveyInterp, c=np.abs(totalTiltInterpY), vmin=0, vmax=7, cmap='bone')
# ax3.plot([-3, 2], [-3, 2], 'k--')
# stats = sb.statsBryant(surveyInterp, crawlInterp)
# ax3.text(-2.75, 0.5, f"RMSE: {stats['RMSEdemeaned']:.2f}[m]\nbias:{stats['bias']:.2f}[m]")
# ax3.set_xlabel('elevation survey')
# ax3.set_ylabel('elevation crawler')
# cbar = plt.colorbar(c)
# cbar.set_label('alonshore residual')
# ax3.set_ylim([subC['elevation_NAVD88_m'].min()-0.5, subC['elevation_NAVD88_m'].max()+0.5])
# ax3.set_xlim([subC['elevation_NAVD88_m'].min()-0.5, subC['elevation_NAVD88_m'].max()+0.5])
#
# plt.tight_layout(rect=[0.01, 0.01, 0.99, 0.95])
# plt.savefig(saveFname)
# plt.close()


#
# # binned residuals and 1:1
# plt.figure();
# ax1 = plt.subplot2grid((2,1), (0,0))
# ax1.plot(surveyE, crawlerE, '.')
# ax1.plot([-3, 3], [-3, 3], 'k-')
# ax1.set_xlabel('survey elevation [m]')
# ax1.set_ylabel('crawler elevation [m]')
# ax2 = plt.subplot2grid((2,1), (1,0))
# ax2.hist(residualsE, bins=100);
# ax2.set_xlabel('residuals [m]')
# plt.savefig(os.path.join(fpath, "binnedResiduals"))
# plt.close()


#
# plt.figure();
# plt.hist2d(surveyE, residualsE, cmap='YlOrBr',  bins=100, norm=LogNorm());
# plt.plot([-3,3], [0,0], 'k-');
# plt.xlabel('depth[m]');
# plt.ylabel('residuals');
# a=plt.colorbar();
# a.set_label('counts')
# plt.title('residuals of depth')
# plt.savefig(os.path.join(fpath, "depthResiduals"))
# plt.close()

#
# plt.figure()
# plt.hist2d(newXE, residualsE, bins=100, cmap='YlOrBr', norm=LogNorm())
# a = plt.colorbar()
# a.set_label('counts')
# plt.plot([50, 275], [0,0], 'k-')
# plt.xlabel('xFRF')
# plt.ylabel('residuals')
# plt.title('residuals of cross-shore')
# plt.savefig(os.path.join(fpath, "crossShoreResiduals"))
# plt.close()
#
# plt.figure()
# ax1 = plt.subplot2grid((1,2), (0,0))
# plt.hist2d(pitchE, residualsE, bins=100, cmap='YlOrBr', norm=LogNorm())
# a = plt.colorbar()
# a.set_label('counts')
# plt.plot([-3, 3], [0,0], 'k-')
# plt.xlabel('pitch')
# plt.ylabel('residuals')
# plt.title('residuals of pitch')
#
# ax1 = plt.subplot2grid((1,2), (0,1))
# plt.hist2d(rollE, residualsE, bins=100, cmap='YlOrBr', norm=LogNorm())
# a = plt.colorbar()
# a.set_label('counts')
# plt.plot([-3, 3], [0,0], 'k-')
# plt.xlabel('roll')
# plt.ylabel('residuals')
# plt.title('residuals of roll')
#
# plt.savefig(os.path.join(fpath, "crossShoreResiduals"))
# plt.close()

###############################
allSurvey, allCrawler = pickle.load(open("allCrawlerComparison.pkl", 'rb'))
plt.style.use('seaborn-paper')
plt.figure()

bins, counts,_, surf = plt.hist2d(allCrawler['speed_over_ground_GPS'], allCrawler['NAVSoln_speed_over_ground'],
                                  bins=[np.arange(0, 1, .01), np.arange(0, 1, .01)], cmap='YlOrBr', norm=LogNorm())
cbar = plt.colorbar(surf)
cbar.set_label('counts')
plt.plot([0,1], [0,1], linestyle=':')
plt.ylim([0,1])
plt.xlim([0, 1])
plt.xlabel('measured speed $m/s$')
plt.ylabel('Navigation Solution speed $m/s$')
stats = sb.statsBryant(allCrawler['speed_over_ground_GPS'], allCrawler['NAVSoln_speed_over_ground'])
print(f" RMSE: {stats['RMSE']:.2f} bias: {stats['bias']:.2f}")
plt.savefig('plots/PaperInfoPlots.eps', format='eps')

############# where's the data from
plt.figure();
ax1 = plt.subplot2grid((2,2), (0,0))
bins, counts, _, surf = ax1.hist2d(allCrawler['xFRF'], allCrawler['yFRF'], bins=[np.arange(50, 275, 5), np.arange(50,
                                                                     1000, 5)],  cmap='YlOrBr', norm=LogNorm())
cbar = plt.colorbar(surf, ax=ax1)
cbar.set_label('counts')
ax1.set_xlabel('xFRF $[m]$')
ax1.set_ylabel('yFRF $[m]$')
ax1.set_title('crawler data collection density')

ax2 = plt.subplot2grid((2,2), (0,1))
bins, counts, _, surf = ax2.hist2d(allSurvey['xFRF'], allSurvey['yFRF'], bins=[np.arange(50, 275, 5), np.arange(50,
                                                               1000, 5)],  cmap='YlOrBr', norm=LogNorm())
cbar = plt.colorbar(surf, ax=ax2)
cbar.set_label('counts')
ax2.set_xlabel('xFRF $[m]$')
ax2.set_ylabel('yFRF $[m]$')
ax2.set_title('survey data collection density')
plt.tight_layout()
plt.savefig('plots/PaperInfoPlots/spatialDataset.eps', format='eps')


for key in statsT:
    if key is not 'residuals':
        print(f"Stat: {key}: {statsT[key]:.2f}")


