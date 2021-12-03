""" script compares inital crawler demo data"""
import pandas as pd
import netCDF4 as nc
from matplotlib import pyplot as plt
import sys, glob, os
from crawlerTools import closest_node, running_mean, closePathOffset
sys.path.append('/home/spike/repos')
from testbedutils import geoprocess as gp
from testbedutils import sblib as sb
import numpy as np

#################################################
headerCount = 4
unitConversion = 1  # 0.3048  # constant for converting from feet (assumed) to meters  1 means units are the same
savePrefix = "plots"   # location to save data
fileDirectory = '/home/spike/repos/bathyCrawler/data/'    # where nav solution data files are
# load FRF LARC survey data
survey = nc.Dataset('http://bones/thredds/dodsC/FRF/geomorphology/elevationTransects/survey/FRF_20190625_1165_FRF_NAVD88_LARC_GPS_UTC_v20190627.nc')
#################################################
# create list of files to process
pathList = sorted(glob.glob(os.path.join(fileDirectory, '*_telemetry.gssbin_OPENINS_NAV_SOLUTION.csv')))
for path in pathList:
    # load crawler data
    data = pd.read_csv(path, header=headerCount)
    # convert to local coordinate system
    crawlCoords = gp.FRFcoord(np.array(data.absolute_position_0), np.array(data.absolute_position_1), coordType='LL')
    ### identify lines of interest for profile plotting
    idxProfileNumber_out = np.argmin(np.abs(np.median(crawlCoords['yFRF']) - np.unique(survey['profileNumber'])))      # find survey line that crawler crawled
    idxSurveyLine_out = np.argwhere(survey['profileNumber'] == np.unique(survey['profileNumber'])[idxProfileNumber_out]).squeeze()  # find all idx of survey line of interest
    ### take offsets to correct crawler elevation to NAVD88
    startPointFRF = (np.median(crawlCoords['xFRF'][0:10]), np.median(crawlCoords['yFRF'][0:10]))
    argSurvey = closest_node(startPointFRF, (survey['xFRF'], survey['yFRF']))
    # find position nearest start point of crawler from survey data [m NAVD88] to offset crawler data from
    # there was no GPS turned on , no RTK
    zOffset = survey['elevation'][argSurvey]
    # zOffset = 0 # assume that vehicle starts around zero
    crawlerZinNAVD88 = -data.absolute_position_2/unitConversion  + zOffset
    ###  Create smoothed crawler profile line
    smoothed_crawlerZ = running_mean(crawlerZinNAVD88, window=500)
    smoothed_crawlerX = running_mean(crawlCoords['xFRF'], window=500)

    ############## make Plot ##############3
    # plot both surveys in xy
    plt.figure()
    plt.suptitle(path, fontsize=8)
    ax1 = plt.subplot2grid((8,8), (0,0), colspan=8, rowspan=3)
    ax1.plot( survey['xFRF'][:], survey['yFRF'][:],'k.', label='LARC')
    ax1.plot(crawlCoords['xFRF'], crawlCoords['yFRF'], 'r.', label='crawler')
    ax1.legend();
    ax1.set_ylabel('alongshore [m]')
    ax1.set_ylim([200, 500])
    # plot in elevation
    ax2 = plt.subplot2grid((8,8), (3,0), colspan=8, rowspan=5, sharex=ax1)
    ax2.plot(survey['xFRF'][idxSurveyLine_out], survey['elevation'][idxSurveyLine_out], 'k.', label='LARC')
    ax2.plot(crawlCoords['xFRF'][:], crawlerZinNAVD88, 'r.', label='crawler')
    ax2.plot(smoothed_crawlerX, smoothed_crawlerZ, 'b.', label='smoothedCrawler')
    ax2.set_xlabel('cross-shore [m]');
    ax2.set_ylabel('alongshore [m]')
    ax2.set_xlim([80, np.max(crawlCoords['xFRF']+25)])
    figSaveName = os.path.join(savePrefix, '_RAWandSMOOTHED_'.join(path.split('/')[-1].split('.')[0:2]))
    plt.savefig(figSaveName);
    plt.close()

################################################################################################################
################################################################################################################
######## long data comparison  #################################################################################
################################################################################################################
################################################################################################################
############ load crawler data
path = '/home/spike/repos/bathyCrawler/data/20190625_201247.208_telemetry.gssbin_OPENINS_NAV_SOLUTION.csv'
# load crawler data
data = pd.read_csv(path, header=headerCount)
# convert to local coordinate system
crawlCoords = gp.FRFcoord(np.array(data.absolute_position_0), np.array(data.absolute_position_1), coordType='LL')
startIdxLineCrawler = 3200
endIdxLineCrawler = 24300
crawlerX_out = crawlCoords['xFRF'][startIdxLineCrawler:endIdxLineCrawler]
crawlerY_out = crawlCoords['yFRF'][startIdxLineCrawler: endIdxLineCrawler]
crawlerZinNAVD88_out = -data.absolute_position_2[startIdxLineCrawler:endIdxLineCrawler]
smoothed_crawlerZ_out = running_mean(crawlerZinNAVD88_out, window=1000)
smoothed_crawlerX_out = running_mean(crawlerX_out, window=1000)

startIdxLineCrawler = 26960
endIdxLineCrawler = 47860
crawlerX_in = crawlCoords['xFRF'][startIdxLineCrawler:endIdxLineCrawler]
crawlerY_in = crawlCoords['yFRF'][startIdxLineCrawler: endIdxLineCrawler]
crawlerZinNAVD88_in = -data.absolute_position_2[startIdxLineCrawler:endIdxLineCrawler]
smoothed_crawlerZ_in = running_mean(crawlerZinNAVD88_in, window=1000)
smoothed_crawlerX_in = running_mean(crawlerX_in, window=1000)

############3#############3#############3#############3#############3#############3#
### identify lines of interest for survey profile plotting
idxProfileNumber_out = np.argmin(np.abs(
    np.median(crawlerY_out) - np.unique(survey['profileNumber'])))  # find survey line that crawler crawled
idxSurveyLine_out = np.argwhere(survey['profileNumber'] == np.unique(survey['profileNumber'])[
    idxProfileNumber_out]).squeeze()  # find all idx of survey line of interest
idxMaxSurvey = np.argmin(np.abs(survey['xFRF'][idxSurveyLine_out] - max(smoothed_crawlerX_out)))
idxMinSurvey = np.argmin(np.abs(survey['xFRF'][idxSurveyLine_out] - min(smoothed_crawlerX_out)))
surveyX_out = survey['xFRF'][idxSurveyLine_out][idxMinSurvey:idxMaxSurvey]
surveyZ_out = survey['elevation'][idxSurveyLine_out][idxMinSurvey:idxMaxSurvey]

idxProfileNumber_in = np.argmin(np.abs(
    np.median(crawlerY_in) - np.unique(survey['profileNumber'])))  # find survey line that crawler crawled
idxSurveyLine_in = np.argwhere(survey['profileNumber'] == np.unique(survey['profileNumber'])[
    idxProfileNumber_in]).squeeze()  # find all idx of survey line of interest
idxMaxSurvey = np.argmin(np.abs(survey['xFRF'][idxSurveyLine_in] - max(smoothed_crawlerX_in)))
idxMinSurvey = np.argmin(np.abs(survey['xFRF'][idxSurveyLine_in] - min(smoothed_crawlerX_in)))
surveyX_in = survey['xFRF'][idxSurveyLine_out][idxMinSurvey:idxMaxSurvey]
surveyZ_in = survey['elevation'][idxSurveyLine_out][idxMinSurvey:idxMaxSurvey]

from scipy import interpolate
## interpolate to equal spacing for comparison
interpX = max(np.median(np.diff(smoothed_crawlerX_in)), np.median(np.diff(surveyX_out)))
f= interpolate.interp1d(smoothed_crawlerX_in, smoothed_crawlerZ_in, bounds_error=False,fill_value=np.nan)
smoothed_crawlerZ_in = f(surveyX_in)
smoothed_crawlerZ_out = np.interp(surveyX_out, smoothed_crawlerX_out, smoothed_crawlerZ_out)

# calculate skill stats
stats_in = sb.statsBryant(surveyZ_in, smoothed_crawlerZ_in)
stats_out = sb.statsBryant(surveyZ_out, smoothed_crawlerZ_out)
############################################################################
plt.figure();
plt.suptitle('Subaqueous profile comparison', fontsize=12)
ax0 = plt.subplot2grid((8, 8), (0, 0), colspan=8, rowspan=3)
ax0.plot(survey['xFRF'][:], survey['yFRF'][:], 'k.', label='LARC')
ax0.plot(crawlCoords['xFRF'], crawlCoords['yFRF'], 'r.', label='crawler')
ax0.set_ylim([0,400])
ax0.set_xlim([0, 750])

ax1 = plt.subplot2grid((8,8), (3,0), colspan=3, rowspan=5)
ax1.plot(crawlerX_out, crawlerZinNAVD88_out, '.',ms=2, label='Raw Crawler')
ax1.plot(surveyX_out, surveyZ_out, 'k.', label='LARC survey')
ax1.plot(surveyX_out, smoothed_crawlerZ_out, 'r.', label='Smoothed Crawler')
ax1.legend()

ax2 = plt.subplot2grid((8,8), (3,3), sharey=ax1, colspan=3, rowspan=5)
ax2.plot(crawlerX_in, crawlerZinNAVD88_in, '.', ms=2, label='Raw Crawler')
ax2.plot(surveyX_in, surveyZ_in, 'k.', label='Larc Survey')
ax2.plot(surveyX_in, smoothed_crawlerZ_in, 'c.', label='Smoothed Crawler')
ax2.legend()

# one 2 one
ax3 = plt.subplot2grid((8,8), (3,6), sharey=ax1, colspan=3, rowspan=5)
ax3.plot(surveyZ_in, smoothed_crawlerZ_in, 'c.', ms=1)
ax3.plot(surveyZ_out, smoothed_crawlerZ_out, 'r.', ms=1)
ax3.plot([-5.5, 0], [-5.5, 0], 'k-')
ax3.text(-5, -1, 'RETURN:\nbias={:.3f} m\nRMSE={:.3f} m'.format(stats_in['bias'], stats_in['RMSE']), fontsize=12)
ax3.text(-4, -5, 'OUTWARD:\nbias={:.3f} m\nRMSE={:.3f} m'.format(stats_out['bias'], stats_out['RMSE']), fontsize=12)
ax1.set_ylabel('Elevation NAVD88 [m]')
ax1.set_xlabel('Local Cross-shore position [m]')
plt.savefig('Long_survey_out and back comparison.png')
plt.close()
################ plot close loop section
plt.figure();
plt.title('Close Loop Test: GPS denied')
plt.plot(crawlCoords['xFRF'], crawlCoords['yFRF'], '.')
plt.ylabel('FRF x Coordinate')
plt.xlabel('FRF y Coordinate')
plt.xlim([95, 100])
plt.ylim([266, 267])
plt.text(96, 266.8, 'track closure offset {:.2f} m'.format(closePathOffset(data)))
plt.savefig('ClosedLoopTest.png')

# plot both surveys in xy
plt.figure()
plt.suptitle(path, fontsize=8)
ax1 = plt.subplot2grid((8, 8), (0, 0), colspan=8, rowspan=3)
ax1.plot(survey['xFRF'][:], survey['yFRF'][:], 'k.', label='LARC')
ax1.plot(crawlCoords['xFRF'], crawlCoords['yFRF'], 'r.', label='crawler')
ax1.legend();
ax1.set_ylabel('alongshore [m]')
ax1.set_ylim([200, 500])
# plot in elevation
ax2 = plt.subplot2grid((8, 8), (3, 0), colspan=8, rowspan=5, sharex=ax1)
ax2.plot(survey['xFRF'][idxSurveyLine_out], survey['elevation'][idxSurveyLine_out], 'k.', label='LARC')
ax2.plot(crawlCoords['xFRF'][:], crawlerZinNAVD88, 'r.', label='crawler')
ax2.plot(smoothed_crawlerX, smoothed_crawlerZ, 'b.', label='smoothedCrawler')
ax2.set_xlabel('cross-shore [m]');
ax2.set_ylabel('alongshore [m]')
ax2.set_xlim([80, np.max(crawlCoords['xFRF'] + 25)])



figSaveName = os.path.join(savePrefix, '_Refined_'.join(path.split('/')[-1].split('.')[0:2]))
plt.savefig(figSaveName);
plt.close()