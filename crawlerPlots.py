from matplotlib import pyplot as plt
import os
import numpy as np
import sys
sys.path.append('')
from testbedutils import sblib as sb

def bathyEnvalopeComparison(GPSfname, data, bathy, **kwargs):
    """makes an envelop for comparison of background data
    
    Args:
        GPSfname: input file name of crawler data
        data: a data frame from crawlerTools.loadCorrectEllipsoid(
        bathy: a dictionary from getdatatestbed

    Keyword Args:
        "fname": the output file name (default: put it in the same folder as the GPSfname)

    Returns:
        None

    """
    fnameBase = os.path.basename(GPSfname).split('.')[0]
    if bathy is None: bathyTime=None
    else: bathyTime = bathy['time'][0].date()
    ## now make plot
    plt.figure()
    plt.suptitle(f"elevation envalope for "
                 f"{data['time'][0].to_pydatetime().strftime('%Y%m%d %H:%M')}\nSurvey Time: "
                 f"{bathyTime}")
    plt.plot(data.xFRF, data.elevation_NAVD88_m, 'x', label='crawler')
    if bathy is not None:
        plt.plot(bathy['xFRF'], bathy['elevation'], '.', label='survey')
    # plt.xlim([0, 250])
    plt.legend()
    plt.xlabel('xFRF')
    plt.ylabel('elevation $NAVD88$ [m]')
    
    fname2 = kwargs.get('fname', os.path.join(os.path.dirname(GPSfname), f'bathyEnvelope_{fnameBase}.png'))
    plt.savefig(fname2)
    plt.close()

def bathyPlanViewComparison(fname, data, bathy, topo, **kwargs):
    """
    
    Args:
        fname:  output file name
        data:
        bathy: a dictionary from getdatatestbed
        topo: topography dictionary from getdatatestbed

    Returns:

    """
    lineNumbers = kwargs.get('lineNumbers', None)
    plotShow = kwargs.get('plotShow', False)
    topoString, surveyString = None, None
    
    if np.size(data) > 0:
        plt.figure(figsize=(8,6))
        if bathy is not None:
            a = plt.scatter(bathy['xFRF'], bathy['yFRF'], c=bathy['elevation'], vmin=-2, vmax=2, label='survey')
            surveyString = f"Survey: {bathy['time'][0].date()}"
        if topo is not None:
            topomap = plt.pcolormesh(topo['xFRF'], topo['yFRF'], topo['elevation_mean'], vmin=-2, vmax=2,
                           label='topo', shading='auto')
            try:
                topoString = f"topo: {topo['time'][0].strftime('%Y%m%d %H:%M')}"
            except(TypeError):
                topoString = f"topo: {topo['time'].strftime('%Y%m%d %H:%M')}"
            # cmap = plt.scatter(data.xFRF, data.yFRF, c=data.elevation_NAVD88_m-offset, marker='x', vmin=-2, vmax=2, label='crawler')
            # cbar = plt.colorbar(topomap)
            # cbar.set_label('elevation NAVD88')
        plt.xlabel('xFRF')
        plt.ylabel('yFRF')
        plt.colorbar(a)
        plt.title(f'crawler comparison crawler Date '
                  f'{data.time.iloc[0].to_pydatetime().strftime("%Y-%m-%dT%H:%M:%SZ")}\n{surveyString} + {topoString}')
        plt.plot(data.xFRF, data.yFRF, '.k', ms=3, label='crawler')
        if lineNumbers is not None:
            plt.plot(np.ones_like(lineNumbers)*80, lineNumbers, 'rX',ms=10, label='Identified Profiles')
        plt.legend()
        ySpan_c = data['yFRF'].max() - data['yFRF'].min()
        ySpan_b = bathy['yFRF'].max() - bathy['yFRF'].min()
        ySpan = np.argmin([ySpan_c, ySpan_b])
        if ySpan == 0:
            plt.ylim([data['yFRF'].min(), data['yFRF'].max()])
        else:
            plt.ylim([bathy['yFRF'].min(), bathy['yFRF'].max()])
        plt.plot()
        plt.xlim([30, 300])
        
        plt.tight_layout()
        plt.savefig(fname)
        if plotShow is False:
            plt.close()
    else:
        print(f"no plot for fname")
    

def singleProfileComparison(savePath,subB, subC):
    """Makes a comparison of a single profile line
    
    Args:
        savePath:
        subB:
        subC:

    Returns:

    """
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
    totalTiltInterpY = np.interp(newX, subC.sort_values(by='xFRF')['xFRF'], np.max([subC.attitude_0,
                                                                                    subC.attitude_1], axis=0))
    #############
    plt.figure()
    plt.suptitle(title)
    ax1 = plt.subplot(211)
    ax1.plot(subB['xFRF'], subB['elevation'], '.', label='survey - raw')
    ax1.plot(subC['xFRF'], subC['elevation_NAVD88_m'], '.', ms=1, label='crawler - raw')
    c = ax1.scatter(newX, crawlInterp, c=totalTiltInterpY, label='crawler - interp')
    ax1.plot(newX, surveyInterp, '.', ms=1, label='survey - interp ')
    cbar = plt.colorbar(c, ax=ax1)
    cbar.set_label('max(pitch,roll)')
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
    c = ax3.scatter(crawlInterp, surveyInterp, c=np.abs(totalTiltInterpY), vmin=0, vmax=7, cmap='bone')
    ax3.plot([-3, 2], [-3, 2], 'k--')
    stats = sb.statsBryant(surveyInterp, crawlInterp)
    ax3.text(-2.75, 0.5, f"RMSE: {stats['RMSEdemeaned']:.2f}[m]\nbias:{stats['bias']:.2f}[m]")
    ax3.set_xlabel('elevation survey')
    ax3.set_ylabel('elevation crawler')
    cbar = plt.colorbar(c)
    cbar.set_label('alonshore residual')
    ax3.set_ylim([subC['elevation_NAVD88_m'].min()-0.5, subC['elevation_NAVD88_m'].max()+0.5])
    ax3.set_xlim([subC['elevation_NAVD88_m'].min()-0.5, subC['elevation_NAVD88_m'].max()+0.5])
    
    plt.tight_layout(rect=[0.01, 0.01, 0.99, 0.95])
    plt.savefig(saveFname)
    plt.close()
    
    return crawlInterp, surveyInterp


def profileCompare(subC, subB, **kwargs):
    """Plotting a comparison with statistics of a single profile line
    
    Args:
        subC: a subset dataframe for a particular profile of interest from the crawler
        subB: a subset dataframe for comparison ground truth, typically crab/larc
        **kwargs:
            fname: save file name
            subC_og: a secondary version of processing (labeled as subC_og)
    Returns:
        statistics dictionary for comparison

    """
    rawCrawler = kwargs.get('rawCrawler', None)
    
    profileNumber = np.unique(subB['profileNumber']).squeeze()
    date = subB['time'].iloc[0].date().strftime("%Y-%m-%d")
    crawlDate = subC['time'].iloc[0].date().strftime("%Y-%m-%d")
    subC_og = kwargs.get('subC_og', None)
    saveFname = kwargs.get('fname', f'SingleProfileCompare_{crawlDate}_{profileNumber}.png')
    plotRaws = kwargs.get('plotRaws', False)
    # first sort things (make sure we're not sorted in time)
    subC = subC.sort_values(by='xFRF', ignore_index=True)#, inplace=True)
    subB = subB.sort_values(by='xFRF', ignore_index=True)# , inplace=True)
    #plot/data bounds
    xStart = np.ceil(max(subC['xFRF'].min(), subB['xFRF'].min()))
    xStop = np.floor(min(max(subC['xFRF']), max(subB['xFRF'])))

    title = f"Comparison with Suvey {date} and\ncrawl date {crawlDate} of profile number {profileNumber}"
    ############ now work through the interpolation  #####
    dx =  0.6       #np.min(np.diff(subB['xFRF'].squeeze()).mean(), np.median(np.diff(subC['xFRF'])) )
    newX = np.linspace(xStart, xStop, np.round((xStop-xStart)/dx).astype(int), endpoint=True)
    crawlInterp = np.interp(newX, subC['xFRF'], subC['elevation_NAVD88_m'])
    surveyInterp = np.interp(newX, subB['xFRF'], subB['elevation'])
    crawlInterpY = np.interp(newX, subC['xFRF'], subC['yFRF'])
    surveyInterpY = np.interp(newX, subB['xFRF'], subB['yFRF'])
    totalPitchInterpX = np.interp(newX, subC['xFRF'], subC.attitude_pitch_deg)
    totalRollInterpY = np.interp(newX, subC['xFRF'], subC.attitude_roll_deg)
    alongshoreResidual = crawlInterpY - surveyInterpY

    if subC_og is not None:
        crawlInterpY_og = np.interp(newX, subC_og.sort_values(by='xFRF')['xFRF'], subC_og.sort_values(by='xFRF')[
            'yFRF'])
        crawlInterp_og = np.interp(newX, subC_og.sort_values(by='xFRF')['xFRF'], subC_og.sort_values(by='xFRF')[
            'elevation_NAVD88_m'])

    # remove flatlined interpolation points here (points with no source data)
    idxDataJump = np.diff(subC.sort_values(by='xFRF')['xFRF']) > 3  # distances above 3m w/o cralwer source data
    for idx in np.argwhere(idxDataJump):
        idx = int(idx)
        jumpDist = np.diff(subC['xFRF'])[idx]  # dist of the interp w/o data in cellCount
        # idx2Remove = np.linspace(idx-1, idx+np.ceil(jumpDist/dx), np.floor(jumpDist/dx).astype(int),
        #                         endpoint=True).astype(int)
        # loc2Remove = subC['xFRF'][idx2Remove]
        idxRemove = (subC['xFRF'][idx] + jumpDist > newX) & (subC['xFRF'][idx] < newX)  # remove indices from below
        # remove those points from interped arrays of interest
        newX = np.delete(newX, idxRemove)
        crawlInterp = np.delete(crawlInterp, idxRemove)
        surveyInterp = np.delete(surveyInterp, idxRemove)
        crawlInterpY = np.delete(crawlInterpY, idxRemove)
        surveyInterpY = np.delete(surveyInterpY, idxRemove)
        totalPitchInterpX = np.delete(totalPitchInterpX, idxRemove)
        totalRollInterpY = np.delete(totalRollInterpY, idxRemove)
        alongshoreResidual = crawlInterpY - surveyInterpY  #
    # calculate stats
    stats = sb.statsBryant(surveyInterp, crawlInterp)

    ####################################################
    surveyMS = 1
    ogMS = 2
    crawlerMS = 5
    yWindow = 15
    ### Now do the plot
    plt.figure(figsize=(12,8))
    plt.suptitle(title)
    ax1 = plt.subplot2grid((2,4), (0,0), colspan=3)  #plt.subplot(211)
    if plotRaws is True:
        #ax1.plot(subB['xFRF'], subB['elevation'], '.', label='survey - raw')
        ax1.plot(subC['xFRF'], subC['elevation_NAVD88_m'], 'kx', ms=surveyMS, label='Crawler - Raw')
    c = ax1.scatter(newX, crawlInterp, c=totalPitchInterpX, s=25, vmin=-8, vmax=8, label='crawler - interp',
                    cmap='Spectral')
    if subC_og is not None:
        ax1.plot(newX, crawlInterp_og, 'b.', ms=ogMS, label='$crawler_{og}$')
    ax1.plot(newX, surveyInterp, 'k.', ms=surveyMS, label='survey-interp')
    cbar = plt.colorbar(c, ax=ax1)
    cbar.set_label('pitch')
    
    ax1.set_xlabel('xFRF [m]')
    ax1.set_ylabel('elevation [m]')
    ax1.set_xlim([xStart-20, xStop+20])
    ax1.set_ylim([subC['elevation_NAVD88_m'].min()-0.5, subC['elevation_NAVD88_m'].max()+0.5])

    ax2 = plt.subplot2grid((2,4), (1,0), colspan=3, sharex=ax1)  #pplt.subplot(223)
    ax2.plot(subB['xFRF'], subB['yFRF'], 'k.',ms=surveyMS, label='survey')
    if plotRaws is True:
        ax2.plot(subC['xFRF'], subC['yFRF'], 'kx', ms=surveyMS, label='Crawler - Raw')
        if subC_og is not None:
            ax2.plot(subC_og['xFRF'], subC_og['yFRF'], 'dC2', ms=surveyMS, label='OG_Raw')
        if rawCrawler is not None:
            ax2.plot(rawCrawler['xFRF'], rawCrawler['yFRF'], 'Xb', ms=ogMS, label='RAW GPS')
    c2 = ax2.scatter(newX, crawlInterpY, c=totalRollInterpY, s=25,vmin=-5, vmax=5, cmap='Spectral', label='crawler')
    # if subC_og is not None:
        # ax2.plot(newX, crawlInterpY_og, 'b.', ms=ogMS, label='$crawler_{og}$')
    ax2.set_xlim([xStart-5, xStop+5])
    ax2.set_xlabel('xFRF')
    ax2.set_ylabel('yFRF')
    ax2.set_ylim([subB['profileNumber'].mean() - yWindow, subB['profileNumber'].mean() + yWindow])
    ax2.legend()
    cbar = plt.colorbar(c2, ax=ax2)
    cbar.set_label('Roll')
    
    ax3 = plt.subplot2grid((2,4), (0, 3), sharey=ax1) # plt.subplot(224)
    c = ax3.scatter(surveyInterp, crawlInterp, c=np.abs(alongshoreResidual), s=crawlerMS, vmin=0, vmax=10,
                    cmap='inferno', label='Translate/Rotate')
    
    # c = ax3.scatter(crawlInterp, surveyInterp, c=np.abs(totalPitchInterpX), vmin=0, vmax=7, cmap='bone')
    ax3.plot([-3, 2], [-3, 2], 'k--')
    ax3.set_xlabel('elevation survey')
    ax3.set_ylabel('elevation crawler')
    cbar = plt.colorbar(c, ax=ax3)
    cbar.set_label('alonshore residual')
    ax3.set_ylim([subC['elevation_NAVD88_m'].min()-0.5, subC['elevation_NAVD88_m'].max()+0.5])
    ax3.set_xlim([subC['elevation_NAVD88_m'].min()-0.5, subC['elevation_NAVD88_m'].max()+0.5])
    
    ax4 = plt.subplot2grid((2,4), (1,3))
    ax4.set_axis_off()
    ax4.text(0, 0.5, f"Profile Statistics:\n\nRMSE: {stats['RMSEdemeaned']:.2f}[m]\nbias:{stats['bias']:.2f}[m]",
             fontsize=12)
    if subC_og is not None:
        c = ax3.scatter(crawlInterp_og, surveyInterp, c=np.abs(alongshoreResidual), marker='x',s=ogMS, vmin=0,
                        vmax=10, cmap='inferno', label='crawler_og')
        # ax3.legend()
        stats_og = sb.statsBryant(surveyInterp, crawlInterp_og)
        ax4.text(0, 0, f'            OG:\nRMSE: {stats_og["RMSEdemeaned"]:.2f}[m]\nbias:{stats_og["bias"]:.2f}[m]')

    plt.tight_layout(rect=[0.01, 0.01, 0.99, 0.98])
    plt.savefig(saveFname)
    plt.close()
    return stats, newX, surveyInterp, crawlInterp, totalPitchInterpX, totalRollInterpY


def speedComparisonPlot(subC, wave, profile, fnameOut, **kwargs):
    ########
    showplot=kwargs.get('showplot', False)
    from scipy import signal
    window= 30
    f, p_xVel = signal.welch(subC['xFRF_velocity'], fs = np.median(np.diff(subC['time']))/np.timedelta64(1,'s'),
                             nperseg=window)
    f, p_yVel = signal.welch(subC['yFRF_velocity'], fs = np.median(np.diff(subC['time']))/np.timedelta64(1,'s'),
                             nperseg=window)
    f, p_zVel = signal.welch(subC['elevation_velocity'], fs = np.median(np.diff(subC['time']))/
                             np.timedelta64(1,'s'), nperseg=window)
    f, p_spd = signal.welch(subC['speed_over_ground_GPS'], fs = np.median(np.diff(subC[                                                                                                              'time']))/
                                                              np.timedelta64(1,'s'), nperseg=window)
    from matplotlib import pyplot as plt
    plt.style.use('seaborn-paper')
    
    plt.figure()
    ax1 = plt.subplot2grid((2,2),(0,0))
    ax1.plot(1/f, p_xVel, label='x Velocity')
    ax1.plot([1/wave['peakf'], 1/wave['peakf']], [0, .5], linewidth=.1)
    ax1.set_title('x spectra')
    ax1.set_xlabel('period [s]')
    
    ax3 = plt.subplot2grid((2,2),(1,0))
    ax3.plot(1/f, p_yVel, label='y Velocity')
    ax3.plot([1/wave['peakf'], 1/wave['peakf']], [0, .5],linewidth=.1)
    ax3.set_title('y spectra')
    ax3.set_xlabel('period [s]')
    ax4 = plt.subplot2grid((2,2),(1,1))
    ax4.plot(1/f, p_zVel, label='z Velocity')
    ax4.plot([1/wave['peakf'], 1/wave['peakf']], [0, .01], linewidth=.1)
    ax4.set_title('z spectra')
    ax4.set_xlabel('period [s]')
    ax2 = plt.subplot2grid((2,2),(0,1))
    a = ax2.scatter(subC['speed_over_ground_GPS'], subC['NAVSoln_speed_over_ground'],
                 c=subC['elevation_NAVD88_m'])
    cbar = plt.colorbar(a)
    cbar.set_label('elevation')
    ax2.plot([0,5], [0,5], 'k--')
    ax2.set_xlabel('GPS based Speed')
    ax2.set_ylabel('NAV Soln based Speed')
    ax2.set_ylim([0, 2])
    ax2.set_xlim([0,2])
    plt.suptitle(f"Velocity Comparison: {wave['time'][0].strftime('%Y-%m-%d')}: profile {profile}\n$H_s$:"
                 f"{wave['Hs'][0]:.1f}, $T_p$: {1/wave['peakf'][0]:.1f}")
    plt.tight_layout()
    plt.savefig(fnameOut)
    if showplot is not True:
        plt.close()