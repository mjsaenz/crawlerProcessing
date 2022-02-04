from matplotlib import pyplot as plt
import os
import numpy as np
from testbedutils import sblib as sb

def bathyEnvalopeComparison(fname, data, bathy):
    """makes an envelop for comparison of background data
    
    Args:
        fname: input file name of crawler data
        data: a data frame from crawlerTools.loadCorrectEllipsoid(
        bathy: a dictionary from getdatatestbed

    Returns:
        None

    """
    fnameBase = os.path.basename(fname).split('.')[0]
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
    
    fname2 = os.path.join(os.path.dirname(fname), fnameBase+"_withLocalObs_XZ.png")
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
            plt.pcolormesh(topo['xFRF'], topo['yFRF'], np.mean(topo['elevation'], axis=0), vmin=-2, vmax=2,
                           label='topo', shading='flat')
            topoString = f"topo: {topo['time'][0].strftime('%Y%m%d %H:%M')}"
            # cmap = plt.scatter(data.xFRF, data.yFRF, c=data.elevation_NAVD88_m-offset, marker='x', vmin=-2, vmax=2, label='crawler')
            cbar = plt.colorbar()
            cbar.set_label('elevation NAVD88')
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
    plotRaws = kwargs.get('plotRaws', True)
    #plot/data bounds
    xStart = max(subC['xFRF'].min(), subB['xFRF'].min())
    xStop = min(max(subC['xFRF']), max(subB['xFRF']))

    title = f"Comparison with Suvey {date} and\ncrawl date {crawlDate} of profile number {profileNumber}"
    ############ now work through the interpolation  #####
    dx =  0.6       #np.min(np.diff(subB['xFRF'].squeeze()).mean(), np.median(np.diff(subC['xFRF'])) )
    newX = np.linspace(xStart, xStop, np.round((xStop-xStart)/dx).astype(int))
    crawlInterp = np.interp(newX, subC.sort_values(by='xFRF')['xFRF'], subC.sort_values(by='xFRF')['elevation_NAVD88_m'])
    surveyInterp = np.interp(newX, subB.sort_values(by='xFRF')['xFRF'], subB.sort_values(by='xFRF')['elevation'])
    crawlInterpY = np.interp(newX, subC.sort_values(by='xFRF')['xFRF'], subC.sort_values(by='xFRF')['yFRF'])
    surveyInterpY = np.interp(newX, subB.sort_values(by='xFRF')['xFRF'], subB.sort_values(by='xFRF')['yFRF'])
    
    alongshoreResidual = crawlInterpY - surveyInterpY
    totalPitchInterpY = np.interp(newX, subC.sort_values(by='xFRF')['xFRF'], subC.sort_values(by='xFRF').attitude_pitch_deg)
    totalRollInterpY = np.interp(newX, subC.sort_values(by='xFRF')['xFRF'], subC.sort_values(
            by='xFRF').attitude_roll_deg)

    if subC_og is not None:
        crawlInterpY_og = np.interp(newX, subC_og.sort_values(by='xFRF')['xFRF'], subC_og.sort_values(by='xFRF')[
            'yFRF'])
        crawlInterp_og = np.interp(newX, subC_og.sort_values(by='xFRF')['xFRF'], subC_og.sort_values(by='xFRF')[
            'elevation_NAVD88_m'])
    
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
    c = ax1.scatter(newX, crawlInterp, c=totalPitchInterpY, s=25, vmin=-8, vmax=8, label='crawler - interp',
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
    c = ax3.scatter(crawlInterp, surveyInterp, c=np.abs(alongshoreResidual), s=crawlerMS, vmin=0, vmax=10,
                    cmap='inferno', label='Translate/Rotate')
    
    # c = ax3.scatter(crawlInterp, surveyInterp, c=np.abs(totalPitchInterpY), vmin=0, vmax=7, cmap='bone')
    ax3.plot([-3, 2], [-3, 2], 'k--')
    ax3.set_xlabel('elevation survey')
    ax3.set_ylabel('elevation crawler')
    cbar = plt.colorbar(c, ax=ax3)
    cbar.set_label('alonshore residual')
    ax3.set_ylim([subC['elevation_NAVD88_m'].min()-0.5, subC['elevation_NAVD88_m'].max()+0.5])
    ax3.set_xlim([subC['elevation_NAVD88_m'].min()-0.5, subC['elevation_NAVD88_m'].max()+0.5])
    stats = sb.statsBryant(surveyInterp, crawlInterp)
    
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