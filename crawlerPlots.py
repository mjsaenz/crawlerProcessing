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

    ## now make plot
    plt.figure()
    plt.suptitle(f"elevation envalope for "
                 f"{data['time'][0].to_pydatetime().strftime('%Y%m%d %H:%M')}\nSurvey Time: {bathy['time'][0].date()}")
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

def bathyPlanViewComparison(fname, data, bathy, topo):
    """
    
    Args:
        fname:  output file name
        data:
        bathy: a dictionary from getdatatestbed
        topo: topography dictionary from getdatatestbed

    Returns:

    """
    topoString, surveyString = None, None
    if np.size(data) > 0:
        plt.figure()
        if bathy is not None:
            plt.scatter(bathy['xFRF'], bathy['yFRF'], c=bathy['elevation'], vmin=-2, vmax=2, label='survey')
            surveyString = f"Survey: {bathy['time'][0].date()}"
        if topo is not None:
            plt.pcolormesh(topo['xFRF'], topo['yFRF'], np.mean(topo['elevation'], axis=0), vmin=-2, vmax=2,
                           label='topo')
            topoString = f"topo: {topo['time'][0].strftime('%Y%m%d %H:%M')}"
            # cmap = plt.scatter(data.xFRF, data.yFRF, c=data.elevation_NAVD88_m-offset, marker='x', vmin=-2, vmax=2, label='crawler')
            cbar = plt.colorbar()
            cbar.set_label('elevation NAVD88')
        plt.xlabel('xFRF')
        plt.ylabel('yFRF')
        
        plt.title(f'crawler comparison crawler Date '
                  f'{data.time.iloc[0].to_pydatetime().strftime("%Y-%m-%dT%H:%M:%SZ")}\n{surveyString} + {topoString}')
        plt.legend()
        # plt.xlim([30, 300])
        # plt.ylim([720, 745])
        plt.plot(data.xFRF, data.yFRF, 'xk', ms=10, label='crawler')
        
        plt.savefig(fname)
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