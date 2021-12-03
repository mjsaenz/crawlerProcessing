from matplotlib import pyplot as plt
import os

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
    plt.suptitle(f'elevation envalope for {fnameBase}')
    plt.plot(data.xFRF, data.elevation_NAVD88_m, 'x', label='crawler')
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
        fname:
        data:
        bathy: a dictionary from getdatatestbed
        topo: topography dictionary from getdatatestbed

    Returns:

    """
    plt.figure()
    if bathy is not None:
        plt.scatter(bathy['xFRF'], bathy['yFRF'], c=bathy['elevation'], vmin=-2, vmax=2, label='survey')
    if topo is not None:
        plt.pcolormesh(topo['xFRF'], topo['yFRF'], np.mean(topo['elevation'], axis=0), vmin=-2, vmax=2,
                       label='topo')
    # cmap = plt.scatter(data.xFRF, data.yFRF, c=data.elevation_NAVD88_m-offset, marker='x', vmin=-2, vmax=2, label='crawler')
    
    cbar.set_label('elevation NAVD88')
    plt.xlabel('xFRF')
    plt.ylabel('yFRF')
    plt.title(f'crawler comparison for {data.time.iloc[0].to_pydatetime().strftime("%Y-%m-%dT%H:%M:%SZ")}')
    plt.legend()
    # plt.xlim([30, 300])
    # plt.ylim([720, 745])
    
    plt.savefig(fname1)
    plt.close()
    cmap = plt.scatter(data.xFRF, data.yFRF, c='k', marker='x', label='crawler')
    cbar = plt.colorbar(cmap)
    fnameBase = os.path.basename(fname).split('.')[0]
    fname1 = os.path.join(os.path.dirname(fname), fnameBase+"_withLocalObs_XY.png")

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
    saveFname = os.path.join(savePath, f'SingleProfileCompare_{profileNumber}.png')
    date = subB['time'][0].date().strftime("%Y-%m-%d")
    title = f"Comparison on {date} of profile number {profileNumber}"

    #############
    plt.figure()
    plt.suptitle(title)
    ax1 = plt.subplot(211)
    ax1.plot(subB['xFRF'], subB['elevation'], '.', label='survey')
    ax1.plot(subC['xFRF'], subC['elevation_NAVD88_m'], '.', label='crawler')
    ax1.legend()
    ax1.set_xlabel('xFRF [m]')
    ax1.set_ylabel('elevation [m]')
    ax1.set_xlim([xStart, xStop])

    ax2 = plt.subplot(223)
    ax2.plot(subB['xFRF'], subB['yFRF'], '.', label='survey')
    ax2.plot(subC['xFRF'], subC['yFRF'], '.', label='crawler')
    ax2.set_xlim([xStart-5, xStop+5])
    ax2.set_xlabel('xFRF')
    ax2.set_ylabel('yFRF')

    ax3 = plt.subplot(224)
    dx =  0.6 #np.min(np.diff(subB['xFRF'].squeeze()).mean(), np.median(np.diff(subC['xFRF'])) )
    newX = np.linspace(xStart, xStop, np.round((xStop-xStart)/dx).astype(int))
    crawlInterp = np.interp(newX, subC['xFRF'], subC['elevation_NAVD88_m'])
    surveyInterp = np.interp(newX, subB['xFRF'].squeeze(), subB['elevation'].squeeze())
    ax3.plot(crawlInterp, surveyInterp, '.')
    ax3.plot([-3, 2], [-3, 2], 'k--')
    stats = sb.statsBryant(surveyInterp, crawlInterp)
    ax3.text(-2.75, 0.5, f"RMSE: {stats['RMSEdemeaned']:.2f}[m]\nbias:{stats['bias']:.2f}[m]")
    ax3.set_xlabel('elevation survey')
    ax3.set_ylabel('elevation crawler')
    plt.tight_layout(rect=[0.01, 0.01, 0.99, 0.95])
    plt.savefig(saveFname)
    plt.close()