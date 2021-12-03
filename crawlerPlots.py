from matplotlib import pyplot as plt
import os

from scratchForComparison import data, fname


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