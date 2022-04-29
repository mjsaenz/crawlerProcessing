""" This code analyses a single profile repeatedly.  this is for the paper """
import matplotlib
matplotlib.use('TkAgg')
import sys, os
import crawlerPlots
import crawlerTools
import datetime as DT
import pickle
import glob
import py2netCDF

###############################

def crawlerProcessingWorkFlow(path2GPSfname, **kwargs):
    """Post processing for crawler data
    
    Args:
        path2GPSfname:  Folder path to look for bin files
        **kwargs:
            "binProcessor": path to bin executable (default="/usr/local/bin/gsbin-log-processor")
            "pickleFname": name for pickle file
    Returns:

    """
    gsBinProcessorLocation = kwargs.get('binProcessor', "/usr/local/bin/gsbin-log-processor")
    fnameBase = os.path.basename(path2GPSfname)
    baseDir = os.path.dirname(path2GPSfname)
    fnameSave = kwargs.get('pickleFname', os.path.join(baseDir, fnameBase, f'processed_{fnameBase}.pkl'))

    crawlerTools.loadAndConvertGSSfiles(gsBinProcessorLocation, path2GPSfname)

    # now pick out the first file with GPS2.csv for processing the day
    GPSfname = glob.glob(os.path.join(path2GPSfname, "*GPS_STAT_2.csv"))[0]
    ###############################################
    ## first load file and correct elipsoid values
    data = crawlerTools.loadAndMergePriorityFiles(GPSfname, verbose=False, combineDays=True)
    data = crawlerTools.convert2FRF(data)                # make FRF coords
    data = crawlerTools.cleanDF(data)                    # clean the data frame
    data, offset = crawlerTools.identfyMastOffset(data)  # generate the mast offset if we have the base
    if offset is None:
        # there's no benchmark and have to assign an offset
        offset = 4.6482 + 0.0843 + 0.32 - 0.0254
        print(f'used default offset {offset}')

    print(f"\n\n Working on date {fnameBase}")
    # # now rotate translate for orientation
    data = crawlerTools.rotateTranslatePoints(data, offset)
    try:
        data = pickle.load(open(fnameSave,'rb'))
        if 'lineNumber' not in data.keys():
            data = crawlerTools.transectSelection(data)
    except FileNotFoundError:
        data = crawlerTools.transectSelection(data)
        pickle.dump(open(fnameSave, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
    print('  TODO: look for pickle file, if exists, open check for line numbers')
        
    data['time'] = [DT.datetime.strptime(data['time'][i], "%Y-%m-%d %H:%M:%S") for i in range(len(data['time']))]
    
    # quick comparison plots for QA/QC
    crawlerPlots.bathyEnvalopeComparison(GPSfname, data, bathy=None)
    fname = os.path.join(os.path.dirname(GPSfname), os.path.join(baseDir, fnameBase, f'planView_{fnameBase}.png'))
    crawlerPlots.bathyPlanViewComparison(fname, data,  bathy=None, topo=None, plotShow=False)


    ############################################################
    pickle.dump(data, open(fnameSave, 'wb'))
    netCDFFileName = os.path.join(path2GPSfname, f"FRF_geomorphology_elevationTransects_survey_{fnameBase}.nc")
    py2netCDF.makenc_generic(netCDFFileName, "yaml_files/transect_Global.yml", "yaml_files/transect_variables.yml",
                             data)
    print('write mast height used for the day as attribute somewhere')
    print('TarFiles When Done')

crawlerProcessingWorkFlow("/home/crawls/gss_logs/20220316")

if __name__ == '__main__':
    
    crawlerProcessingWorkFlow()