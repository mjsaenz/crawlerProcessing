""" This code analyses a single profile repeatedly.  this is for the paper """
import math
import matplotlib
matplotlib.use('TkAgg')
import sys, os
import matplotlib.pyplot as plt
import numpy as np
import crawlerPlots
from testbedutils import geoprocess as gp
import crawlerTools
import datetime as DT
import pickle
import glob
###############################
import py2netCDF
#
# dateString = '20220113' #'20211020'
# savePath = "plots/DUNEXcomparisons"
# GPSfname = "/data/20220322/20220322_192544.302_telemetry.gssbin_GPS_STAT_2.csv"

# path2GPSfname =
def crawlerProcessingWorkFlow(path2GPSfname, **kwargs):
    """
    
    Args:
        path2GPSfname:  Folder path to look for bin files
        **kwargs:
            "binProcessor": path to bin executable (default="/usr/local/bin/gsbin-log-processor")
    
    Returns:

    """
    
    gsBinProcessorLocation = kwargs.get('binProcessor', "/usr/local/bin/gsbin-log-processor")
    flist = glob.glob(os.path.join(path2GPSfname , "*.gssbin"))
    # find all of the raw GSS bin files then unpack them
    for fname in flist:
        os.system(f"{gsBinProcessorLocation} --logfile {fname}")

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

    start = data['time'].iloc[0].to_pydatetime()
    print(f"\n\n Working on date {start}")
    # # now rotate translate for orientation
    data = crawlerTools.rotateTranslatePoints(data, offset)
    data = crawlerTools.transectSelection(data)
    print('  TODO: look for pickle file, if exists, open check for line numbers')
    print ('  TODO Save pickle file after processing')
    
    data['time'] = [DT.datetime.strptime(data['time'][i], "%Y-%m-%d %H:%M:%S") for i in range(len(data['time']))]
    
    # quick comparison plots for QA/QC
    crawlerPlots.bathyEnvalopeComparison(GPSfname, data, bathy=None)
    fname = os.path.join(os.path.dirname(GPSfname), ''.join(os.path.basename(GPSfname).split('.')[0])+"_withLocalObs_XY.png")
    crawlerPlots.bathyPlanViewComparison(fname, data,  bathy=None, topo=None, plotShow=False)
    ############################################################
    print(' TODO: make netCDF files here!')
    fnameSave = 'data/repeatabilityTest2DataCrawler.pickle'
    pickle.dump(data,open(fnameSave, 'wb'))
    netCDFFileName = os.path.join(path2GPSfname, f"FRF_geomorphology_elevationTransects_survey_{dateString}.nc")
    py2netCDF.makenc_generic(netCDFFileName, "yaml_files/transect_Global.yml", "yaml_files/transect_variables.yml",
                             data)
    print('write mast height used for the day as attribute somewhere')
    print('TarFiles When Done')

crawlerProcessingWorkFlow("/home/crawls/gss_logs/20220316")

if __name__ == '__main__':
    
    crawlerProcessingWorkFlow()