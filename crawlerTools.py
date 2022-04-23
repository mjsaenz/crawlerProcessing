import heapq
import math
import numpy as np
import pandas as pd
import scipy.signal
from matplotlib import pyplot as plt
from pygeodesy import geoids
import glob
import os
from testbedutils import sblib as sb
from scipy.spatial import cKDTree
from testbedutils import geoprocess as gp
import statistics
from tkinter import *
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
import pandas as pd
import pyproj
import numpy as np
import crawlerFunctions as cf
from datetime import datetime
import utm
from pygeodesy import geoids
import glob
import os
from getDataFRF import getObs
import math
import pickle
import xlsxwriter


def transectSelection(data):
    """
      Args:
          data: dataframe containing crawler transect data, to be modified with isTransect and profileNumber columns

      Returns:
          data: input dataframe with columns isTransect and profileNumber added, isTransect is a boolean denoting
          whether a point is part of a transect, profileNumber is a float to designate the transect a point is a part of
          typically the mean FRFy coordinate of the transect, points not part of a transect are assigned a profileNumber
          of Nan
    """
    # added columns for isTransect boolean and profileNumber float to data dataframe
    data["isTransect"] = [False]*data.shape[0]
    data["profileNumber"] = [float("nan")]*data.shape[0]
    # create copy of data for display, points are removed from the frame once identified as part of a transect
    dispData = data.copy(deep=True)
    # main loop for identifying transects, continues to allow for selections while user inputs y/Y
    transectIdentify = input("Do you want to select a transect? (Y/N):")
    while transectIdentify == "Y" or transectIdentify == "y":

        print("To identify a transect, please place a single point at the start and end of the transect with left click")
        print("Right click to erase the most recently selected point. Middle click (press the scroll wheel) to save.")
        print("Points have saved when they no longer appear on the graph, close the graph window to proceed.")
        print("Remember to remove points used in zooming and panning with right click.")
        print("If more or less than 2 points are selected, no changes will be made")
        print("Each graph's colorscale represents the y axis of the other graph, i.e. the colorscale of the xy graph is time, and vice versa")
        print("Select the transect using only 1 graph at a time")
        # displays plots of two subplots, one with x vs y colored in time and one with time vs y colored in x
        fig, axs = plt.subplots(2)
        fig.suptitle("Transects xFRF (top) and time (bottom) vs yFRF ")
        axs[0].scatter(dispData["xFRF"], dispData["yFRF"], c=dispData["time"], cmap='hsv', s=1)
        axs[0].set(xlabel="FRF Coordinate System X (m)", ylabel = "FRF Coordinate System Y (m)")

        axs[1].scatter(dispData["UNIX_timestamp"], dispData["yFRF"], c=dispData["xFRF"], cmap='hsv', s=1)
        axs[1].set(xlabel="UNIX Timestamp (seconds)", ylabel = "FRF Coordinate System Y (m)")
        nodes = plt.ginput(-1, 0)
        print("Selected Points: ")
        print(nodes)

        # ginput returns list of tuples of selected coordinates, each is in its graph's proper scale
        if len(nodes) == 2:
            # false means ycoord is yFRF, true means UNIX Timestamp
            isTime = [False, False]
            isTime[0] = nodes[0][0] > 1500
            isTime[1] = nodes[1][0] > 1500
            if isTime[0] == isTime[1]:
                endpts = []
                # each node is matched to the closest point in the dispData dataframe
                for x in range(len(nodes)):
                    curr = nodes[x]
                    prevDist = float('inf')
                    closest = tuple()
                    for y in range(dispData.shape[0]):
                        if isTime[x]:
                            dist = math.sqrt(
                                (dispData["UNIX_timestamp"][y] - curr[0]) ** 2 + (dispData["yFRF"][y] - curr[1]) ** 2)
                        else:
                            dist = math.sqrt((dispData["yFRF"][y] - curr[1]) ** 2 + (dispData["xFRF"][y] - curr[0]) ** 2)
                        if dist < prevDist:
                            prevDist = dist
                            closest = (dispData["xFRF"][y], dispData["yFRF"][y])
                    endpts.append(closest)

                # identify endpoints within dispdata frame
                isEndPt = []
                for x in range(dispData.shape[0]):
                    if (dispData["xFRF"][x], dispData["yFRF"][x]) in endpts:
                        isEndPt.append(True)
                    else:
                        isEndPt.append(False)
                dispData["endPt"] = isEndPt
                # endPt column identifies where each transect starts and stops

                # identify transect within dispdata
                isTransect = []
                betweenNodes = False
                for x in range(dispData.shape[0]):
                    if dispData["endPt"][x] and not betweenNodes:
                        # first node in time of transect
                        betweenNodes = True
                        isTransect.append(True)
                    elif dispData["endPt"][x] and betweenNodes:
                        # last node in time of transect
                        betweenNodes = False
                        isTransect.append(True)
                    else:
                        isTransect.append(betweenNodes)
                dispData["isTransect"] = isTransect

                # assign id to current transect
                currTransect = dispData.loc[dispData["isTransect"] == True]
                # remove newly assigned transect from display dataframe
                dispData = dispData.loc[dispData["isTransect"] == False]
                dispData = dispData.reset_index(drop=True)
                meanY = statistics.mean(currTransect["yFRF"])
                print("Close the window to continue.")
                plt.figure()
                plt.hist(currTransect["yFRF"])
                plt.title("FRFy coords of selected transect")
                plt.show()
                print("Mean FRFy coord of selected transect: ", meanY)
                transectID = float(input("What profile number would you like to assign this transect? (float type): "))
                currTransect['profileNumber'] = currTransect['profileNumber'].replace([float("nan")], transectID)

                print("Updating dataframe...")
                # update primary dataframe
                # search once to find first timestamp, iterate afterwards
                startTime = currTransect["UNIX_timestamp"].iloc[0]
                endTime = currTransect["UNIX_timestamp"].iloc[currTransect.shape[0] - 1]
                firstI = 0
                for y in range(data.shape[0]):
                    if data["UNIX_timestamp"].iloc[y] == startTime:
                        firstI = y
                        break

                for x in range(currTransect.shape[0]):
                    data.loc[x + firstI, "profileNumber"] = transectID
                    data.loc[x + firstI, "isTransect"] = True
            else:
                # ignore selected points if from different plots
                print("Selected points from different plots. Discarding selected points.")
        else:
            # ignore selected points if more or less than 2 selected
            print("Selected more or less than 2 points. Discarding selected points.")

        # display selected transects overlayed over all points, colored by profile number
        print("Displaying current progress. Close the window to continue.")
        transectsOnly = data.loc[data["isTransect"] == True]
        plt.figure()
        plt.scatter(data["xFRF"], data["yFRF"], c="black", s=1)
        plt.scatter(transectsOnly["xFRF"], transectsOnly["yFRF"], c=transectsOnly["profileNumber"].to_list(),
                    cmap='hsv', s=1)
        cbar = plt.colorbar()
        plt.xlabel("FRF Coordinate System X (m)")
        plt.ylabel("FRF Coordinate System Y (m)")
        cbar.set_label('Transect Number')
        plt.title("Current Progress")
        plt.show()
        transectIdentify = input("Do you want to select another transect? (Y/N):")

    # prompts for saving charts, excel and pickle
    title = input("What would you like to title the charts?: ")
    filenames = input("What would you like to name the files? (Type null to not save file): ")
    if filenames != "null":
        transectsOnly = data.loc[data["isTransect"] == True]
        print("Close the window to continue.")
        plt.figure()
        plt.scatter(transectsOnly["xFRF"], transectsOnly["yFRF"], c=transectsOnly["profileNumber"].to_list(),
                    cmap='hsv', s=1)
        cbar = plt.colorbar()
        plt.xlabel("FRF Coordinate System X (m)")
        plt.ylabel("FRF Coordinate System Y (m)")
        cbar.set_label('Transect Number')
        plt.title(title)
        plt.savefig(filenames + ".png")
        plt.show()

        print("Close the window to continue.")
        plt.figure()
        plt.scatter(data["xFRF"], data["yFRF"], c="black", s=1)
        plt.scatter(transectsOnly["xFRF"], transectsOnly["yFRF"], c=transectsOnly["profileNumber"].to_list(),
                    cmap='hsv', s=1)
        cbar = plt.colorbar()
        plt.xlabel("FRF Coordinate System X (m)")
        plt.ylabel("FRF Coordinate System Y (m)")
        cbar.set_label('Transect Number')
        plt.title("Identified Transects vs All Points")
        plt.savefig(filenames + "Overlayed.png")
        plt.show()

        data.to_excel(filenames + ".xlsx", engine='xlsxwriter')
        data.to_pickle(filenames + ".pkl")
    return data


def searchPointsInRadius(groundTruth, searchPoints, radius=5, **kwargs):
    """
    
    Args:
        groundTruth: Loop over these points, find closest corresponding point in search point
        searchPoints: search these points for corresponding
        radius: if point located outside of this radius, skip the point

    Keyword Args:
        removeDuplicates: will go through and remove points that are saved to the output dataframe that have same
        groundtruth point
        
    Returns:

    """
    onlylinePoints = kwargs.get('searchOnlyLinePoints', False)
    removeDuplicates = kwargs.get('removeDuplicates', True)
    #intialize output dataframe
    out = pd.DataFrame()
    # make the points to search a KD Tree
    if onlylinePoints:
        
        point_tree = cKDTree(np.array([searchPoints['xFRF'][(searchPoints['profileNumber'] !=0)].to_numpy(), \
                                       searchPoints['yFRF'][(searchPoints['profileNumber'] !=0)].to_numpy()]).T)
    else:
        point_tree = cKDTree(np.array([searchPoints['xFRF'].to_numpy(), searchPoints['yFRF'].to_numpy()]).T)
    seedArray = np.array([groundTruth['xFRF'].to_numpy(), groundTruth['yFRF'].to_numpy()])
    for xy in range(seedArray.shape[1]):
        pointsInRadus = point_tree.data[point_tree.query_ball_point([seedArray[0,xy], seedArray[1, xy]], radius,
                                                                    return_sorted=True)]
        if pointsInRadus.size == 0:
            continue
        else:
            # find the closest point in the list with radius - all directions
            # closestPointIndex = sb.closestRadialNode((seedArray[0, xy], seedArray[1, xy]), pointsInRadus)
            # closet point with corresponding x-location
            closestPointIndex = np.argmin(np.abs(seedArray[0, xy] - pointsInRadus[:,0]))
            
            # find the distance of the point
            dist = math.dist((seedArray[0,xy], seedArray[1, xy]), pointsInRadus[closestPointIndex])
            # now log
            point2Log = searchPoints.query(f'xFRF == {pointsInRadus[closestPointIndex][0]} & yFRF =='
                                f' {pointsInRadus[closestPointIndex][1]}').copy()
            point2Log.loc[point2Log.index,'radius'] = dist
            point2Log.loc[point2Log.index, 'xFRFmatchPoint'] = seedArray[0, xy] # pointsInRadus[closestPointIndex][0]
            point2Log.loc[point2Log.index, 'yFRFmatchPoint'] = seedArray[1, xy] # pointsInRadus[closestPointIndex][1]
            out = out.append(point2Log)
            
            # print(f'closest point to x: {seedArray[0,xy].astype(float):.2f} y: {seedArray[1, xy].astype(float):.2f} is '
            #       f'x: {pointsInRadus[closestPointIndex][0].astype(float):.2f} '
            #       f'y: {pointsInRadus[closestPointIndex][1].astype(float):.2f} '
            #       f'with distance of {math.dist((seedArray[0,xy], seedArray[1, xy]), pointsInRadus[closestPointIndex]):.2f}')
            # print(f'           LOG   x: {point2Log["xFRFmatchPoint"].item():.2f} y: '
            #       f'{point2Log["yFRFmatchPoint"].item():.2f} is x: {point2Log.xFRF.item():.2f} y'
            #       f': {point2Log.yFRF.item():.2f} with distance of'
            #       f' {point2Log["radius"].item():.2f}')
    out.reset_index(inplace=True)
    if removeDuplicates is True:
        for i in range(out.shape[0]): #out.index.uniqe():
            mask = (out['xFRFmatchPoint'].iloc[i] == out['xFRFmatchPoint']) & (out['yFRFmatchPoint'].iloc[i] == out['yFRFmatchPoint'])
            if mask.sum() > 1:
                print(f'duplicate on idx {i}')
                dupIdx = np.argwhere(mask.to_numpy()==True).squeeze()
                deleteIdx = np.delete(dupIdx, out.iloc[dupIdx]['radius'].argmin())
                out.drop(labels=deleteIdx, axis=0)
                out['xFRFmatchPoint'].iloc[i], out['yFRFmatchPoint'].iloc[i]
                
        out.reset_index()
        
    return out
    #
    # for xy in zip(xGround, yGround):
    #     idx = sb.closestRadialNode(xy, data)
    #     spanBtwPoints = math.dist(xy, (bathy['xFRF'].iloc[idx], bathy['yFRF'].iloc[idx]))
    #     if spanBtwPoints <= searchRadius:
    #         print('stuff')
    
def loadAndMergeFiles(path2SingleFile, verbose=True):
    """This function loads all Greensea file types.
    
    function will first go out to find similar csv's with the same time stamp.  then it will correct the GGA geoid
    with NAD83 - geoid 2012B. then it will interpolate all of the files to the GPS time stamp (1HZ), then it will
    compbine all of the data to the NAV Soln data frame
    
    Args:
        path2SingleFile: a single file (full extension)

    Returns:
        single data frame with all csv's combined
        
    """
    print('WARNING: Load and Merge Files has no control of variable names')
    # first search the path for all files
    flist = sorted(glob.glob(os.path.join(os.path.dirname(path2SingleFile), os.path.basename(path2SingleFile).split(
            '.')[0] + "*.csv")))
    if len(flist) == 0: return None
    # then load NAv solution file
    GPSfname = flist.pop(np.argwhere(["GPS_STAT_2" in f for f in flist]).squeeze())
    data = loadCorrectEllipsoid(GPSfname, geoidFile='data/g2012bu8.bin', plot=False)
    for fname in flist:
        if verbose: print(f'loading {fname}')
        if 'SPARTON' in fname or "OPENINS_COMPASS_STAT" in fname:
            if verbose: print(f'skipping {fname}')
            # skipping thes
            continue
        tempdf = pd.read_csv(fname, header=4, error_bad_lines=False)
        if "GPS_STAT_1" in fname:
            tempdf = tempdf.add_prefix("GPS1_")
            tempdf['UNIX_timestamp'] = tempdf.pop("GPS1_UNIX_timestamp")
        
        dataOut = interpDataFrames(data.UNIX_timestamp, tempdf, verbose=verbose)
        data = data.merge(dataOut, how='left', on="UNIX_timestamp")
        if 'gga_fix_quality' not in data.keys():
            print(f"ERROR HERE: {fname}")
    
    return data
def convert2FRF(data):
    coords = gp.FRFcoord(data.longitude.to_numpy(), data.latitude.to_numpy())
    data['xFRF'] = coords['xFRF']
    data['yFRF'] = coords['yFRF']
    return data

def loadAndMergePriorityFiles(path2SingleFile, verbose=True, combineDays=True):
    """This function loads all Greensea file types.
    
    function will first go out to find similar csv's with the same time stamp.  then it will correct the GGA geoid
    with NAD83 - geoid 2012B. then it will interpolate all of the files to the GPS time stamp (1HZ), then it will
    compbine all of the data to the NAV Soln data frame
    
    Args:
        path2SingleFile: a single file (full extension)
        verbose: will print verbose loading information (default=True)
        combineDays: will look at single file name and find the day, then look for other files in the same day to
        merge (default=True)
        
    Returns:
        single data frame with all csv's combined
        
    """
    listofFiles = ["SPARTON", "OPENINS_NAV_SOLUTION"]
    # first search the path for all files
    if combineDays is True:
        dateString = os.path.basename(path2SingleFile).split("_")[0]
        searchString = os.path.basename(path2SingleFile).split("_telemetry")[-1]
        allFilesCollectedToday = glob.glob(os.path.join(os.path.dirname(path2SingleFile), f"*{dateString}*{searchString}"))
    else:
        allFilesCollectedToday = [path2SingleFile]
    
    # run loop for all files in the day
    for path2SingleFile in allFilesCollectedToday:
        # first generate a file list to load all associated files with the GPS 2 file that is priority
        flist = glob.glob(os.path.join(os.path.dirname(path2SingleFile), os.path.basename(path2SingleFile).split(
            '.')[0] + "*GPS*.csv"))
        # could loop through below with different key's for specific files interested in loading
        for key in listofFiles:
            flist.extend(glob.glob(os.path.join(os.path.dirname(path2SingleFile), os.path.basename(path2SingleFile).split(
                '.')[0] + f"*{key}*.csv")))
        if len(flist) == 0: return None
        # then load NAv solution file
        GPSfname = flist.pop(np.argwhere(["GPS_STAT_2" in f for f in flist]).squeeze())
        data = loadCorrectEllipsoid(GPSfname, geoidFile='data/g2012bu8.bin', plot=False)
        # sub-loop for support files for inital file
        for fname in sorted(flist):
            if verbose: print(f'loading {fname}')
            tempdf = pd.read_csv(fname, header=4, error_bad_lines=False)
            
            if "OPENINS_NAV_SOLUTION" in fname: ## this is primary solution
              # attitude_0 ==> attitude_roll_deg; attitude_1 ==> attitude_pitch_deg; attitude_2 ==> attitude_heading_deg
                tempdf=tempdf.add_prefix("NAVSoln_")
                tempdf['UNIX_timestamp'] = tempdf.pop('NAVSoln_UNIX_timestamp')
                tempdf['attitude_roll_deg'] = tempdf['NAVSoln_attitude_0']
                tempdf['attitude_pitch_deg'] = tempdf['NAVSoln_attitude_1']
                tempdf['attitude_heading_deg'] = tempdf['NAVSoln_attitude_2']
                
            elif 'SPARTON_AHRSM2' in fname:
                tempdf = tempdf.add_prefix("SPARTON_")
                tempdf['UNIX_timestamp'] = tempdf.pop("SPARTON_UNIX_timestamp")  # need to rename unix time stamp

            elif "GPS_STAT_1" in fname and "GPS_STAT_2" in path2SingleFile:
                tempdf = tempdf.add_prefix("GPS1_")
                tempdf['UNIX_timestamp'] = tempdf.pop("GPS1_UNIX_timestamp")
            
            dataOut = interpDataFrames(data.UNIX_timestamp, tempdf, verbose=verbose)
            data = data.merge(dataOut, how='left', on="UNIX_timestamp")
            if 'gga_fix_quality' not in data.keys():
                print(f"ERROR HERE: {fname}")
                
            
        if path2SingleFile == allFilesCollectedToday[0]: # if its the first file of the list
            allDataOut = data
        else:
            # allDataOut = allDataOut.append(data, verify_integrity=Flase)
            allDataOut = pd.concat([allDataOut, data], ignore_index=True, axis=0)
        # now clean up, sort and reset the index
        allDataOut.sort_values('UNIX_timestamp', inplace=True, ignore_index=True)
        allDataOut.reset_index(inplace=True, drop=True)
    return allDataOut


def interpDataFrames(timeStamp2Interp, df, verbose=False): #  = IMUdf , timeStamp2Interp = data['UNIX_timestamp']
    out = pd.DataFrame()
    for key in df.keys():
        if key != 'UNIX_timestamp':
            try:
                out[key] = np.interp(timeStamp2Interp, df.UNIX_timestamp.astype(float), df[key])
            except TypeError:  # typically strings that are all the same
                if verbose: print(f'   {key} not included')
                continue
                # if df[key].all(): #[0] == df[key]).all():
                #     out[key] = df[key][:len(timeStamp2Interp)]
            except ValueError:
                if verbose: print(f'   {key} not included')
                continue
            finally:
                out['UNIX_timestamp'] = timeStamp2Interp
    return out

def identfyMastOffset(data, **kwargs):
    """ searches for location of "beanch mark" assumed to be at the FRF, then will isolate that to measure the GPS
    antenna's height from ground.
    
    Args:
        data: pandas data frame of loaded crawler data
        **kwargs:
            "plotting": turns plutton on/off (default=False)

    Returns:

    """
    plotting = kwargs.get('plotting', False)
    coordWindow = 1 # meter in y and x to search for "benchmarkpoint"
    # xWindow = [4, 6] # meters in the x direction
    #yWindow = [515, 517]
    # FRF benchmark in lat/lon
    benchmarkLocation = [-75.751422, 36.182032]
    benchmarkElevation = 5.222 #
    # Identify points that are west of 10 m in FRF x
    coord = gp.FRFcoord(benchmarkLocation[0], benchmarkLocation[1], coordType='LL')
    df_benchMark = data.query(f'xFRF > {coord["xFRF"]-coordWindow} & xFRF < {coord["xFRF"]+coordWindow} '
                              f'& yFRF < {coord["yFRF"]+coordWindow} & yFRF > {coord["yFRF"]-coordWindow}')
    if df_benchMark.empty:
        return data, None
    # now split between the two samples (start/end)
    idxSplit = np.argwhere(np.diff(df_benchMark['UNIX_timestamp']) > 5).squeeze().astype(int) + 1
    start = df_benchMark.iloc[:int(idxSplit)]
    end = df_benchMark.iloc[int(idxSplit):]
    startMed = np.median(start['elevation_NAVD88_m'])
    endMed =  np.median(end['elevation_NAVD88_m'])

    endOffset = endMed - benchmarkElevation
    startOffset = startMed - benchmarkElevation
    offset = np.mean((startOffset, endOffset))
    dataOut = data.query(f'xFRF < {coord["xFRF"]-coordWindow} & xFRF > {coord["xFRF"]+coordWindow} '
               f'& yFRF > {coord["yFRF"]+coordWindow} & yFRF < {coord["yFRF"]-coordWindow}')
    data.drop(df_benchMark.index, axis=0).reset_index(inplace=True)    # greater than 5 seconds
    if plotting == True:
        plt.style.use('seaborn-pastel')
        fig = plt.figure(figsize=(15,5)) #axs = plt.subplots(1,3, figsize=(15,5))
        fig.suptitle(f'todays data {data.time.iloc[0].strftime("%Y-%m-%d")}')
        ax1 = plt.subplot2grid((2,2), (0,0), colspan=2)
        ax2 = plt.subplot2grid((2,2), (1,0))
        ax3 = plt.subplot2grid((2,2), (1,1), sharey=ax2)
        ax1.plot(data['xFRF'], data['yFRF'], '.')
        ax1.plot(df_benchMark['xFRF'], df_benchMark['yFRF'], 'rd', label='identified benchmark')
        ax1.legend()
        ax2.plot(start['time'],  start['elevation_NAVD88_m'] - np.ones_like(start['time'],
                                                                       dtype=float) * benchmarkElevation, 'd')
        ax3.plot(end['time'],  end['elevation_NAVD88_m'] - np.ones_like(endOffset, dtype=float) * benchmarkElevation, 'd')
        ax2.plot(start['time'], np.ones_like(start['time'],dtype=float) * startOffset, label='median offset')
        ax3.plot(end['time'], np.ones_like(end['time'],dtype=float) * endOffset , label='median offset')
        ax2.plot(start['time'], np.ones_like(start['time'], dtype=float) * offset, 'k-', label='Total Offset')
        ax3.plot(end['time'], np.ones_like(end['time'], dtype=float) * offset, '-k', label='Total Offset')
        ax2.set_title('preBeach Runs')
        ax3.set_title('post Beach Runs')
        ax2.set_ylabel('mast offset [m]')
        ax2.legend()
        ax3.legend()
        plt.tight_layout()
    return data, offset

def RotateTranslate(vYaw, vPitch, vRoll,  Gxw, Gyw, Gzw, Gxv=0, Gyv=0, Gzv=0, **kwargs):
    """Rotation and Translation of GPS coordinates from top of crawler mast to vehicle origin (base).  Process
    assumes that data from vehicle/IMU are in North East Down (NED) coordinate system and
    
    Args:
        vYaw: Vehicle heading in global coordinates (degrees)
        vPitch: Vehicle pitch (about y-axis) in local coordinate system (degrees)
        vRoll: Vehicle roll (about x-axis) in local coordinate system (degrees)
        Gxw: X position of the GPS sensor in World coordinates (assumes rectilinear coordinate system, e.g. FRF,
        stateplane or UTM)
        Gyw: Y position of the GPS sensor in World coordinates (assumes rectilinear coordinate system, e.g. FRF,
        stateplane or UTM)
        Gzw: Z position of the GPS sensor in World coordinates - vertical datum
        Gxv: X position in vehicle coordinate system of the GPS sensor (Default=0)
        Gyv: Y position in vehicle coordinate system of the GPS sensor (Default=0)
        Gzv: Z position in vehicle coordinate system of the GPS sensor (keep in mind that positive is down,
        value is likely negative) (Default=0)
        **kwargs:

    Returns:
        newX, newY, newZ corrected global coordinate values Gxw, Gyw, Gzw
    
    Notes:
        kudos to Brittany Bruder for help on this one

    """
    
    verbose=kwargs.get('verbose', True)
    #Rotate GPS 90+Yyaw CW about Zg so Yv and Yw match.
    Rz= np.matrix([[np.cos(np.deg2rad(90+vYaw)), -np.sin(np.deg2rad(90+vYaw)), 0],
                   [np.sin(np.deg2rad(90+vYaw)),  np.cos(np.deg2rad(90+vYaw)), 0],
                   [                          0,                            0, 1]])
    
    # Rotate GPS Yg 180-Pitch CW so the Xv/Xg and Zv/Zg match
    Ry=np.matrix([[np.cos(np.deg2rad(180-vPitch)),  0, np.sin(np.deg2rad(180-vPitch))],
                  [                             0,  1,                              0],
                  [-np.sin(np.deg2rad(180-vPitch)), 0, np.cos(np.deg2rad(180-vPitch))]])
    
    # Rotate GPS Xg Roll CCW so that Yv/Yg and Yv/Yg match
    Rx=np.matrix([[1,                            0,                             0],
                  [0,   np.cos(np.deg2rad(-vRoll)),  -np.sin(np.deg2rad(-vRoll))],
                  [0,   np.sin(np.deg2rad(-vRoll)),  np.cos(np.deg2rad(-vRoll))]])
    
    #Translation Matrices - moving World coords to GPS position
    # so you can rotate around it (world to GPS) -- and back
    # Tw2g= np.matrix([[1, 0, 0, -Gxw],
    #                  [0, 1, 0, -Gyw],
    #                  [0, 0, 1, -Gzw]])
    #
    Tg2w= np.matrix([[1, 0, 0, Gxw],
                     [0, 1, 0, Gyw],
                     [0, 0, 1, Gzw]])
    # create matricies:
    # order of rotation goes right to left
    # Rg2v = Rx @ Ry @ Rz
    # Rv2g = Rg2v.T
    Rv2g = np.matrix(Rx@Ry@Rz).T
    
    # Translation of the vehicle to the GPS
    Tv2g = np.matrix([[1, 0, 0, -Gxv],
                      [0, 1, 0, -Gyv],
                      [0, 0, 1, -Gzv]])
    
    # Vg = Rv2g *np.matrix([Gxv, Gyv, Gzv]).T
    # Tg2v = np.matrix([[1, 0, 0, -Vg.item(0)],
    #                   [0, 1, 0, -Vg.item(1)],
    #                   [0, 0, 1, -Vg.item(2)]])
    # where on vehicle we're interested in projecting the GPS coordinate
    pos = np.matrix([0, 0, 0, 1]).T
    
    sub = Rv2g @ Tv2g @ pos
    Cw = Tg2w @ np.matrix([sub.item(0), sub.item(1), sub.item(2), 1]).T
    
    if verbose: print(f"new Position x: {Cw.item(0):.2f}, y: {Cw.item(1):.2f}, z: {Cw.item(2):.2f}")
    
    return Cw.item(0), Cw.item(1), Cw.item(2)

def rotateTranslatePoints(data, offset, **kwargs):
    """takes a crawler data frame and rotates translates all of the points, based on the RotateTranslate function
    Args:
        data: crawler data frame
        offset: scalar value between GPS antenna and "ground", the code will take care of direction convention
    
    Returns:
        Assigns new values to a returned data instance with keys
            elevation_NAVD88_m_corrected, yFRF_corrected, xFRF_corrected
         
    """
    verbose=kwargs.get('verbose', False)
    data.rename(columns={'xFRF': 'xFRF_GPS', 'yFRF': 'yFRF_GPS', 'elevation_NAVD88_m': 'elevation_NAVD88_m_GPS'},
                inplace=True)
    if ('attitude_pitch_deg' not in data.keys()) & ('attitude_roll_deg' not in data.keys()) & ('attitude_heading_deg'
            not in data.keys()):
        print('  NO GOOD IMU Data to rotate/translate GPS values with')
        return None
    for idx in data.index: #range(data.shape[0]):
        # print(idx)
        x = data['xFRF_GPS'].iloc[idx]
        y = data['yFRF_GPS'].iloc[idx]
        z = data['elevation_NAVD88_m_GPS'].iloc[idx]
        pitch_i = data['attitude_pitch_deg'].iloc[idx]
        roll_i = data.attitude_roll_deg.iloc[idx]
        yaw_i =  data.attitude_heading_deg.iloc[idx]
        #translate = np.ones_like(roll) * offset
        newX, newY, newZ = RotateTranslate( vYaw=yaw_i, vPitch=pitch_i, vRoll=roll_i, Gxv=0, Gyv=0, Gzv=-offset,
                                            Gxw=x, Gyw=y, Gzw=z, verbose=verbose)
        data.at[idx, 'xFRF'] = newX
        data.at[idx, 'yFRF'] = newY
        data.at[idx, 'elevation_NAVD88_m'] = newZ
    
    return data

def TranslateOnly_Wrong(data, verticalOffset):
    """Translates the measured antenna elevations to the ground vertially (not accouting for pitch/roll of mast
    
    Args:
        data: a data frame (generated by loadCorrectElipsoid)
        verticalOffset: total offset between antenna centroid and ground

    Returns:
        corrected vertical elevation
        
    """
    data.elevation_NAVD88_m = data.elevation_NAVD88_m - verticalOffset
    return data

def loadCorrectEllipsoid(fname, geoidFile='g2012bu8.bin', plot=False):
    """ This function loads the csv file that is output by the greensea software and will correct the GGA string
    elevations and geiod to the 2012B geoid that is commonly used at the FRF.  It does this by taking the elevation
    from the GGA string, subtracting the GGA string geoid, leaving you with raw ellipsoid value.  That ellipsoid
    value is then corrected with the 2012B geoid that is commonly used at the FRF.
    
    NOTE: if this is the first time you're using this, you'll likely have to go get the geoid bin file.  Code was
        developed using the uncompressed bin file.  It is unclear if the pygeodesy library requires the bin file to
        be uncompressed.  https://geodesy.noaa.gov/GEOID/GEOID12B/GEOID12B_CONUS.shtml
        
    Args:
        fname: input file to be corrected
        geoidFile: location of your geoidfile (default='g2012bu8.bin')
        plot: show an xy scatter plot with elevations
        
    Returns:
         a pandas object with original and corrected values
         
    """
    
    data = pd.read_csv(fname, header=4)
    data['ellipsoid'] = data.gga_altitude_m + data.gga_height_geoid_m
    instance = geoids.GeoidG2012B(geoidFile)
    if (data.longitude == 0).all() & (data.latitude == 0).all():
        return None
    geoidHeight = instance.height(data.latitude, data.longitude)
    data['elevation_NAVD88_m'] = data.ellipsoid - geoidHeight
    if plot is True:
        data.plot.scatter('longitude', 'latitude', c='elevation_NAVD88_m', cmap='ocean')
    return data

def cleanDF(data, acceptableFix=4):
    """cleans data frame by removing columns with all zeros and converts time, and duplicate frames
    
    """
    if data is None:  # sometimes the data frame loaded is bad
        return data
    data = data.loc[:, ~data.columns.duplicated()]  # remove duplicate columns
    for key in data.keys():
        if (data[key] == 0).all() and key != 'gga_fix_quality':
            data = data.drop(columns=key)
    
    data['time'] = pd.to_datetime(data['UNIX_timestamp'], unit='s')
    data = data.query(f'gga_fix_quality == {acceptableFix}')  # some files oddly have this key
    print(f' FOR PAPER: cleaned acceptable RTK fix with {acceptableFix}')
    data.reset_index(inplace=True, drop=True)  # fixes any index issues
    return data


def closest_node(node, nodes):
    nodes = np.asarray(nodes).T
    try:
        deltas = nodes - node
    except ValueError:
        deltas = nodes.T - node
    dist_2 = np.einsum('ij,ij->i', deltas, deltas)
    return np.argmin(dist_2)


def running_mean(data, window):
    """found running mean function on the internet, untested

    Args:
      data: data to run mean
      window: window over which to take mean

    Returns:
      meaned data

    Reference:
        stack question https://stackoverflow.com/questions/13728392/moving-average-or-running-mean
    """
    return pd.Series(data).rolling(window=window).mean().iloc[window - 1:].values


def closePathOffset(data):
    startX = data['xFRF'][0]
    startY = data['yFRF'][0]
    endX = data['xFRF'][-1]
    endY = data['yFRF'][-1]

    return np.sqrt(np.abs(startX - endX)**2 - np.abs(startY - endY))


def identifyCrawlerProfileLines(data, angleWindow=25, **kwargs):
    """ Function identifies crawler profile lines by taking top two headings (out/back) and then applying a window
    Args:
        data: is the assumed data frame of the crawler processing.  Looks for keys 'attitude_heading_deg' and
        'x/yFRF' for plotting
        angleWindow: Window over which to consider part of the profile line. the difference from the top two headings (default=25)
    Keyword Args:
        "plot": turn plot  on or off (default=True)
        'consecutivePointThresh': connect points that are smaller than this and meet other criteria (default=50).
        default value equates to about a 1.5m
        'lineLengthThreshold': minimum length that a profile line can be in points (defaul=75)
        "lineNumbers": pre=assign line numbers to look at
        "alongshoreSearchWindow": used in combination with line Numbers for finding matching points.
        
    Returns
        data frame
    """
    plotting=kwargs.get('plot', True)
    fname = kwargs.get('fname', 'ProfileLineCrawler.png')
    lineNumbers = kwargs.get("lineNumbers", None)
    lineAngles = kwargs.get('lineAngles', None)
    Thresh4ConcurrentLine = kwargs.get('consecutivePointThresh', 40)   # number of consecutive points required for a line
    lineLengthThreshold = kwargs.get('lineLengthThreshold', 150)
    alongshoreWin = kwargs.get("alongshoreSearchWindow", 5)
    if lineAngles is None:  # we're not given direct input
        counts, bins, _  = plt.hist(data['attitude_heading_deg'], bins=20)
        print('use np.histogram')
        plt.close()
        val1, val2 = heapq.nlargest(2, counts)
        angle1 =  bins[np.argwhere(counts==val1).squeeze()]
        angle2 = bins[np.argwhere(counts==val2).squeeze()]
        assert (angle1 %angle2) > np.diff(bins).mean() and (angle2 %angle1) > np.diff(bins).mean() , \
            "peak angles are next to each other"
    else:
        angle1, angle2= lineAngles
    if lineNumbers is not None:
        totalIdx =  data['yFRF'] == 20000000 # some totally impossible scenario to initalize a df
        
        for line in lineNumbers:
            mask = (data['yFRF'] > (line - alongshoreWin)) & (data['yFRF'] < (line+alongshoreWin))
            maskDir1 = (data['attitude_heading_deg'] > angle1 - angleWindow) & (data['attitude_heading_deg'] < angle1
                                                                               + angleWindow)
            maskDir2 = (data['attitude_heading_deg'] > angle2 - angleWindow) & (data['attitude_heading_deg'] < angle2
                                                                                + angleWindow)
            totalMask = mask & (maskDir1 | maskDir2)
            if plotting: plt.plot(data['xFRF'][totalMask], data['yFRF'][totalMask], '.')
            totalIdx[mask] = line
        
        if plotting:
            plt.plot(data['xFRF'], data['yFRF'], 'k.', ms=1)
            plt.savefig(os.path.dirname(fname) + '/profileLineID.png'); plt.close()
    else:
        backAngle = np.max([angle1, angle2])
        outAngle = np.min([angle1, angle2])
        outIdx = (data['attitude_heading_deg'] <= outAngle + angleWindow) & \
                 (data['attitude_heading_deg'] >= outAngle - angleWindow)
        inIdx = (data['attitude_heading_deg'] <= backAngle + angleWindow)  & \
                (data['attitude_heading_deg'] >= backAngle - angleWindow)
        totalIdx = outIdx | inIdx  # combining
        peakX, _ = scipy.signal.find_peaks(np.diff(totalIdx))
        # find start/end of lines (correct start of line index (+1)
        peakXpos = np.argwhere(np.diff(totalIdx.astype(int)) > 0).squeeze() + 1
        peakXneg = np.argwhere(np.diff(totalIdx.astype(int)) < 0).squeeze()
        assert len(peakXneg) == len(peakXpos), 'Error: found more outs than backs!'
        
        startIdxOfShortWindows = np.argwhere(peakXpos - peakXneg < Thresh4ConcurrentLine).squeeze()
        if np.ndim(startIdxOfShortWindows) == 0: startIdxOfShortWindows = np.expand_dims(startIdxOfShortWindows, axis=0)
        # take the start of each of the windows, then set the boolean locations to True and delete the location for the
        # ones that are below the
        for idx in startIdxOfShortWindows:
            # corrected = slice(peakX[idx], peakX[idx] + np.diff(peakX)[idx] + 1) #
            mySlice = slice(peakXneg[idx], peakXpos[idx])
            print("    IDprofiles: add logic for turning around")
            totalIdx[mySlice] = True  # these are actually continuing a line
    
        newPeakX,  _ = scipy.signal.find_peaks(np.diff(totalIdx)) # np.delete(peakX, startIdxOfShortWindows)  # remove the
        newPeakX = np.append(0, newPeakX)
        # corrected points from the line start array
        for pp, peak in enumerate(newPeakX):
            if pp == len(newPeakX)-1:
                endSeg = len(totalIdx)
            else:
                endSeg = newPeakX[pp+1]
            #if totalIdx[newPeakX[ppprint(endSeg-peak)
            if (endSeg - peak < lineLengthThreshold) & (totalIdx[newPeakX[pp]] == False):
                totalIdx[peak: endSeg+1] = False
        
        # now remove weirdly short profile line segments that are true
        newPeakX,  _ = scipy.signal.find_peaks(np.diff(totalIdx)) # np.delete(peakX, startIdxOfShortWindows)  # remove the
        # QA/QC plot
        # plt.figure();
        # plt.suptitle(f'LINE ID: lineLengthMin: {Thresh4ConcurrentLine};\nconnect segs shorter than {lineLengthThreshold}')
        # plt.subplot(211)
        # plt.scatter(data['xFRF'], data['yFRF'], c=totalIdx, marker='x')
        # # plt.plot(data['xFRF'][peakX], data['yFRF'][peakX], 'rx', label='og peaks', ms=10)
        # #plt.plot(data['xFRF'][peakX[startIdxOfShortWindows]], data['yFRF'][peakX[startIdxOfShortWindows]], 'kx', ms=10,
        # #          label='peaksRemoved')
        # plt.plot(data['xFRF'][newPeakX], data['yFRF'][newPeakX], 'bX', label='newPeaks')
        # plt.legend()
        # #
        # plt.subplot(212)
        # plt.plot(totalIdx)
        # plt.plot(peakX,np.ones_like(peakX), 'rx', label='og peaks')
        # plt.plot(newPeakX, np.ones_like(newPeakX), 'bx', label='newPeaks')
        # plt.legend()

        newPeakX = np.append(0, newPeakX)
        # now find the median y-position and assign that as the line number
        for idx, peak in enumerate(newPeakX):
            if idx == len(newPeakX)-1:
                mySlice = slice(newPeakX[idx], len(data))
            else:
                mySlice = slice(newPeakX[idx]+1, newPeakX[idx+1])
            if data['profileNumber'][mySlice].all() == True:
                data.loc[mySlice, 'profileNumber'] = np.median(data['yFRF'][mySlice]).astype(int)


    data['profileNumber'] = totalIdx
    
    if plotting is True:
        data['profileNumber'].unique()
    
        plt.figure()
        # plt.scatter(data['xFRF'], data['yFRF'], c=totalIdx)
        plt.scatter(np.ma.array(data['xFRF'], mask=data['profileNumber'] == False), np.ma.array(data['yFRF'],
                    mask=data['profileNumber'] == False), c=np.ma.array(data['profileNumber'], mask=data[
                    'profileNumber'] == False))
        plt.colorbar()
        plt.title(f'Identified profiles lines from crawler (yellow) \n{data["time"][0].strftime("%Y-%m-%d")}')
        plt.savefig(fname)
        plt.close()
    return data


def calculate3Dspeed(t, x, y, z, **kwargs):
    """ calculates speed in three dimensions with input positional timeseries.  if you want to calculate in 2D,
    send zeros/ones for the positonal argument you're not interested in.
    
    Args:
        t: time
        x: x position
        y: y position
        z: z position
        **kwargs:
            "velocities": if True, will also return velocities in x,y,z order (default=False)

    Returns:
        a list of speeds
        if "velocities" = True, then speed, xVelocities, yVelocities, zVelocities

    """
    assert len(x) == len(y) == len(z) == len(t), "t, x,y,z need to be the same length"
    velocities = kwargs.get('velocities', False)
    xVelo, yVelo, zVelo, speed = [0], [0], [0], [0]
    for i in range(len(x)-1):
        if t[i+1] == t[i]:
            speed.append(0)
            xVelo.append(0)
            yVelo.append(0)
            zVelo.append(0)
        else:
            yVelo.append(y[i+1] - y[i])
            xVelo.append(x[i+1] - x[i])
            zVelo.append(z[i+1] - z[i])
            speed.append(np.sqrt( xVelo[-1]**2 + yVelo[-1]**2 + zVelo[-1]**2)/(t[i+1]-t[
                i]).total_seconds())
    if velocities is True:
        return speed, xVelo, yVelo, zVelo
    else:
        return speed