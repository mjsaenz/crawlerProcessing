"""This script will calculate the differece between two points for calibration.  it requires two input files, one for
 RTK and another for the "absolute" position.  Script assumes calibration proceedure to start recording of data, then
 pause to collect numerious points at the beginning, then travel (without RTK gps) then stop, and turn RTK on at end.
 wait long enough for RTK to settle in and collect a significant number of points.  Start and end locations should be
 approximately 2 minutes at each to differentiate from normal "travel" data collection. """
import glob
import os
import sys

import pandas as pd
from matplotlib import pyplot as plt

sys.path.append('/home/c2i/repos')
from testbedutils import geoprocess as gp
import numpy as np


def find_start_end(counts, bins):
    """ finds the two highest histogram values in the counts then averages the bins (counts are between bin values) to
    create "start" and "end" locations.  Lat/Lon was development environment, but should operate on any coordinate
    system.
    Args:
        'counts'(list): a list of counts for each of the bin values
        'bins' (list): a list of bin values.  Should be of length len(counts)+1
    Returns:
        coordinate values (averaged between bin boundaries)
        counts found to identify each coordinate value
    """
    idxs = sorted([(x, i) for (i, x) in enumerate(counts)], reverse=True)[:2]
    vals = ((bins[idxs[0][1]] + bins[idxs[0][1] + 1]) / 2, (bins[idxs[1][1]] + bins[idxs[1][1] + 1]) / 2)
    # print("counts {} for loc {}".format( (idxs[0][0], idxs[1][0]), vals))

    return vals, (idxs[0][0], idxs[1][0])


def get_distance_between_pts(pos1, pos2):
    """Function will take two points and calculate distance over ground between them.

    This process will convert to north Carolina state plane meters to calculatea process.

    Args:
        pos1(tuple): assumed to be longitude/latitude of length two
        pos2(tuple): assumed to be longitude/latitude of length two

    Returns:
        Euclidian distance between two points

    """
    out1 = gp.FRFcoord(pos1[0], pos1[1])
    out2 = gp.FRFcoord(pos2[0], pos2[1])
    return np.sqrt((out1['StateplaneE'] - out2['StateplaneE']) ** 2 + (out1['StateplaneN'] - out2['StateplaneN']) ** 2)


def calibration_comparison(NavSolnFname, RTKfname):
    ####################################################################################
    # fname = "20200826_202745.117_telemetry.gssbin_OPENINS_NAV_SOLUTION.csv"
    # RTKfname = "20200826_202745.117_telemetry.gssbin_GPS_STAT_2.csv"

    # load data into data frames
    df_rtk = pd.read_csv(RTKfname, header=4)
    df_rtk['time'] = pd.to_datetime(df_rtk['UNIX_timestamp'], unit='s')  # convert time to datetimee

    df_enc = pd.read_csv(NavSolnFname, header=4)
    df_enc['time'] = pd.to_datetime(df_enc['UNIX_timestamp'], unit='s')  # convert time to datetime

    ##################################################################################
    binCount = 1000
    ######################################################################################################
    plt.figure(figsize=(12, 8))
    ax21 = plt.subplot2grid((8, 8), (4, 4), colspan=4, rowspan=4)
    enc_lon = ax21.hist(df_enc['absolute_position_0'], bins=binCount, label='Nav Soln')
    enc_lon_position, enc_lon_counts = find_start_end(enc_lon[0], enc_lon[1])
    ax21.plot(enc_lon_position, enc_lon_counts, 'rd')
    rtk_lon = ax21.hist(df_rtk['longitude'], bins=binCount, label='RTK')
    rtk_lon_position, rtk_lon_counts = find_start_end(rtk_lon[0], rtk_lon[1])
    ax21.plot(rtk_lon_position, rtk_lon_counts, 'kx')
    ax21.set_title('longitude counts')
    ax21.legend()
    ax21.semilogy()

    ax22 = plt.subplot2grid((8, 8), (0, 4), colspan=4, rowspan=4)
    enc_lat = ax22.hist(df_enc['absolute_position_1'], bins=binCount, label='Nav Soln')
    enc_lat_position, enc_lat_counts = find_start_end(enc_lat[0], enc_lat[1])
    ax22.plot(enc_lat_position, enc_lat_counts, 'rd')
    rtk_lat = ax22.hist(df_rtk['latitude'], bins=binCount, label='RTK')
    rtk_lat_position, rtk_lat_counts = find_start_end(rtk_lat[0], rtk_lat[1])
    ax22.plot(rtk_lat_position, rtk_lat_counts, 'kx')
    ax22.set_title('latitude counts')
    ax22.legend()
    ax22.semilogy()

    ax33 = plt.subplot2grid((8,8), (0,0), colspan=4, rowspan=8)
    ax33.plot(df_enc['absolute_position_0'], df_enc['absolute_position_1'], '.', label='Nav Soln')
    ax33.plot(df_rtk['longitude'], df_rtk['latitude'], '.', label='RTK')
    ax33.plot(enc_lon_position[0], enc_lat_position[0], 'rd', ms=10, label='Nav Soln')
    ax33.plot(enc_lon_position[1], enc_lat_position[1], 'rd', ms=10)
    ax33.plot(rtk_lon_position[0], rtk_lat_position[0], 'kx', ms=20, label='RTK')
    ax33.plot(rtk_lon_position[1], rtk_lat_position[1], 'kx', ms=20)
    ax33.legend()
    ax33.set_aspect('equal', 'box')

    distStart = get_distance_between_pts(pos1=(rtk_lon_position[0], rtk_lat_position[0]),
                                         pos2=(enc_lon_position[0], enc_lat_position[0]))
    distEnd = get_distance_between_pts(pos1=(rtk_lon_position[1], rtk_lat_position[1]),
                                       pos2=(enc_lon_position[1], enc_lat_position[1]))
    totalDistance = get_distance_between_pts(pos1=(rtk_lon_position[0], rtk_lat_position[0]),
                                             pos2=(rtk_lon_position[1], rtk_lat_position[1]))
    TotalDistanceNav = get_distance_between_pts(pos1=(enc_lon_position[0], enc_lat_position[0]),
                                                pos2=(enc_lon_position[1], enc_lat_position[1]))

    ax33.set_title('Total RTK Distance Traveled = {:.1f}m\n& NAV Distance Traveled = {:.1f}m\n'.format(totalDistance,
                                                                                                     TotalDistanceNav))
    ax33.text((enc_lon_position[0] + rtk_lon_position[0]) / 2, (enc_lat_position[0] + rtk_lat_position[0]) / 2,
              "Start sep\n{:.3f}".format(distStart))
    ax33.text((enc_lon_position[1] + rtk_lon_position[1]) / 2, (enc_lat_position[1] + rtk_lat_position[1]) / 2,
              "Start sep\n{:.3f}".format(distEnd))

    plt.tight_layout()
    outFname = '.'.join(NavSolnFname.split('.')[:2]) + "_CalibrationResults.png"
    plt.savefig(outFname)
    print('Saved Output Here: {}'.format(outFname))
    plt.close()


if __name__ == "__main__":
    args = sys.argv[1:]

    navFlist = sorted(glob.glob(os.path.join(args[0], '*OPENINS_NAV_SOLUTION*.csv')))
    rtkFlist = sorted(glob.glob(os.path.join(args[0], '*GPS_STAT_2*.csv')))
    assert len(navFlist) > 1, "Script didn't find appropriate number of files in the folder"
    assert len(rtkFlist) > 1, "Script didn't find appropriate number of files in the folder"
    # establish how i want to loop (for multiple files)
    if len(args) > 1 and args[1].lower() == 'all':
        plotCount = len(navFlist)

    elif len(args) == 1 or args[1].lower() == 'one':
        # this is for a single plot(the most recent data) -- default behavior
        plotCount = 1

    else:
        print('This script requires a folder with csv files for RTK gps and Open INS nav solutions ')

    # now run loop
    for i in range(int(plotCount)):
        if plotCount == 1 and args[1].lower() == 'one':
            # find the specific file in file Lists to use
            inNav = [i for i in navFlist if args[-1] in i][0]
            inRTK = [i for i in rtkFlist if args[-1] in i][0]
        elif plotCount == 1:
            inNav = navFlist[-1]
            inRTK = rtkFlist[-1]

        else:
            inRTK = rtkFlist[i]
            inNav = navFlist[i]
        # first check to make sure i have proper files
        assert inRTK.find('GPS_STAT_2') > 0, "script didn't find the RTK file"
        assert inNav.find('OPENINS_NAV_SOLUTION') > 0, "script didn't find the proper NAV file"
        assert inNav.split('gss')[0] == inRTK.split('gss')[0], "script isn't using files from the same time"
        # run above calibration script
        calibration_comparison(inNav, inRTK)
