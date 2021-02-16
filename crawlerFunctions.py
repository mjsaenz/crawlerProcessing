import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import math 
from scipy.interpolate import griddata
from scipy import interpolate
from sklearn.metrics import mean_squared_error as rmse
from scipy import signal
import geoprocess as gp

def addFRFxy(lat,lon,depth):
    # convert to FRF coords
    x=y=[]
    for i in range(0,lat):
        coordsOut = gp.FRFcoord(lon[i], lat[i], coordType='LL')
        x.append(np.asarray(coordsOut['xFRF']))
        y.append(np.asarray(coordsOut['yFRF']))
    return x,y

def split_runs(nav_df,start_x=110):
    '''
    Determines the start and end indices in nav_df
    start_x = low estimate of where to start analysis in x direction (approx shoreline position)
    '''
    splits = []
    xi=nav_df.x
    xi[xi<start_x]=0.0
    for i in range(1,len(xi)-10): 
        if ((xi[i] == 0) & ((xi[i+1]>0) or (xi[i-1]>0))): splits.append(i)
    newsplits = []
    istart = 0
    for isplit in splits:
        if (istart < isplit): 
            x_max = xi[istart:isplit].max()
            if (x_max > 200): newsplits.extend([istart, (istart+list(xi[istart:isplit]).index(x_max))])
            istart = isplit

    newsplits.append(splits[(len(splits)-1)])
    return newsplits

def apply_filter_nav(nav_df,T=40,fs=0.04):
    '''
    Butterworth IIR filter to remove waves from Nav Soln
    T = period to decimate below (1/cutoff frequency essentially)
    fs = sampling rate

    Returns filtered signal
    '''
    #make nav length odd for filter implementation
    if (len(nav_df.depth) % 2) == 0: nav_df = nav_df[:-1]

    # Filter NAV w/ LPF, butterworth IIR
    fs = 1/fs # sampling f in Hz
    cutoff = 1/T # Hz, decimate T s period or lower
    nyq = 0.5 * fs # Nyquist f in Hz
    normal_cutoff = cutoff / nyq
    b, a = signal.butter(2, normal_cutoff, btype='low', analog=False)
    nav_filtered_depth = signal.lfilter(b, a, nav_df.depth)
    nav_df.depth = nav_filtered_depth
    return nav_df

def adjust_WL(nav_df,tide):
    '''
    Adjust depths in nav_df to local with averaged WL across the individual run
    tide = dataframe of NOAA tide (water level) from entire day, ideally sampled at 1 hr or less
    '''
    start_t= nav_df.UNIX_timestamp[0]
    end_t = nav_df.UNIX_timestamp[len(nav_df.UNIX_timestamp)-1] # google why you cant index at -1
    tide_run = tide.WL[(tide.unix <= end_t) & (tide.unix >= start_t)]
    avg_WL = sum(tide_run)/len(tide_run)
    nav_df['depth_adj'] = nav_df.depth - avg_WL
    return nav_df

def load_nav_rtk(f1,f2,tide):
    '''
    f1 = RTK file in csv format
    f2 = NAV SOLN file in csv format
    Loads files and filters and adjusts NAV SOLN to NAVD88
    Returns two data frames: RTK and NAV SOLN
    '''
    rtk_df = pd.read_csv(f1, header=0)
    nav_df = pd.read_csv(f2, header=0)
    nav_df.depth = nav_df.depth*(-1)

    nav_df = apply_filter_nav(nav_df)
    nav_df = adjust_WL(nav_df,tide)
    
    return rtk_df,nav_df

def grid_crawler(x,y,z,larc_df,up,x_start=130,y_transect=183,rtk=False):
    '''
    Resamples the depth at every 0.5m in x direction and on y transect
    x = x data
    y = y data
    z = depth
    x_start = beginning of x grid, approximately where shoreline position is
    y_transect = transect of analysis
    rtk = False, whether or not to adjust for RTK offset
    Returns data frame gridded on specified y transect
    '''    
    # Grid data
    grid_df = pd.DataFrame()
    x_grid = np.arange(x_start, max(x), 0.5)
    grid_df['x'] = x_grid
    if up == 0:
        z = z[::-1]
        x = x[::-1]

    # Use same x grid for larc and nav/rtk
    z_grid = np.interp(x_grid,x,z)
    larc_grid = griddata((larc_df.x,larc_df.y), larc_df.depth, (x_grid,183), method='cubic')

    # Adjust rtk to offset value
    if rtk:
        rtkadj = z_grid[0] - larc_grid[0]
        z_grid = z_grid-rtkadj
    grid_df['larc'] = larc_grid
    grid_df['depth'] = z_grid
    
    grid_df.dropna(inplace = True)
    return grid_df

def calculate_rmse(grid_df):
    # Calculate overall RMSEs and RMSE at 50 m increments for NAV data
    larc_grid = grid_df.larc
    d = grid_df.depth
    x = grid_df.x
    overall = rmse(larc_grid, d, squared=False)
    seg_1 = rmse(larc_grid[(x < 180)],d[(x < 180)], squared=False)
    seg_2 = rmse(larc_grid[(x >= 180) & (x < 230)], d[(x >= 180) & (x < 230)], squared=False)
    seg_3 = rmse(larc_grid[(x >= 230)], d[(x >= 230)], squared=False)
    return [overall,seg_1,seg_2,seg_3]