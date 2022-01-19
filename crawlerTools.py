import numpy as np
import pandas as pd
from pygeodesy import geoids
import glob
import os

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
    # first search the path for all files
    flist = glob.glob(os.path.join(os.path.dirname(path2SingleFile), os.path.basename(path2SingleFile).split('.')[0] + "*.csv"))
    # then load NAv solution file
    GPSfname = flist.pop(np.argwhere(["GPS_STAT_2" in f for f in flist]).squeeze())
    data = loadCorrectEllipsoid(GPSfname, geoidFile='data/g2012bu8.bin', plot=False)
    for fname in flist:
        if verbose: print(f'loading {fname}')
        tempdf = pd.read_csv(fname, header=4, error_bad_lines=False)
        if "GPS_STAT_1" in fname:
            tempdf = tempdf.add_prefix("GPS1_")
            tempdf['UNIX_timestamp'] = tempdf.pop("GPS1_UNIX_timestamp")
        dataOut = interpDataFrames(data.UNIX_timestamp, tempdf, verbose=verbose)
        data = data.merge(dataOut, how='left', on="UNIX_timestamp")
        if 'gga_fix_quality' not in data.keys():
            print(f"ERROR HERE: {fname}")
    
    return data



def interpDataFrames(timeStamp2Interp, df, verbose=False): #  = IMUdf , timeStamp2Interp = data['UNIX_timestamp']
    out = pd.DataFrame()
    for key in df.keys():
        if key != 'UNIX_timestamp':
            try:
                out[key] = np.interp(timeStamp2Interp, df.UNIX_timestamp.astype(float), df[key])
            except TypeError:  # typically strings that are all the same
                if verbose: print(f'{key} not included')
                continue
                # if df[key].all(): #[0] == df[key]).all():
                #     out[key] = df[key][:len(timeStamp2Interp)]
            except ValueError:
                if verbose: print(f'{key} not included')
                continue
            finally:
                out['UNIX_timestamp'] = timeStamp2Interp
    return out


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
    Rz= np.matrix([[np.cos(np.deg2rad(90+vYaw)), -np.sin(np.rad2deg(90+vYaw)), 0],
                   [np.sin(np.deg2rad(90+vYaw)),  np.cos(np.rad2deg(90+vYaw)), 0],
                   [                          0,                            0, 1]])
    
    # Rotate GPS Yg 180-Pitch CW so the Xv/Xg and Zv/Zg match
    Ry=np.matrix([[np.cos(np.deg2rad(180-vPitch)),  0, np.sin(np.deg2rad(180-vPitch))],
                  [                             0,  1,                              0],
                  [-np.sin(np.deg2rad(180-vPitch)), 0, np.cos(np.deg2rad(180-vPitch))]])
    
    # Rotate GPS Xg Roll CCW so that Yv/Yg and Yv/Yg match
    # SB question: why neg sin of neg roll???
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

def rotateTranslatePoints(data):
    """takes a crawler data frame and rotates translates all of the points, based on the RotateTranslate function
    Args:
        data: crawler data frame
    
    Returns:
        Assigns new values to a returned data instance with keys
            elevation_NAVD88_m_corrected, yFRF_corrected, xFRF_corrected
         
    """
    for idx in range(data.shape[0]):
        x = data.xFRF.iloc[idx]
        y = data['yFRF'].iloc[idx]
        z = data['elevation_NAVD88_m'].iloc[idx]
        pitch_i = data.attitude_pitch_deg.iloc[idx]
        roll_i = data.attitude_roll_deg.iloc[idx]
        yaw_i =  data.attitude_heading_deg.iloc[idx]
        #translate = np.ones_like(roll) * offset
        newX, newY, newZ = RotateTranslate( vYaw=yaw_i, vPitch=pitch_i, vRoll=roll_i, Gxv=0, Gyv=0, Gzv=-offset,
                                            Gxw=x, Gyw=y, Gzw=z)
        data.at[idx, 'xFRF_corrected'] = newX
        data.at[idx, 'yFRF_corrected'] = newY
        data.at[idx, 'elevation_NAVD88_m_corrected'] = newZ
    
    return data

def TranslateOnly_Wrong(data, verticalOffset, pitch=0, roll=0):
    """Rotates and translates the measured antenna elevations to the ground
    
    Args:
        data: a data frame (generated by loadCorrectElipsoid)
        verticalOffset: total offset between antenna centroid and ground
        pitch: from imu
        roll:  from imu

    Returns:
        corrected vertical elevation
        
    Notes:
        this is just a place holder for the matrix version
        
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