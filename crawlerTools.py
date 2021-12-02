import pandas as pd
from pygeodesy import geoids

def correctEllipsoid(fname, geoidFile='g2012bu8.bin', plot=False):
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
    geoidHeight = instance.height(data.latitude, data.longitude)
    data['elevation_NAVD88_m'] = data.ellipsoid - geoidHeight
    if plot is True:
        data.plot.scatter('longitude', 'latitude', c='elevation_NAVD88_m', cmap='ocean')
    
    return data

def cleanDF(data):
    """cleans data frame by removing columns with all zeros and converts time"""
    for key in data.keys():
        if (data[key] == 0).all():
            data.drop(columns=key, inplace=True)
    data['time'] = pd.to_datetime(data['UNIX_timestamp'], unit='s')
    return data