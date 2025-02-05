
# script will read in hypack raw file and key on "MSG" line containing the raw GPS message
# will fill 3 lists with values - time(UTC HHMMss.ss) latitude,longitude,ellipsoid
# edit inputFile with path to your hypack raw file
import pandas as pd


def readHypackRAW(fname):
    # open hypack raw file and key on "msg" line containing raw GPS string
    fc = open(fname)
    
    time_RTN,lat_RTN,long_RTN,ellipsoid_RTN = [],[],[],[]
    
    
    for line in fc:
        if line.startswith("MSG"):
            msg = line.split(',')
            if  msg[8]== '3': # msg[0][4]== '1' and the 3 is indicator that point is RTK fixed
                time_RTN.append(float(msg[2])) # extracting UTC time
                lat_degree_RTN = msg[4][0:2]
                lat_minute_RTN = msg[4][2::]
                lat_decimal_min_RTN = float(lat_minute_RTN)/60 + float(lat_degree_RTN)
                dd_lat_RTN = lat_decimal_min_RTN
                lat_RTN.append(dd_lat_RTN)
                long_degree_RTN = msg[6][0:3]
                long_minute_RTN = msg[6][3::]
                long_decimal_min_RTN = float(long_minute_RTN)/60 + float(long_degree_RTN)
                long_decimal_min = long_decimal_min_RTN *-1
                long_RTN.append(-long_decimal_min_RTN)
                ellipsoid_RTN.append(float(msg[11][3::]))
    fc.close()
    
    # now place into pandas data frame
    df = pd.DataFrame(list(zip(time_RTN, lat_RTN, long_RTN, ellipsoid_RTN)), columns=['time', 'latitude', 'longitude',
                                                                                 'ellipsoid'])
    return df


