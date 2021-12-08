import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import crawlerFunctions as cf

file_changes = ['20190624_185609.342','20190625_170616.362', '20190625_172934.497']
                # done: '20190626_150153.830','20190626_143524.431', 
                # error: '20190625_201247.208','20190625_192717.445', '20190624_172623.177',
path_in = 'C:/Users/Maile/Desktop/FRF/'
path_out = 'C:/Users/Maile/Desktop/FRF/'
# larc_file_changes = ['FRF_20190625_1165_FRF_NAVD88_LARC_GPS_UTC_v20190627.csv']

for i,c in enumerate(file_changes):
    #load files
    f1 = path_in + c + '_telemetry.gssbin_OPENINS_GPS_STAT.csv'
    f2 = path_in + c + '_telemetry.gssbin_OPENINS_NAV_SOLUTION.csv'
    
    # filenames to save
    f1out = path_out + 'rtk_df_xyAdded' + c + '.csv'
    f2out = path_out + 'nav_df_xyAdded' + c + '.csv'

    rtk_df = pd.read_csv(f1, header=0)
    nav_df = pd.read_csv(f2, header=4)

    # Make sure that the header titles in the csv files match
    x,y = cf.addFRFxy(rtk_df.lat,rtk_df.lon,rtk_df.depth)
    nav_df['x'] = x
    nav_df['y'] = y

    rtk_df.to_csv(f1out)
        

