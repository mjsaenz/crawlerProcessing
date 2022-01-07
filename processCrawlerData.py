"""This script was developed by Maile for inital comparison between crawler data and observations.  Using derived
bathy from pressure and RTK comparison against LARC data .... Offset from GPS antenna to ground is unknown at this
time"""
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import warnings
import crawlerFunctions as cf
warnings.filterwarnings("ignore")

file_changes = ['20201117_165351.848','20201117_162048.480','20201117_154249.573','20201117_150844.037']
file_append = 0
fig, axs = plt.subplots(2,sharex=True, sharey=True)

# load tide and larc ground truth
f_tide = 'C:/Users/Maile/Desktop/FRF/NOAATide20201117.csv'
f_larc = 'C:/Users/Maile/Desktop/FRF/xyAdded/LARC_xyAdded2020.csv'
larc_df = pd.read_csv(f_larc, header=0)
larc_df = larc_df.iloc[::-1]
tide = pd.read_csv(f_tide, header=0)

# Initialize RMSE lists
rmse_nav_up = []
rmse_nav_back = []
rmse_rtk = rmse_nav = rmse_130 = rmse_180 = rmse_230 = []

# Larger loop is to go through each csv file
for ind, c in enumerate(file_changes):
    print(ind)
    # Load files and filter/ adjust WL
    f1 = 'C:/Users/Maile/Desktop/FRF/xyAdded/rtk_df_xyAdded' + c + '.csv'
    f2 = 'C:/Users/Maile/Desktop/FRF/xyAdded/nav_df_xyAdded' + c + '.csv'
    rtk_df,nav_df = cf.load_nav_rtk(f1,f2,tide)

    # Find indices of run starts and ends
    nav_splits = cf.split_runs(nav_df)
    rtk_splits = cf.split_runs(rtk_df)

    indstart = 0
    # Inner loop is to split into runs and process each run separately
    for isplit in range(1,len(nav_splits)):
        rtk_fout = 'C:/Users/Maile/Desktop/FRF/Crawler/rtk_df_gridded_' + str(file_append) + '.csv'
        nav_fout = 'C:/Users/Maile/Desktop/FRF/Crawler/nav_df_gridded_' + str(file_append) + '.csv'

        # Get index of run start and end, and create subset of data
        start_i = nav_splits[indstart]
        end_i = nav_splits[isplit]
        nav_x = nav_df.x[start_i:end_i]
        nav_d = nav_df.depth_adj[start_i:end_i]
        # Determine if run is "up" or "back" (eg: away from shore vs towards shore)
        if nav_df.depth_adj[start_i] > nav_df.depth_adj[end_i]: up = 1
        else: up =0
        nav_y = nav_df.y[start_i:end_i]

        start_i = rtk_splits[indstart]
        end_i = rtk_splits[isplit]
        rtk_x = rtk_df.x[start_i:end_i]
        rtk_d = rtk_df.depth[start_i:end_i]
        rtk_y = rtk_df.y[start_i:end_i]
        indstart = isplit

        # Grid the individual run depths to 0.5 m increments in x direction, on specified y-transect
        nav_grid_df = cf.grid_crawler(nav_x,nav_y,nav_d,larc_df,up)
        rtk_grid_df = cf.grid_crawler(rtk_x,rtk_y,rtk_d,larc_df,up,rtk=True)

        # If complete run, calculate RMSE and add to plots
        if (max(nav_grid_df.x) > 230) & (min(nav_grid_df.x) < 180):

            # Save individual runs as file
            #nav_grid_df.to_csv(nav_fout)

            # Calculate RMSEs
            rmse_list_rtk = cf.calculate_rmse(rtk_grid_df)

            '''
            # Uncomment to split plot by NAV out and back
            # Determine which plot to add to (up or back)
            if (up == 1):
                axs[0].plot(nav_grid_df.x,nav_grid_df.depth)
                rmse_nav_up.append(rmse_list[0])
            else:
                axs[0].plot(nav_grid_df.x,nav_grid_df.depth)
                rmse_nav_back.append(rmse_list[0])
            '''

            # Uncomment to split subplot by NAV and RTK
            axs[0].plot(rtk_grid_df.x,rtk_grid_df.depth)
            rmse_rtk.append(rmse_list_rtk[0])
            axs[1].plot(nav_grid_df.x,nav_grid_df.depth)
            rmse_list = cf.calculate_rmse(nav_grid_df)
            rmse_130.append(rmse_list[1])
            rmse_180.append(rmse_list[2])

            
            # Add LARC data to plot at beginning
            if (ind == 3) & (up==1):
                axs[0].plot(nav_grid_df.x,nav_grid_df.larc,'k:',label= 'LARC (ground truth)')
                axs[1].plot(nav_grid_df.x,nav_grid_df.larc,'k:',label= 'LARC (ground truth)')
        else:
            print('File ' + c + ' has an incomplete run')
        
        file_append += 1

# Format and save plot
axs[0].set_title('RTK',style = 'italic')
axs[1].set_title('Nav Soln',style = 'italic')
title_name = 'Crawler Bathy Estimates vs LARC, y = 183'
fig.suptitle(title_name)
plt.xlabel('x (m) in FRF Coords')
axs[0].set(ylabel='Depth (m)')
axs[1].set(ylabel='Depth (m)')
axs[0].legend(loc="upper right")
axs[1].legend(loc="upper right")

# Format Upper Subplot
rmsertk = str(int(np.mean(rmse_rtk)*1000)/1000)
print(rmsertk)
axs[0].vlines(130, -5, 0,colors='y', linestyles='-.', label='130 m')
axs[0].annotate('Overall RMSE = ' +rmsertk+'m',(175,-5),ha='right',fontsize='smaller')
axs[0].vlines(180, -5, 0, colors='y', linestyles='-.', label='180 m')
axs[0].vlines(230, -5, 0, colors='y', linestyles='-.', label='230 m')

# Format Lower Subplot
rmse130 = str(int(np.mean(rmse_130)*1000)/1000)
print(rmse130)
axs[1].vlines(130, -5, 0, colors='y', linestyles='-.', label='130 m')
axs[1].annotate('x<180 RMSE = '+rmse130+'m',(175,-5),ha='right',fontsize='smaller')
rmse180 = str(int(np.mean(rmse_180)*1000)/1000)
print(rmse180)
axs[1].annotate('x<230 m, RMSE = '+rmse180+'m',(225,-5),ha='right',fontsize='smaller')
axs[1].vlines(180, -5, 0, colors='y', linestyles='-.', label='180 m')
axs[1].vlines(230, -5, 0, colors='y', linestyles='-.', label='230 m')

# Save plot to files
outFname = 'C:/Users/Maile/Desktop/CombinedCrawlerBathyPlotRTKNAV.png'
plt.savefig(outFname)
plt.show()