import pandas as pd
import glob, os, sys
from matplotlib import pyplot as plt
import numpy as np

#######################################################################################
# load files
timeString = "20200903_193927.463"
f1 = "/home/c2i/gss_logs/{}_telemetry.gssbin_GPS_STAT_2.csv".format(timeString)
f2 = "/home/c2i/gss_logs/{}_telemetry.gssbin_OPENINS_NAV_SOLUTION.csv".format(timeString)
f3 = "/home/c2i/gss_logs/{}_telemetry.gssbin_C2I_SPC_PCOMMS_STAT.csv".format(timeString)
f4 = "/home/c2i/gss_logs/{}_telemetry.gssbin_OPENINS_IMU_STAT.csv".format(timeString)
## Crow island data
f1 = "/home/c2i/repos/dataProcess/CrowIslandLongRuns/20200812_153930.704_telemetry.gssbin_GPS_STAT_2.csv"
f2 = "/home/c2i/repos/dataProcess/CrowIslandLongRuns/20200812_153930.704_telemetry.gssbin_OPENINS_NAV_SOLUTION.csv"
f3 = "/home/c2i/repos/dataProcess/CrowIslandLongRuns/20200812_153930.704_telemetry.gssbin_C2I_SPC_PCOMMS_STAT.csv"
f4 = ""
#######################################################################################
# plot scatter colored by Z elevations
rtk_df = pd.read_csv(f1, header=4)
# 'UNIX_timestamp', 'LCM_event_timestamp', 'count_publish', 'fix_time', 'fix_valid', 'gga_altitude_m',
# 'gga_fix_quality', 'gga_height_geoid_m', 'gga_horizontal_dop', 'gga_num_satellites_tracked', 'gga_valid', 'gll_fix_quality',
# 'gll_valid', 'gsa_3D_fix', 'gsa_horizontal_dop', 'gsa_num_satellite_prns', 'gsa_positional_dop', 'gsa_valid',
# 'gsa_vertical_dop', 'gsv_derived_edop', 'gsv_derived_hdop', 'gsv_derived_ndop', 'gsv_derived_pdop', 'gsv_derived_vdop',
# 'gsv_gps_downing_error', 'gsv_gps_easting_error', 'gsv_gps_northing_error', 'gsv_measured_ure',
# 'gsv_num_satellites_in_view', 'gsv_num_tracked', 'gsv_valid', 'latitude', 'longitude', 'rmc_course_over_ground_deg',
# 'rmc_magnetic_declination', 'rmc_speed_over_ground_knots', 'rmc_valid', 'unix_time'
nav_df = pd.read_csv(f2, header=4)
# ['UNIX_timestamp', 'LCM_event_timestamp', 'absolute_position_0', 'absolute_position_1', 'absolute_position_2',
# 'absolute_position_ok', 'altitude_above_bottom', 'altitude_ok', 'attitude_0', 'attitude_1', 'attitude_2',
# 'attitude_acceleration_0', 'attitude_acceleration_1', 'attitude_acceleration_2', 'attitude_dot_0', 'attitude_dot_1',
# 'attitude_dot_2', 'attitude_ok', 'count_publish', 'course_over_ground', 'depth', 'depth_ok', 'dvl_error',
# 'have_bottom_lock', 'imu_error', 'initial_lonlat_fix_0', 'initial_lonlat_fix_1', 'last_lonlat_fix_0',
# 'last_lonlat_fix_1', 'relative_acceleration_0', 'relative_acceleration_1', 'relative_acceleration_2',
# 'relative_position_0', 'relative_position_1', 'relative_position_2', 'relative_position_dot_0',
# 'relative_position_dot_1', 'relative_position_dot_2', 'relative_position_ok', 'speed_ok', 'speed_over_ground',
# 'unix_time', 'vehicle']
incl_df = pd.read_csv(f3, header=4)
# ['UNIX_timestamp', 'LCM_event_timestamp', 'SPC_ERR', 'SPC_INCLO_X' 'SPC_INCLO_Y', 'SPC_PRES', 'count_publish', 'sender_id'],
imu_df = pd.read_csv(f4, header=4)
#['UNIX_timestamp', 'LCM_event_timestamp', 'acceleration_valid', 'acceleration_x_m_sec2', 'acceleration_y_m_sec2',
# 'acceleration_z_m_sec2', 'attitude_heading_deg', 'attitude_pitch_deg', 'attitude_roll_deg', 'attitude_valid',
# 'count_publish', 'count_publish_error', 'id_device', 'magnetometer_x', 'magnetometer_y', 'magnetometer_z',
# 'num_analogs', 'rotational_vel_heading_deg_sec', 'rotational_vel_pitch_deg_sec', 'rotational_vel_roll_deg_sec',
# 'rotational_vel_valid', 'status', 'temperature']

#######################################################################################
# break up variables
tDz = incl_df.UNIX_timestamp         # inclenometer timestamps
incl = incl_df.SPC


_INCLO_X           # inclenometer values
tDz = tDz[~np.isnan(incl)]           # remove inclenometer values that are nans
incl = incl[~np.isnan(incl)]         # remove inclenometer values that are nans
tRawSpeed = nav_df.UNIX_timestamp    # nav time stamp
rawSpeed = nav_df.speed_over_ground  # nav speed
inclOffset = -1.8340
# time average speed to interval of inclenometer
avgspd = []
for ii, dz in enumerate(incl):
    try:
        idx = (tRawSpeed < tDz.values[ii + 1]) & (tRawSpeed >= tDz.values[ii])
    except (IndexError, KeyError):
        idx = tRawSpeed >= tDz.values[ii]
    if idx.any():
        avgspd.append(rawSpeed[idx].mean())  # no nans

# integrate measured slope to elevation
#     dt * speed * slope
#     use average speed, measured inclenometer values, and timestamp from inclenometer
elevation, rise, dt = [0], [0], [0]
for ii, tt in enumerate(tDz):
    if ii < tDz.size-1:
        dt.append(tDz.values[ii+1] - tDz.values[ii])
        rise.append((np.sin(np.deg2rad(incl.values[ii]-inclOffset)))*avgspd[ii]*dt[-1])
        elevation.append(elevation[-1]+rise[-1])
#######################################################################################
plt.figure()
# plt.plot(tDz.values, np.array(rise)*-1, '.', ms = 2, label='rise/timestep [m]')
plt.plot(tDz.values, np.array(elevation)+rtk_df['gga_altitude_m'][0], '.', label='calculated elevation [m]\noffset = {}'.format(inclOffset))
# plt.plot(tDz.values, avgspd, '.', label='vehicle speed (from Nav Soln)[m/s]')
# plt.plot(tDz.values, incl.values, '--.', ms=2, label='Raw Incl Vals [rise/run]')
# plt.plot(tDz.values, dt, '.', label='dt')
plt.plot(rtk_df['UNIX_timestamp'], rtk_df['gga_altitude_m'], '.', label='gga_altitude [m]')
plt.legend()
plt.title('Comparison of calculations')

# plt.figure()
# plt.scatter(rtk_df.longitude, rtk_df.latitude, c=rtk_df.gga_altitude_m, vmin=-0.1, vmax=0.5);
# cbar = plt.colorbar();
# cbar.set_label('RTK')
# plt.scatter(nav_df.absolute_position_0, nav_df.absolute_position_1, c=nav_df.absolute_position_2, vmin=-0.1, vmax=0.5);
# cbar = plt.colorbar();
# cbar.set_label('Nav Soln')
#
