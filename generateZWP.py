#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 15:09:02 2025

@author: ahyder
"""

import planetmapper
from planetmapper.kernel_downloader import download_urls
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
from astropy.io import fits
import time
from scipy.signal import correlate
import scipy as sp
import pytz
import datetime
import csv
import os

from zonalwindfunctions import ZWP_Class

long_res, lat_res = 360, 180
long_array, lat_array = np.linspace(360, 0,long_res), np.linspace(-90, 90,lat_res)
R_J = 71492e3 # m


main_path = 'spex/'
file_paths_list_51 = [] # with threshold set to 51 min (0.85)
t_diff_list_51 = []

with open(main_path+"pairs_2024_51min.csv", "r") as f:
    reader = csv.reader(f)
    for row in reader:
        selected_pair, t_diff_of_pair = [row[0], row[1]], float(row[2])
        t_diff_list_51.append(t_diff_of_pair)
        file_paths_list_51.append([selected_pair[0], selected_pair[1]])
        
for i in range(len(file_paths_list_51)):
    print(file_paths_list_51[i], t_diff_list_51[i]*3600)


r0_array = np.array([170.8, 162.7, 162.7, 162.7, 162.7, 162.7, 162.7, 162.7, 198.25, 198.25])

# full path to file
sample_file_ind = 0
file_list = [file_paths_list_51[sample_file_ind][0], file_paths_list_51[sample_file_ind][1]]


#%%

# needed only to populate date_meta_list

# Number of images
num_images = 2
ncols = math.ceil(math.sqrt(num_images))
nrows = math.ceil(num_images / ncols)

fig1, axs1 = plt.subplots(
nrows=nrows, ncols=ncols, dpi=150)
axs1 = axs1.flatten()

date_meta_list = []

for ax1, f in zip(axs1, file_list):
    print(f, type(f))
    
    hdul = fits.open(f)
    hdr = hdul[0].header
    print(hdr['DATE_OBS'], hdr['TIME_OBS'])
    date_tag = hdr['DATE_OBS']
    time_tag = hdr['TIME_OBS']
    utc_str = date_tag + 'T' + time_tag

    date_meta_list.append(utc_str)

    disc_x, disc_y = hdr['CX'], hdr['CY'] # obtained automatically so the disc will be around the correct point
    print(disc_x, type(disc_x))
    print()


    # navigational data made using the GUI - doing it here results in the NASA IRTF being set as observer.
    # Any spice available for it though? Would make things easier!
    observation = planetmapper.Observation(f, 
                                           target='jupiter', 
                                           utc=utc_str, 
                                           observer='EARTH')

    # Run the GUI to fit the disc interactively
    #observation.run_gui() # good for fitting the disc directly
    disc_r0 = r0_array[sample_file_ind]
    # disc_r0 = 198.25


    ax1.imshow(observation.data[0, ...], origin='lower', cmap='hot')
    observation.set_disc_params(x0=disc_x, y0=disc_y, r0=disc_r0)
    observation.plot_wireframe_xy(ax1, formatting={
                                      'grid': {'linestyle': '-', 'linewidth': 0.5, 'alpha': 0.3, 'color':'w'},
                                      'prime_meridian': {'linewidth': 1, 'color': 'r'},
                                      'terminator': {'linewidth': 0.5, 'color':'w'},
                                      'equator': {'linewidth': 2, 'color':'w'},
                                      'limb_illuminated': {'linewidth': 0.5, 'color':'w'},},
                                    indicate_equator=True,
                                )
    ax1.set_title(utc_str)
    # plt.show()
fig1.tight_layout()


plt.show()

#%%

# needed only to populate cmap_meta_list and obs_meta_list

fig3, axs3 = plt.subplots(
nrows=2, ncols=1, dpi=150)
axs3 = axs3.flatten()

cmap_meta_list = []
obs_meta_list = []

'''
This block only works with non-interactive plotting for some stupid reason! It has to do with backend Qt5Agg and the plt.show() command.
It's possible that the observation object messes things up due to how the internal plotting of planetmapper is defined.
Just leave it alone and keep it inline. Use contourf in a different cell if needed ffs...
'''
# keep this here so this cells work well - otherwise planetmapper has a seizure. # KEEP THIS LINE
#%matplotlib inline 

for ax3, f in zip(axs3, file_list):
    
    hdul = fits.open(f)
    hdr = hdul[0].header
    print(hdr['DATE_OBS'], hdr['TIME_OBS'])
    date_tag = hdr['DATE_OBS']
    time_tag = hdr['TIME_OBS']
    utc_str = date_tag + 'T' + time_tag
    disc_x, disc_y = hdr['CX'], hdr['CY']
    
    observation = planetmapper.Observation(f, 
                                            target='jupiter', 
                                            utc=utc_str, 
                                            observer='EARTH')
    observation.set_disc_params(x0=disc_x, y0=disc_y, r0=disc_r0)   
    obs_meta_list.append(observation)

    mapped_data = observation.get_mapped_data('linear')  # get the mapped data
    cmap_meta_list.append(mapped_data)
    
    # Display mapped image and add a useful annotation
    observation.imshow_map(mapped_data[0, ...], ax3, cmap='hot')
    # ax3.contourf(np.linspace(360, 0,360), np.linspace(-90, 90,180), mapped_data[0, ...], 256, cmap='hot')
    # ax3.set_title(hdr['TIME_OBS']+'UTC')
    # ax3.set_ylim(-60, 60)
    ax3.set_title(utc_str)

fig3.tight_layout()

plt.show()

#%%

for f in ['spex/2024feb5/jc118152.cmap.fits.gz']:
    print(f, type(f))
    
    hdul = fits.open(f)
    hdr = hdul[0].header
    print(hdr['DATE_OBS'], hdr['TIME_OBS'])
    date_tag = hdr['DATE_OBS']
    time_tag = hdr['TIME_OBS']
    utc_str = date_tag + 'T' + time_tag
print(hdul[0].data, hdul[0].data.shape)

# Boolean Intersection Map (BIM)
test_bool_intersect = ZWP_Class.getValidLatLonRange(cmap_meta_list[0][0, ...], cmap_meta_list[1][0, ...])


# # Boolean Intersection Map (BIM)
# test_bool_intersect = ZWP_Class.getValidLatLonRange(cmap_meta_list[0][0, ...], cmap_meta_list[1][0, ...])
test_bool_cmap0 = ZWP_Class.getValidLatLonRange(cmap_meta_list[0][0, ...], np.ones((lat_res,long_res)))
test_bool_cmap1 = ZWP_Class.getValidLatLonRange(cmap_meta_list[1][0, ...], np.ones((lat_res,long_res)))

plt.figure(figsize=(6,2), dpi=150)
plt.contourf(long_array, lat_array, test_bool_intersect, 1, cmap='gray', alpha=0.5)
plt.gca().invert_xaxis()
plt.grid(True, 'both', ls=':', alpha=0.5)
# plt.colorbar()
plt.xlabel('Planetographic longitude (W)'), plt.ylabel('Planetographic latitude')
plt.title('Valid coordinates for overlapped regions in the two images - BIM')
plt.show()

print((np.sum(test_bool_intersect) / test_bool_intersect.size) * 100, ' % of the map is valid for use in a ZWP generation')
print('Original map percentages: (Makes sense as it is approximately half of the planet)')
print((np.sum(test_bool_cmap0) / test_bool_cmap0.size) * 100, ' %')
print((np.sum(test_bool_cmap1) / test_bool_cmap1.size) * 100, ' %')



#%%


orig_ex_long = 50
test_shift_long = 10

test_long_shift_1, test_map_shift_1 = ZWP_Class.shiftMap(long_array, cmap_meta_list[0][0, ...], test_shift_long)
test_long_shift_2, test_map_shift_2 = ZWP_Class.shiftMap(long_array, cmap_meta_list[1][0, ...], test_shift_long)
_, intersect_shift_map = ZWP_Class.shiftMap(long_array, test_bool_intersect, test_shift_long)


plt.figure(figsize=(7,3), dpi=150)
plt.contourf(long_array, lat_array, cmap_meta_list[0][0, ...], 256, cmap='gray', alpha=0.5)
plt.contour(long_array, lat_array, test_map_shift_1, 256, cmap='hot', alpha=0.5)
plt.vlines(orig_ex_long, min(lat_array), max(lat_array), color='c', ls='dashed', lw=0.75, label='Sample Longitude')
plt.vlines(orig_ex_long + test_shift_long, min(lat_array), max(lat_array), color='c', lw=0.75, label='Shift by given amount')
plt.gca().invert_xaxis()
plt.grid(True, 'both', ls=':', alpha=0.5)
plt.title('Testing shifting procedure relative to original Western coordinate by ' + str(test_shift_long) + '˚ (22nd)')
plt.xlabel('Planetographic longitude (W)'), plt.ylabel('Planetographic latitude')
plt.legend(loc='lower left')

plt.figure(figsize=(7,3), dpi=150)
plt.contourf(long_array, lat_array, cmap_meta_list[1][0, ...], 256, cmap='gray', alpha=0.5)
plt.contour(long_array, lat_array, test_map_shift_2, 256, cmap='hot', alpha=0.5)
plt.vlines(orig_ex_long, min(lat_array), max(lat_array), color='c', ls='dashed', lw=0.75, label='Sample Longitude')
plt.vlines(orig_ex_long + test_shift_long, min(lat_array), max(lat_array), color='c', lw=0.75, label='Shift by given amount')
plt.gca().invert_xaxis()
plt.grid(True, 'both', ls=':', alpha=0.5)
plt.title('Testing shifting procedure relative to original Western coordinate by ' + str(test_shift_long) + '˚ (22nd)')
plt.xlabel('Planetographic longitude (W)'), plt.ylabel('Planetographic latitude')
plt.legend(loc='lower left')
# plt.show()

plt.figure(figsize=(7,3), dpi=150)
plt.contourf(long_array, lat_array, test_bool_intersect, 1, cmap='gray', alpha=0.75)
plt.contour(long_array, lat_array, intersect_shift_map, 1, cmap='hot', alpha=0.75)
plt.vlines(orig_ex_long, min(lat_array), max(lat_array), color='c', ls='dashed', lw=0.75, label='Sample Longitude')
plt.vlines(orig_ex_long + test_shift_long, min(lat_array), max(lat_array), color='c', lw=0.75, label='Shift by given amount')
plt.gca().invert_xaxis()
plt.grid(True, 'both', ls=':', alpha=0.5)
plt.title('Testing shifting procedure relative to original Western coordinate by ' + str(test_shift_long) + '˚ (22nd)')
plt.xlabel('Planetographic longitude (W)'), plt.ylabel('Planetographic latitude')
plt.legend(loc='lower left')
plt.show()


#%%

'''
Pre- and post margin application longitudes:
'''

boolean_mask = ZWP_Class.defineLongThresh(test_bool_intersect, long_array, lat_array)  # longitudes from intersection
boolean_mask_after_margin = ZWP_Class.applyMargin(boolean_mask, long_array, margin=5)  # (1, 360) boolean
print(boolean_mask)
# print()

# pre-margin
longitude_masked = np.ma.masked_array(long_array, np.logical_not(boolean_mask))
print(longitude_masked) 
print()

# post-margin
longitude_masked_margin = np.ma.masked_array(long_array, np.logical_not(boolean_mask_after_margin))
print(longitude_masked_margin)
print(bool(boolean_mask_after_margin[0]))
# print(shiftMap(long_array, boolean_mask_after_margin, 0))


#%%


plt.figure(figsize=(7,3), dpi=150)
plt.contourf(long_array, lat_array, test_bool_intersect, 1, cmap='gray', alpha=0.5)
plt.hlines(60, min(long_array), max(long_array), color='k', ls='dashed', lw=0.75)
plt.hlines(-60, min(long_array), max(long_array), color='k', ls='dashed', lw=0.75)

pre_margin_indices  = np.where(abs(np.diff(boolean_mask)) == 1)[0]
post_margin_indices = np.where(abs(np.diff(boolean_mask_after_margin)) == 1)[0]

plt.vlines(long_array[pre_margin_indices[0]], min(lat_array), max(lat_array), color='k', ls='dashed', lw=0.75, label='pre-margin')
plt.vlines(long_array[pre_margin_indices[1]], min(lat_array), max(lat_array), color='k', ls='dashed', lw=0.75)

plt.vlines(long_array[post_margin_indices[0]], min(lat_array), max(lat_array), color='c', ls='dashed', lw=0.75, label='post-margin')
plt.vlines(long_array[post_margin_indices[1]], min(lat_array), max(lat_array), color='c', ls='dashed', lw=0.75)

plt.gca().invert_xaxis()
plt.legend(loc='upper right')
plt.grid(True, 'both', ls=':', alpha=0.5)
# plt.colorbar()
plt.xlabel('Planetographic longitude (W)'), plt.ylabel('Planetographic latitude')
plt.title('BIM')
plt.tight_layout()
plt.show()


#%%



critical_lat = 50
shifted_f1, shifted_f2, shifted_BIM, margin_applied_in_func_shifted, final_longitudes_masked_shifted, corr_1D = ZWP_Class.getFullLatitudinalCorrelation(
                            cmap_meta_list[0][0, ...], cmap_meta_list[1][0, ...], 
                            long_array, lat_array, shift_long=90, crit_lat=critical_lat, 
                            margin=3, max_delta_long=8, force_shift=False, lat_res_tol=1e-2)

#%%

delta_t = ZWP_Class.getTemporalDifference(file_list)

zwp_5microns = -ZWP_Class.getSpeedatLat(corr_1D, delta_t) # minus is necessary to match the other observations of the wind speed.
smooth_zwp_5microns = ZWP_Class.getSmoothedProfile(zwp_5microns, 7, 2)

hst_wind_1D = np.genfromtxt('Simon+2015_ZWP_HST.csv', delimiter=",", dtype=None, encoding=None, names=True)

y_lim_range = 50 #critical_lat - critical_lat*0.05 # -5% for aesthetics...
plt.figure(figsize=(3,6), dpi=150)
# plt.plot(corr_1D / 20.5219444, lat_array)
# plt.plot(corr_1D * 71492e3 / 73879, lat_array)
plt.plot(-corr_1D, lat_array, lw=0.5)
plt.xlabel('∆longitude')
plt.ylabel('Planetographic Latitude')
# plt.xlim(-1000, 1000)
# plt.ylim(-y_lim_range, y_lim_range)
plt.show()
print(delta_t)


#%%

fig, axs = plt.subplots(
nrows=1, ncols=3, figsize=(8, 2), dpi=150)
axs[0].contourf(long_array[margin_applied_in_func_shifted.astype(bool)], lat_array, shifted_f1[:, margin_applied_in_func_shifted.astype(bool)], 256, alpha=1)
# axs[0].contour(long_array[margin_applied_in_func_shifted.astype(bool)], lat_array, shifted_f2[:, margin_applied_in_func_shifted.astype(bool)], 128, cmap='bone', alpha=0.25)
axs[1].contourf(long_array[margin_applied_in_func_shifted.astype(bool)], lat_array, shifted_f2[:, margin_applied_in_func_shifted.astype(bool)], 256)
axs[0].grid(True, 'both', ls=':', alpha=0.7), axs[1].grid(True, 'both', ls=':', alpha=0.7)
axs[0].set_ylim(-y_lim_range, y_lim_range), axs[1].set_ylim(-y_lim_range, y_lim_range)
axs[0].invert_xaxis(), axs[1].invert_xaxis()

axs[2].plot(hst_wind_1D['x'], hst_wind_1D['y'], color='b', lw=0.5, label='HST 2015')
axs[2].plot(zwp_5microns, lat_array, marker='s', ms=2, color='k', lw=0.5, label='IRTF 2024')
axs[2].plot(smooth_zwp_5microns, lat_array, color='r', lw=0.9, label='Smoothed')
axs[2].set_xlabel('Zonal velocity, m/s')
# plt.xlim(-1000, 1000)
axs[2].set_ylim(-y_lim_range, y_lim_range)
axs[2].grid(True, 'both', ls=':', alpha=0.7)
# axs[2].legend(loc='lower right', framealpha=0.75)

axs[0].set_ylabel('Planetographic Latitude')
axs[0].set_xlabel('Planetographic Longitude'), axs[1].set_xlabel('Planetographic Longitude')

plt.show()

#%%


plt.figure(figsize=(3,6), dpi=175)
plt.plot(hst_wind_1D['x'], hst_wind_1D['y'], color='b', lw=0.5, label='HST 2015')
plt.plot(zwp_5microns, lat_array, marker='s', ms=2, color='k', lw=0.5, label='IRTF 2024')
plt.plot(smooth_zwp_5microns, lat_array, color='r', lw=0.9, label='Smoothed')
plt.ylim(-y_lim_range, y_lim_range)
plt.grid(True, 'both', ls=':', alpha=0.7)
plt.legend(loc='upper right', framealpha=0.75)
plt.xlabel('Zonal velocity, m/s'), plt.ylabel('Planetographic Latitude')
plt.show()


#%%


print(cmap_meta_list[0].shape)
ZWP_Class.checkMapResolution(cmap_meta_list[0], np.ones((200, 200)))
