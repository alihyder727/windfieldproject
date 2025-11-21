#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 10:04:09 2025

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
from scipy.interpolate import RegularGridInterpolator
import scipy as sp
import pytz
import datetime
import csv
import os

from zonalwindfunctions import ZWP_Class


#%%

long_res, lat_res = 360, 180
long_array, lat_array = np.linspace(360, 0,long_res), np.linspace(-90, 90,lat_res)
R_J = 71492e3 # m


main_path = 'spex/'
file_paths_list_51 = [] # with threshold set to 51 min (0.85)
t_diff_list_51 = []

with open(main_path+"pairs_2024_51min_cmap.csv", "r") as f:
    reader = csv.reader(f)
    for row in reader:
        selected_pair, t_diff_of_pair = [row[0], row[1]], float(row[2])
        t_diff_list_51.append(t_diff_of_pair)
        file_paths_list_51.append([selected_pair[0], selected_pair[1]])
        
for i in range(len(file_paths_list_51)):
    print(file_paths_list_51[i], t_diff_list_51[i])


r0_array = np.array([170.8, 162.7, 162.7, 162.7, 162.7, 162.7, 162.7, 162.7, 198.25, 198.25])

# full path to file
sample_file_ind = 0
file_list = [file_paths_list_51[sample_file_ind][0], file_paths_list_51[sample_file_ind][1]]
print()
print(file_list)
print()

#%%

'''
Super important index list! It is necessary for compatiblity of the maps.
0 - pair is already compatible.
1 - first map has a higher resolution and needs to be unresolved.
2 - second map has a higher resolution and needs to be unresolved.
'''

smart_index_list = np.zeros(len(file_paths_list_51))

for i_f, f in enumerate(file_paths_list_51):
    hdul0 = fits.open(f[0])
    hdul1 = fits.open(f[1])
    print(f[0].split('spex/')[-1], f[1].split('spex/')[-1])
    smart_index_list[i_f] = ZWP_Class.checkMapResolution(hdul0[0].data, hdul1[0].data)
    hdul0.close()
    hdul1.close()
    
#%%

print(smart_index_list)

for i_f, f in enumerate(file_paths_list_51):
    
    hdul0 = fits.open(f[0])
    hdul1 = fits.open(f[1])
    
    try:
        if smart_index_list[i_f] == 1:
            long_res_hi = hdul0[0].header['NAXIS1']
            lat_res_hi  = hdul0[0].header['NAXIS2']
            data_res_hi = hdul0[0].data

        elif smart_index_list[i_f] == 2:
            long_res_hi = hdul1[0].header['NAXIS1']
            lat_res_hi  = hdul1[0].header['NAXIS2']
            data_res_hi = hdul1[0].data

        else:
            continue  # Skip to the next iteration
            
        # Generate grid for interpolation (ensure correct order)
        lon_grid = np.linspace(360, 0, long_res_hi)
        lat_grid = np.linspace(-90, 90, lat_res_hi)
        print('long_res_hi, lat_res_hi: ', long_res_hi, lat_res_hi)
        
        # Create interpolation function
        aux_interp = RegularGridInterpolator((lat_grid, lon_grid), data_res_hi, 
                                             method='linear', bounds_error=False, fill_value=None)
        
        # Ensure lat_array and long_array are correctly broadcasted
        lon_interp, lat_interp = np.meshgrid(long_array, lat_array, indexing='ij')
        
        # Interpolate to lower resolution
        data_res_low = aux_interp((lat_interp, lon_interp))
        
        plt.figure(2, dpi=150)
        plt.subplot(211)
        plt.imshow(data_res_hi, origin='lower')
        plt.title(str(data_res_hi.shape), fontsize=8)
        plt.subplot(212)
        plt.imshow(data_res_low.T, origin='lower')
        plt.title(str(data_res_low.T.shape), fontsize=8)
        plt.show()
        plt.tight_layout()

    finally:
        hdul0.close()
        hdul1.close()
        
    

    


#%%

test_compatibility_func = ZWP_Class.makeMapsCompatible(file_paths_list_51, smart_index_list)

print(len(test_compatibility_func[0]))

for i in range(len(test_compatibility_func)):
    print(test_compatibility_func[i][0].shape, test_compatibility_func[i][1].shape, int(smart_index_list[i]))


#%%

test_i = 4
plt.figure(2, dpi=150)
plt.subplot(211)
plt.imshow(test_compatibility_func[test_i][0], cmap='gist_heat', origin='lower')
plt.title(str(test_compatibility_func[test_i][0].shape), fontsize=8)
plt.subplot(212)
plt.imshow(test_compatibility_func[test_i][1], cmap='gist_heat', origin='lower')
plt.title(str(test_compatibility_func[test_i][1].shape), fontsize=8)
plt.show()
plt.tight_layout()

#%%

print(len(test_compatibility_func))