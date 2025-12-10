#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 13:32:22 2025

@author: ahyder
"""

# import planetmapper
from planetmapper.kernel_downloader import download_urls
# import matplotlib
import matplotlib.pyplot as plt
import numpy as np
# import math
from astropy.io import fits
# import time
from scipy.signal import correlate
from scipy.interpolate import interp2d
from scipy.interpolate import RegularGridInterpolator
import scipy as sp
# import pytz
import datetime
# import csv
# import os
from scipy import fftpack, ndimage, stats, signal

R_J   = 71492e3 # m
R_J_p = 66854e3 # m

oblateness = (R_J - R_J_p) / R_J # Jupiter oblateness
eccentricity = np.sqrt((R_J**2 - R_J_p**2) / R_J**2) # Jupiter eccentricity

class ZWP_Class:
    def __init__(self):
        pass

    @staticmethod
    def closestIndex2Lat(latarr_aux, lat_i):
        # take in a specific latitude, and spit out the index and value in the given latitude array that is closest to it (helps in indexing and binning later)
        closest_value_in_array = latarr_aux[min(range(len(latarr_aux)), key=lambda i: abs(latarr_aux[i]-lat_i))]
        return np.where(np.array(latarr_aux)==closest_value_in_array)[0][0], closest_value_in_array
    
    @staticmethod
    def outputnaming(file_paths_list, year_str, test_i):
        # spit out the file path naming this will be used as the basis for outputs
        fn1, fn2 = file_paths_list[test_i][0].split('/')[-1].split('.cmap')[0], file_paths_list[test_i][1].split('/')[-1].split('.cmap')[0]
        naming_str = year_str + '_' + fn1 + '_' + fn2
        return naming_str


    @staticmethod
    def planeto_Centric2Graphic(R_J, R_J_p, centric_lat_array):
        # convert planetocentric latitudinal coordinates to planetographic
        centric_lat_array_rad = centric_lat_array * np.pi/180
        ratio_radii = (R_J_p/R_J)**2
        return np.arctan(np.tan(centric_lat_array_rad) / ratio_radii)
    
    @staticmethod
    def planeto_Graphic2Centric(R_J, R_J_p, graphic_lat_array):
        # convert planetographic latitudinal coordinates to planetocentric
        graphic_lat_array_rad = graphic_lat_array * np.pi/180
        ratio_radii = (R_J_p/R_J)**2
        return np.arctan(np.tan(graphic_lat_array_rad) * ratio_radii)
        
    
    @staticmethod
    def checkMapResolution(map1, map2, print_out=False):
        check_ind = 0
        if map1.size == map2.size:
            if print_out:
                print('Resolution is the same; check_ind = ', check_ind, ' and shapes: ', map1.shape, map2.shape)
                print()
        else:
            if map1.size > map2.size:
                check_ind = 1
                if print_out:
                    print('Please interpolate MAP1 to a lower resolution to continue; check_ind = ', check_ind, ' and shapes: ', map1.shape, map2.shape)
                    print()
            else:
                check_ind = 2
                if print_out:
                    print('Please interpolate MAP2 to a lower resolution to continue; check_ind = ', check_ind, ' and shapes: ', map1.shape, map2.shape)
                    print()
        return check_ind
    
    @staticmethod
    def makeMapsCompatible(file_paths_list, smart_index_list, long_res_base=360, lat_res_base=180, emission_angle_correction=False, execution_path_bool=False, fits_str='.mu.fits.gz'):
        '''
        Get the cylindrical map list (binary pairs), and make sure each pair has a consistent resolution. Either base (180, 360) or doubled.
        '''
        compatible_cmaps = []
        long_array_base, lat_array_base = np.linspace(360, 0,long_res_base), np.linspace(-90, 90,lat_res_base)        

        # Update if need be but try to keep things where they are. execution_path is the directory in which /prvt/ exists.
        execution_path = '/Users/ahyder/Desktop/windfieldproject'

        for i_f, f in enumerate(file_paths_list):
            
            if execution_path_bool:
                hdul0 = fits.open(execution_path + f[0])
                hdul1 = fits.open(execution_path + f[1])
            else:
                hdul0 = fits.open(f[0])
                hdul1 = fits.open(f[1])
            
            # in case emission angle should be divided out:
            if emission_angle_correction:
                if execution_path_bool:
                    emm0 = fits.open((execution_path + f[0]).split('.cmap')[0] + fits_str)
                    emm1 = fits.open((execution_path + f[1]).split('.cmap')[0] + fits_str)
                else:
                    emm0 = fits.open(f[0].split('.cmap')[0] + fits_str)
                    emm1 = fits.open(f[1].split('.cmap')[0] + fits_str)
                
            try:
                if smart_index_list[i_f] == 1:
                    long_res_hi = hdul0[0].header['NAXIS1']
                    lat_res_hi  = hdul0[0].header['NAXIS2']
                    data_res_hi = hdul0[0].data
                    data_from_other_map = hdul1[0].data
                    # divide through by the emission angle after the nans have been removed there as well
                    if emission_angle_correction:
                        emmi_res_hi, emmi_from_other_map = emm0[0].data, emm1[0].data
        
                elif smart_index_list[i_f] == 2:
                    long_res_hi = hdul1[0].header['NAXIS1']
                    lat_res_hi  = hdul1[0].header['NAXIS2']
                    data_res_hi = hdul1[0].data
                    data_from_other_map = hdul0[0].data
                    # divide through by the emission angle after the nans have been removed there as well
                    if emission_angle_correction:
                        emmi_from_other_map, emmi_res_hi = emm0[0].data, emm1[0].data
        
                else:
                    # make np.nan where map is zero (no data)
                    aux0_data, aux1_data = hdul0[0].data, hdul1[0].data
                    aux0_data[aux0_data == 0], aux1_data[aux1_data == 0] = np.nan, np.nan
                    
                    if emission_angle_correction:
                        aux0_emmi, aux1_emmi = emm0[0].data, emm1[0].data
                        aux0_emmi[aux0_emmi == 0], aux1_emmi[aux1_emmi == 0] = np.nan, np.nan
                        aux0_data, aux1_data = aux0_data/aux0_emmi, aux1_data/aux1_emmi

                    compatible_cmaps.append([aux0_data, aux1_data])
                    continue  # Skip to the next iteration
                    
                # Generate grid for interpolation (ensure correct order)
                lon_grid = np.linspace(360, 0, long_res_hi)
                lat_grid = np.linspace(-90, 90, lat_res_hi)
                # print('long_res_hi, lat_res_hi: ', long_res_hi, lat_res_hi)
                
                # Create interpolation function
                aux_interp = RegularGridInterpolator((lat_grid, lon_grid), data_res_hi, 
                                                     method='linear', bounds_error=False, fill_value=None)
                # Ensure lat_array and long_array are correctly broadcasted
                lon_interp, lat_interp = np.meshgrid(long_array_base, lat_array_base, indexing='ij')
                
                # Interpolate to lower resolution - ensure to transpose due to how 'ij' indexing is ordered in the interpolation
                data_res_low = aux_interp((lat_interp, lon_interp))
                
                # make np.nan where map is zero (no data)
                data_res_low[data_res_low == 0] = np.nan
                data_from_other_map[data_from_other_map == 0] = np.nan
                
                if emission_angle_correction:
                    # Create interpolation function for emmision angle, then replace zeros and divide the data by emission values.
                    emm_interp = RegularGridInterpolator((lat_grid, lon_grid), emmi_res_hi, 
                                                         method='linear', bounds_error=False, fill_value=None)
                    emmi_res_low = emm_interp((lat_interp, lon_interp))
                    emmi_res_low[emmi_res_low == 0] = np.nan
                    emmi_from_other_map[emmi_from_other_map == 0] = np.nan
                    data_res_low = data_res_low/emmi_res_low
                    data_from_other_map = data_from_other_map/emmi_from_other_map
                
                # follow the original order of the input file list...
                if smart_index_list[i_f] == 1:
                    compatible_cmaps.append([data_res_low.T, data_from_other_map])
                elif smart_index_list[i_f] == 2:
                    compatible_cmaps.append([data_from_other_map, data_res_low.T])
                # smart_index_list[i_f] == 0 case already handled in previous else.
        
            finally:
                hdul0.close()
                hdul1.close()
        return compatible_cmaps
    
    @staticmethod
    def getTemporalDifference(file_list):
        '''
        Returns the total seconds elapsed between observations
        '''
        d_obj = []
        for i_f, f in enumerate(file_list):
            hdul = fits.open(f)
            hdr = hdul[0].header
            date_obs_str = hdr['DATE_OBS'].split("-")
            time_obs_str = hdr['TIME_OBS'].split(":")
            d_obj.append(datetime.datetime(int(date_obs_str[0]), int(date_obs_str[1]), int(date_obs_str[2]), 
                                           int(time_obs_str[0]), int(time_obs_str[1]), int(np.round(float(time_obs_str[2])))))
        return (d_obj[-1] - d_obj[0]).total_seconds()

    @staticmethod
    def getDeltaXfromLat(delta_longitude, lat_array, lat_value=None):
        '''
        Get the latitudinally dependent ∆x from ∆-longitude
        '''
        if lat_value:
            print('Speed at latitude: ', lat_value)
            # R_J*(1 - oblateness*np.sin(lat_array*np.pi/180)**2)
        radius_in_lat = R_J*(1 - oblateness*np.sin(lat_array*np.pi/180)**2) # radius in meters at specific latitude
        return (delta_longitude / 360) * 2*np.pi*radius_in_lat # in meters

    @staticmethod
    def getSpeedatLat(delta_longitude, lat_array, delta_t, lat_value=None):
        if lat_value:
            print('Speed at latitude ', lat_value, '˚: ', ZWP_Class.getDeltaXfromLat(delta_longitude, lat_array)/delta_t, ' m/s')
        else:
            return ZWP_Class.getDeltaXfromLat(delta_longitude, lat_array)/delta_t # m/s

    @staticmethod
    def getValidLatLonRange(cyl_map1, cyl_map2):
        
        # Ensure maps have the same shape
        assert cyl_map1.shape == cyl_map2.shape, "Maps must have the same dimensions"
        
        # Create an intersection nd-array where data exists in both maps
        intersection = (~np.isnan(cyl_map1) & ~np.isnan(cyl_map2)).astype(int)
        
        return intersection

    @staticmethod
    def shiftMap(unshifted_long, unshifted_map, shift_long):
    
        mean_step = int(np.round(shift_long / abs(unshifted_long[1] - unshifted_long[0])))
        shifted_long = np.roll(unshifted_long, -mean_step, axis=0)
    
        if len(unshifted_map.shape) == 2:
            shifted_map = np.roll(unshifted_map, -mean_step, axis=1)
        elif len(unshifted_map.shape) == 1:
            shifted_map = np.roll(unshifted_map, -mean_step, axis=0)
    
        # Ensure the longitude array is wrapped properly to maintain order
        shifted_long[shifted_long < 0] += 360
        shifted_long[shifted_long > 360] -= 360
    
        return shifted_long, shifted_map
    
    @staticmethod
    def defineLongThresh(BIM, long_array, lat_array, crit_lat = 60, lat_res_tol = 1e-2):
        '''
        Obtain the indices of the longitude array that can be used for ZWP generation.
        This function supplements the definition of the BIM, as it limits the latitude range.
        '''
         # threshold should be modified if higher resolution data is used.
        lat_init_ind      = np.where(np.isclose(lat_array, -crit_lat, lat_res_tol))
        lat_endd_ind      = np.where(np.isclose(lat_array, crit_lat, lat_res_tol)) 
        lat_delta_indices = [lat_init_ind[0][0], lat_endd_ind[0][0]+1] # start and end indices for critical range
        
        # define the longitude index mask that will be populated in the loop; this is a boolean array
        longThresh_mask = np.zeros(len(long_array))
        
        # maximum value of summed longitude within the critical latitudes:
        maxLatSumAllowedinCrit = len(BIM[lat_delta_indices[0]:lat_delta_indices[1], 0])
        
        for i in range(len(long_array)):
            sum_at_given_long = np.sum(BIM[lat_delta_indices[0]:lat_delta_indices[1], i])
            if (int(sum_at_given_long) == maxLatSumAllowedinCrit):
                longThresh_mask[i] = True
            else:
                longThresh_mask[i] = False
    
        # the mask must be inverted to get the appropriate longitudes using np.logical_not
        return longThresh_mask # 1D array
    
    @staticmethod
    def applyMargin(longThresh_mask, long_array, margin=2):
        '''
        Use some margin to limit how close you get to the NaNs.
        Optional to use, but may limit the uncertainty in the correlations generated from data too close to the limb of the disk.
        '''
    
        margin_applied = np.copy(longThresh_mask) # shallow copy of primary boolean mask
        delta_margin_in_indices = int(np.round(margin / abs(long_array[1] - long_array[0])))
        start_margin = np.where(abs(np.diff(longThresh_mask)) == 1)[0][0]
        stop_margin  = np.where(abs(np.diff(longThresh_mask)) == 1)[0][1] + 1 # +1 so the index is on True, consistently.
        print('start and stop margins (prior to shifting): ', start_margin, stop_margin)
        
        # For values within the margin, set to False for each of the wrapping cases:
        if (longThresh_mask[0] == True) & (longThresh_mask[-1] == True):
            # check if the map is wrapped:
            # [1, 1, 0, ..., 0, 1, 1], e.g.
            margin_applied[start_margin:start_margin - delta_margin_in_indices:-1]  = False
            margin_applied[stop_margin:stop_margin + delta_margin_in_indices:1]     = False
        
        else:
            # case where the longitudes were not wrapped:
            # [0, 0, 1, ..., 1, 0, 0] or [1, 1, 1, ..., 0, 0] or [0, 0, ..., 0, 1, 1], e.g.
            margin_applied[start_margin:start_margin + delta_margin_in_indices:1]  = False
            margin_applied[stop_margin:stop_margin - delta_margin_in_indices:-1]     = False
    
        return margin_applied # 1D array

    @staticmethod
    def getInterpolatedArray(x, y, y_new):
        f_lin = sp.interpolate.interp1d(x, y, kind = 'linear', axis = - 1, copy = True, fill_value = 'extrapolate', assume_sorted = True)
        return f_lin(y_new)

    @staticmethod
    def getSmoothedProfile(unsmooth_1D, lat_array, filter_type='savgol', jet_width=5, poly_order=2, mode='interp'):
        '''
        Apply the Savitzky-Golay filter to the final 1D profile
        '''
        delta_lat = abs(lat_array[1] - lat_array[0])
        window_len = int(np.ceil(jet_width/delta_lat) // 2 * 2 + 1) # window len will change based on resolution and will be odd.
        if filter_type=='savgol':
            return sp.signal.savgol_filter(unsmooth_1D, window_len, poly_order, mode=mode)
        elif filter_type=='median':
            return sp.ndimage.median_filter(unsmooth_1D, size=window_len)
        else:
            raise ValueError("Select appropriate filter type: 'savgol', 'median'")

    @staticmethod
    def generateErrors(shifted_f1, shifted_f2, estimated_lag, lat_array, lat_delta_indices, margin_applied_in_func_shifted, max_lag, delta_long_step):
        
        # use the mean lag for the full pair (remove the sliding window dependence if it exists)
        if len(estimated_lag.shape) == 2:
            # estimated_lag is in ∆-longitudes so divide by delta_long_step to get back to ∆-pixels
            estimated_lag_across_lat = np.mean(estimated_lag, axis=1) / delta_long_step
        else:
            estimated_lag_across_lat = estimated_lag / delta_long_step
            
        error_array = np.zeros(len(lat_array))
        fwhm_array  = np.zeros(len(lat_array))
        
        for l in range(len(lat_array[lat_delta_indices[0]:lat_delta_indices[1]])):
            A = shifted_f1[l+lat_delta_indices[0], margin_applied_in_func_shifted.astype(bool)]
            B = shifted_f2[l+lat_delta_indices[0], margin_applied_in_func_shifted.astype(bool)]
            
            # shift A forward by ∆-pix//2:
            A_prime = np.roll(A, shift = -int(estimated_lag_across_lat[l+lat_delta_indices[0]] / 2))
            # shift B backwards by ∆-pix//2:
            B_prime = np.roll(B, shift = int(estimated_lag_across_lat[l+lat_delta_indices[0]] / 2))
            
            # Standardize the primes for correlation test...
            A_prime_mean, A_prime_std = np.mean(A_prime), np.std(A_prime)
            B_prime_mean, B_prime_std = np.mean(B_prime), np.std(B_prime)
            A_prime_norm = (A_prime - A_prime_mean) / A_prime_std
            B_prime_norm = (B_prime - B_prime_mean) / B_prime_std

            # conventional normalization if needed. Uncomment if using the following and comment above...
            # A_prime_norm = A_prime/np.nanmax(A_prime)
            # B_prime_norm = B_prime/np.nanmax(B_prime)
            
            # Compute cross-correlation between the half-shifted arrays
            correlation_result = correlate(A_prime_norm, B_prime_norm, mode='full', method='auto')
            lags = np.arange(-len(A_prime) + 1, len(B_prime)) # total amount of possible lags that can be investigate prior to max_lag limit
            
            # Apply max lag threshold
            valid_indices = (lags >= -max_lag) & (lags <= max_lag)
            correlation_limited = correlation_result[valid_indices]
            lags_limited = lags[valid_indices]
            
            # Find the lag where correlation peaks (error in estimated shift)
            error_in_estimated_lag_across_lat_in_pix = lags_limited[np.argmax(correlation_limited)]
            
            # Convert to ∆-longitude uncertainty and populate output array
            error_array[l+lat_delta_indices[0]] = error_in_estimated_lag_across_lat_in_pix * delta_long_step
            
            # Compute Full Width at Half Maximum (FWHM) as a more appropriate measure of confidence
            half_max = np.max(correlation_limited) / 2
            fwhm_lags = lags_limited[correlation_limited > half_max] # indices where the correlation is above the FWHM
            # 'fwhm' is the pixel range of the peak correlation - shows how unsure it is by capturing the extent of the spread
            fwhm_in_pix = np.max(fwhm_lags) - np.min(fwhm_lags) if len(fwhm_lags) > 1 else 0 # Zero if peak is very well defined
            # fwhm_in_degrees
            fwhm_array[l+lat_delta_indices[0]] = fwhm_in_pix * delta_long_step

        return error_array, fwhm_array
    
    @staticmethod
    def fullCorrelation(shifted_f1, shifted_f2, lat_array, lat_delta_indices, margin_applied_in_func_shifted, final_longitudes_masked_shifted, max_lag, max_corr_lim, delta_long_in_lat, delta_long_step):

        # Main correlation loop...sort of like a kai-square minimization.
        # Within these latidude range, the values are defined and away from the NaNs
        print('Entered sub function')
        for l in range(len(lat_array[lat_delta_indices[0]:lat_delta_indices[1]])):
            # the .astype(bool) just converts the [0, 0, 1, 1] to [False, False, True, True]
            shifted_f1_valid = shifted_f1[l+lat_delta_indices[0], margin_applied_in_func_shifted.astype(bool)]
            shifted_f2_valid = shifted_f2[l+lat_delta_indices[0], margin_applied_in_func_shifted.astype(bool)]
            # valid_longitudes = final_longitudes_masked_shifted[margin_applied_in_func_shifted.astype(bool)]
            
            # Compute mean and standard deviation of valid values (for single lat)
            f1_mean, f2_mean = np.mean(shifted_f1_valid), np.mean(shifted_f2_valid)
            f1_std, f2_std   = np.std(shifted_f1_valid), np.std(shifted_f2_valid)
            
            # Standardize valid data
            f1_normalized = (shifted_f1_valid - f1_mean) / f1_std
            f2_normalized = (shifted_f2_valid - f2_mean) / f2_std

            # Normalize valid data using conventional normalization method. Comment if not using and vice versa with above...
            # f1_normalized = shifted_f1_valid / np.nanmax(shifted_f1_valid)
            # f2_normalized = shifted_f2_valid / np.nanmax(shifted_f2_valid)
    
            # Perform 1D cross-correlation between the two normalized maps
            corr = correlate(f1_normalized, f2_normalized, mode='full', method='auto')
            lags = np.arange(-len(shifted_f1_valid) + 1, len(shifted_f2_valid))
    
            # Limit the lag range using the physical limit argument
            lag_limit_indices = (lags >= -max_lag) & (lags <= max_lag)
            corr_limited      = corr[lag_limit_indices]
            lags_limited      = lags[lag_limit_indices]
            
            # Find the lag with maximum correlation
            max_corr_index = np.argmax(corr_limited) # index of best correlative strength within the physically limited lag range
            best_lag       = lags_limited[max_corr_index]
            
            # save the maximum correlation coefficient as a measure of confidence
            pearson_corr = np.corrcoef(f1_normalized, f2_normalized)[0, 1]
            max_corr_lim[l+lat_delta_indices[0]] = pearson_corr #np.max(corr_limited)
    
            # Compute the angular shift; convert from pixel lag to a physical lag
            delta_long_in_lat[l+lat_delta_indices[0]] = best_lag * delta_long_step  # Assumes uniform spacing
        
        return delta_long_in_lat, max_corr_lim
    
    
    @staticmethod
    def SlidingWindowCorrelation(sliding_window_len, window_step, shifted_f1, shifted_f2, lat_array, lat_delta_indices, margin_applied_in_func_shifted, final_longitudes_masked_shifted, max_lag, delta_long_step):
        '''
        Enter if using the sliding window approach. A smaller window is defined within each image pair using 
        margin_applied_in_func_shifted. Each window will result in an individual correlation
        '''
        window_step_in_pixels        = int(np.ceil(window_step / delta_long_step))
        sliding_window_len_in_pixels = int(np.ceil(sliding_window_len / delta_long_step))  
        # print(window_step_in_pixels, sliding_window_len_in_pixels)
        
        # indices of where the full longitude range is defined post-margin and shift.
        loc = np.where(abs(np.diff(margin_applied_in_func_shifted)) == 1)[0]
        # how many shifts are possible given the sliding window size and the size of the slide
        num_of_possible_shifts = (loc[-1] - loc[0] - sliding_window_len_in_pixels) // window_step_in_pixels
        
        # in_func_counter always ends up being 110 (for sliding_window_len=30) as both, sliding_window_len_in_pixels and 
        # window_step_in_pixels are increased when the resolution is increased. So the amount of steps is the same.
        # +1 is needed because num_of_possible_shifts is 109 (excludes the final step so it must be added explicitly).
        delta_long_in_lat_slide = np.zeros((len(lat_array), num_of_possible_shifts + 1))
        max_corr_lim_slide      = np.zeros((len(lat_array), num_of_possible_shifts + 1))
        
        in_func_counter = 0 # serves as the primary index for the num_of_possible_shifts since we have used a while loop
        RHS_bound = loc[0] # right bound of sliding window
        
        print('possible shifts :', num_of_possible_shifts, 'window len in pix', sliding_window_len_in_pixels, 'window step in pix', window_step_in_pixels, 'loc 0, -1, and diff. ', loc[0], loc[-1], loc[-1]-loc[0])
        
        while RHS_bound + sliding_window_len_in_pixels <= loc[-1]:
            aux_bool = np.zeros_like(margin_applied_in_func_shifted)
            aux_bool[RHS_bound + 1:RHS_bound + 1 + sliding_window_len_in_pixels] = 1
            RHS_bound += window_step_in_pixels
            # print(in_func_counter, RHS_bound)
            
            # Main correlation loop...sort of like a kai-square minimization - being used for each window here
            # Within these latidude range, the values are defined and away from the NaNs
            for l in range(len(lat_array[lat_delta_indices[0]:lat_delta_indices[1]])):
                # the .astype(bool) just converts the [0, 0, 1, 1] to [False, False, True, True]
                shifted_f1_valid = shifted_f1[l+lat_delta_indices[0], aux_bool.astype(bool)]
                shifted_f2_valid = shifted_f2[l+lat_delta_indices[0], aux_bool.astype(bool)]
                # valid_longitudes = final_longitudes_masked_shifted[margin_applied_in_func_shifted.astype(bool)]
                
                # Compute mean and standard deviation of valid values (for single lat)
                f1_mean, f2_mean = np.mean(shifted_f1_valid), np.mean(shifted_f2_valid)
                f1_std, f2_std   = np.std(shifted_f1_valid), np.std(shifted_f2_valid)
                
                # Standardize valid data
                f1_normalized = (shifted_f1_valid - f1_mean) / f1_std
                f2_normalized = (shifted_f2_valid - f2_mean) / f2_std

                # Normalize valid data using conventional normalization method. Comment if not using and vice versa with above...
                # f1_normalized = shifted_f1_valid / np.nanmax(shifted_f1_valid)
                # f2_normalized = shifted_f2_valid / np.nanmax(shifted_f2_valid)
        
                # Perform 1D cross-correlation between the two normalized maps
                corr = correlate(f1_normalized, f2_normalized, mode='full', method='auto')
                lags = np.arange(-len(shifted_f1_valid) + 1, len(shifted_f2_valid))
        
                # Limit the lag range using the physical limit argument
                lag_limit_indices = (lags >= -max_lag) & (lags <= max_lag)
                corr_limited      = corr[lag_limit_indices]
                lags_limited      = lags[lag_limit_indices]
                
                # Find the lag with maximum correlation
                max_corr_index = np.argmax(corr_limited) # index of best correlative strength within the physically limited lag range
                best_lag       = lags_limited[max_corr_index]
                
                # save the maximum correlation coefficient as a measure of confidence
                pearson_corr = np.corrcoef(f1_normalized, f2_normalized)[0, 1]
                # max_corr_lim_slide.append([pearson_corr]) #np.max(corr_limited)
                max_corr_lim_slide[l+lat_delta_indices[0], in_func_counter] = pearson_corr #np.max(corr_limited)
        
                # Compute the angular shift; convert from pixel lag to a physical lag
                delta_long_in_lat_slide[l+lat_delta_indices[0], in_func_counter] = best_lag * delta_long_step  # Assumes uniform spacing
            
            # update the index for num_of_possible_shifts (axis = 1 for delta_long_in_lat_slide and max_corr_lim_slide)
            in_func_counter += 1
        print('final in_func_counter: ', in_func_counter)
        return delta_long_in_lat_slide, max_corr_lim_slide
    
    @staticmethod
    def detrend_ignore_nan(x):
        x = np.asarray(x, float)
        mask = ~np.isnan(x)
        if mask.sum() < 2:
            return x  # not enough points to detrend; return as-is
        # detrend valid part
        x_valid = signal.detrend(x[mask], type='linear')
        out = x.copy()
        out[mask] = x_valid
        return out

    @staticmethod
    def spectral_entropy_1d(x, eps=1e-12):
        # power spectrum (one-sided)
        X = np.abs(fftpack.fft(x))**2
        X = X[:len(X)//2]
        P = X / (X.sum() + eps)
        H = -np.sum(P * np.log(P + eps))
        # normalize entropy to [0,1] by dividing by log(N)
        return H / np.log(len(P) + eps)

    @staticmethod
    def hf_power_fraction(x, cutoff_frac=0.5):
        # fraction of power above cutoff_frac * Nyquist (0..1)
        N = len(x)
        X = np.abs(fftpack.fft(x))**2
        Xh = X[:N//2]
        cutoff_idx = int(np.floor(cutoff_frac * len(Xh)))
        if cutoff_idx < 1:
            return 0.0
        return Xh[cutoff_idx:].sum() / (Xh.sum() + 1e-12)

    @staticmethod
    def laplacian_energy(x):
        # mean energy of second difference - actually not really a laplacian energy per se...
        d2 = np.diff(x, n=2)
        return np.mean(d2**2)

    @staticmethod
    def total_variation(x):
        return np.mean(np.abs(np.diff(x)))

    @staticmethod
    def mad(x):
        return np.median(np.abs(x - np.median(x)))

    @staticmethod
    def robust_scale_norm(arr):
        # map arr to roughly 0..1 using median/IQR then sigmoid
        med = np.median(arr)
        q1, q3 = np.percentile(arr, [25,75])
        iqr = max(q3 - q1, 1e-12)
        z = (arr - med) / iqr
        # use a sigmoid to compress
        return 1.0 / (1.0 + np.exp(-z))
    
    @staticmethod
    def irregularity_scores(image, smooth_sigmas=(0.0, 1.0, 3.0), hf_cutoff_frac=0.4, lat_smooth_sigma=1.0):
        """
        image : 2D numpy array of shape (n_lat, n_lon)
        smooth_sigmas : sequence of gaussian sigmas to try (multi-scale)
        returns: score vector length n_lat (0..1)
        """
        n_lat, n_lon = image.shape
        # store metrics per-latitude
        m_std = np.zeros(n_lat)
        m_mad = np.zeros(n_lat)
        m_tv = np.zeros(n_lat)
        m_lap_ms = np.zeros(n_lat)
        m_hf = np.zeros(n_lat)
        m_kurt = np.zeros(n_lat)
        m_peakcount = np.zeros(n_lat)

        for i in range(n_lat):
            row = image[i, :].astype(float)
            # baseline detrend (optional) to remove large-scale gradient:
            # row_d = signal.detrend(row, type='linear')
            row_d = ZWP_Class.detrend_ignore_nan(row)

            # simple metrics:
            m_std[i] = np.std(row_d)
            m_mad[i] = ZWP_Class.mad(row_d)
            m_tv[i] = ZWP_Class.total_variation(row_d)
            m_kurt[i] = stats.kurtosis(row_d, fisher=True, bias=False)  # Fisher => 0 for Gaussian

            # peak count as simple heuristic
            peaks, props = signal.find_peaks(row_d, height=np.std(row_d)*0.5, distance=max(1, n_lon//50))
            m_peakcount[i] = len(peaks)

            # laplacian energy at multiple scales (to capture narrow and broad)
            lap_vals = []
            for sigma in smooth_sigmas:
                if sigma > 0.0:
                    r_s = ndimage.gaussian_filter1d(row_d, sigma=sigma, mode='reflect')
                else:
                    r_s = row_d
                lap_vals.append(ZWP_Class.laplacian_energy(r_s))
            m_lap_ms[i] = np.max(lap_vals)  # max across scales

            # HF power fraction
            m_hf[i] = ZWP_Class.hf_power_fraction(row_d, cutoff_frac=hf_cutoff_frac)

        # Normalize each metric robustly to 0..1
        n_std = ZWP_Class.robust_scale_norm(m_std)
        n_mad = ZWP_Class.robust_scale_norm(m_mad)
        n_tv = ZWP_Class.robust_scale_norm(m_tv)
        n_lap = ZWP_Class.robust_scale_norm(m_lap_ms)
        n_hf = ZWP_Class.robust_scale_norm(m_hf)
        # kurtosis can be negative; map via absolute then norm
        n_kurt = ZWP_Class.robust_scale_norm(np.abs(m_kurt))
        n_peaks = ZWP_Class.robust_scale_norm(m_peakcount)

        # Combine: give more weight to metrics that detect spikes (lap, tv, hf, kurt)
        weights = {
            'std': 0.15,
            'mad': 0.10,
            'tv': 0.15,
            'lap': 0.25,
            'hf': 0.2,
            'kurt': 0.1,
            'peaks': 0.05
        }
        # compute weighted sum - this is kind of arbitrary but worth playing around with...test each individually by adjusting specific weights...
        combined = (weights['std'] * n_std +
                    weights['mad'] * n_mad +
                    weights['tv'] * n_tv +
                    weights['lap'] * n_lap +
                    weights['hf'] * n_hf +
                    weights['kurt'] * n_kurt +
                    weights['peaks'] * n_peaks)

        # optional smooth across latitudes to reduce single-lat noise
        if lat_smooth_sigma > 0:
            combined = ndimage.gaussian_filter1d(combined, sigma=lat_smooth_sigma, mode='reflect')

        # clip to 0..1
        combined = np.clip(combined, 0.0, 1.0)
        return combined
    
    @staticmethod
    def getFullLatitudinalCorrelation(cyl_map1, cyl_map2, long_array, lat_array, shift_long=90, crit_lat=60, margin=3, max_delta_long=10, force_shift=False, sliding_window_boolean=False, sliding_window_len=30, window_step=1, lat_res_tol=1e-2):
        '''
        Main function for zonal wind generation. Parameters defined below.
    
        ### Inputs:
        
        cyl_map1, cyl_map2:     Cylindrical maps (reduced and outputted from the planetmapper object)
        long_array, lat_array:  Full longitude and latitude of the images.
        shift_long:             Define map shift so that the correlation function can run properly. Defaults to 90˚.
        crit_lat:               Latitudinal range for the correlation. It will go from -crit_lat to +crit_lat. Defaults to 60˚.
                                Careful not to use just any value; often the maps have a varying latitudinal resolution so be specific.
        margin:                 How far you want to be from the NaNs. Also helps in limiting one to lower emission angles. Defaults to 3˚.
        max_delta_long:         Physical limit in longitude to how far correlative strengths should be checked. Defaults to 8˚.
        force_shift:            If True, force a shift by the amount defined by shift_long. Default is False.
        lat_res_tol:            Relative tolerance for how the critical latitudes are found. Typically needs to be varied depending on the 
                                resolution of the map, and specific latitude.
                               
        sliding_window_boolean: Using the sliding window technique as devised by Johnson+2018.
        sliding_window_len:     The window length in longitudinal degrees.
        window_step:            The sliding step size in longitudinal degrees.
    
        ### Returns:
    
        shifted_f1, shifted_f2:          Cylindrical maps shifted by shift_long to the West.
        shifted_BIM:                     Cylindrical boolean intersection map shifted by shift_long to the West.
        margin_applied_in_func_shifted:  Shifted 1D boolean array of margin-limited longitudes that define the validity zone for the correlation.
        final_longitudes_masked_shifted: Shifted 1D longitude array that is margin-limited.
        
        if not sliding_window_boolean:
            delta_long_in_lat:               ∆ in longitude for each latitude; result of the pixel lag determined via maximum correlative strength.
            max_corr_lim:                    Some measure of correlation strength (pearson coefficient, e.g.).
        
        if sliding_window_boolean:
            delta_long_in_lat_slide:         Same as delta_long_in_lat but for each window - so this output is 2D
            max_corr_lim_slide:              Same as max_corr_lim but for each window - so this output is 2D
        '''
        
        # first, check if the lat/long arrays have the appropriate size for the incoming maps. Sometimes, both are 360x720 so the hi-res is used.
        if cyl_map1.size > 180*360:
            long_array = np.linspace(360, 0, 720)
            lat_array = np.linspace(-90, 90, 360)
        else:
            pass
        
        BIM_in_func             = ZWP_Class.getValidLatLonRange(cyl_map1, cyl_map2) # (180, 360) Boolean intersection map
        longThresh_mask_in_func = ZWP_Class.defineLongThresh(BIM_in_func, long_array, lat_array, crit_lat = crit_lat) # longitudes from intersection
        margin_applied_in_func  = ZWP_Class.applyMargin(longThresh_mask_in_func, long_array, margin=margin) # (1, 360) boolean
    
        # (1, 360) array with the final useful longitudes in the map
        # make mask
        inverted_margin_array   = np.logical_not(margin_applied_in_func)
        # apply mask to create final useful longitudes, also (1, 360)
        final_longitudes_masked = np.ma.masked_array(long_array, inverted_margin_array)
    
        shift_necessity_bool = bool(margin_applied_in_func[0]) and bool(margin_applied_in_func[-1])
    
        if shift_necessity_bool or force_shift:
            if shift_necessity_bool:
                print('Shift applied due to map wrapping.')
                print()
            else:
                print('Shift applied by force. Proceed with caution.')
                print()
            # enter only if a shift is needed or you're deciding to force one for testing purposes.
            # apply linear shift (default is 90 but needs to be set for each observation cycle; this can be automated later)
            _, shifted_f1  = ZWP_Class.shiftMap(long_array, cyl_map1, shift_long) # 2D
            _, shifted_f2  = ZWP_Class.shiftMap(long_array, cyl_map2, shift_long) # 2D
            _, shifted_BIM = ZWP_Class.shiftMap(long_array, BIM_in_func, shift_long) # 2D
            _, margin_applied_in_func_shifted  = ZWP_Class.shiftMap(long_array, margin_applied_in_func, shift_long) # 1D
            _, final_longitudes_masked_shifted = ZWP_Class.shiftMap(long_array, final_longitudes_masked, shift_long) # 1D
        else:
            print('No shift applied.')
            print()
            # apply linear shift (default is 90 but needs to be set for each observation cycle; this can be automated later)
            _, shifted_f1  = ZWP_Class.shiftMap(long_array, cyl_map1, 0) # 2D
            _, shifted_f2  = ZWP_Class.shiftMap(long_array, cyl_map2, 0) # 2D
            _, shifted_BIM = ZWP_Class.shiftMap(long_array, BIM_in_func, 0) # 2D
            _, margin_applied_in_func_shifted  = ZWP_Class.shiftMap(long_array, margin_applied_in_func, 0) # 1D
            _, final_longitudes_masked_shifted = ZWP_Class.shiftMap(long_array, final_longitudes_masked, 0) # 1D
            
        # post-shift margin bounds:
        print('start and stop margins (after possible shift): ', np.where(abs(np.diff(margin_applied_in_func_shifted)) == 1)[0][0], np.where(abs(np.diff(margin_applied_in_func_shifted)) == 1)[0][1] + 1)
        # All maps, longitudes, and boolean arrays are now shifted by the same amount. If no shift was needed, only margins have been applied.
        # Now, the correlation can be performed by ignoring the NaNs around the data, but no wrapping/circular cross-correlation is needed.
        
        # get indices to limit the latitudinal extent:
        lat_init_ind      = np.where(np.isclose(lat_array, -crit_lat, lat_res_tol))
        lat_endd_ind      = np.where(np.isclose(lat_array, crit_lat, lat_res_tol)) 
        lat_delta_indices = [lat_init_ind[0][0], lat_endd_ind[0][0]+1] # start and end indices for critical range
    
        num_of_lat = len(lat_array)
        
        # define the lags within the region of applicability (where margin_applied_in_func_shifted
        # is 1, so the sum is the total number of valid longitudes).
    
        # previous definition is more general. The following is better for our case, as we know the physical limit of the lag.
        delta_long_in_lat = np.zeros(num_of_lat)
        max_corr_lim = np.zeros(num_of_lat)
        delta_long_step   = abs(long_array[150] - long_array[149])
        max_lag           = int(max_delta_long / delta_long_step) # maximum number of pixels allowed to be shifted to test 1D similarity.
        
        # returns must be inside the if-loop as the outputs delta_long_in_lat_slide and max_corr_lim_slide are of different shapes than
        # when sliding_window_boolean=False...
        if sliding_window_boolean:
            '''Enter sliding window mode - additional slide parameters must be defined'''
            delta_long_in_lat_slide, max_corr_lim_slide = ZWP_Class.SlidingWindowCorrelation(sliding_window_len, window_step, shifted_f1, shifted_f2, 
                                                                                             lat_array, lat_delta_indices, margin_applied_in_func_shifted, 
                                                                                             final_longitudes_masked_shifted, max_lag, delta_long_step)
            
            error_array, fwhm_array = ZWP_Class.generateErrors(shifted_f1, shifted_f2, delta_long_in_lat_slide, lat_array, lat_delta_indices, margin_applied_in_func_shifted, max_lag, delta_long_step)
            
            return shifted_f1, shifted_f2, shifted_BIM, margin_applied_in_func_shifted, final_longitudes_masked_shifted, delta_long_in_lat_slide, max_corr_lim_slide, error_array, fwhm_array
        
        else:
            '''Sub-function that calculates the correlation over the entire longitudinal range (no sliding window)'''
            delta_long_in_lat, max_corr_lim = ZWP_Class.fullCorrelation(shifted_f1, shifted_f2, lat_array, lat_delta_indices, 
                                                                        margin_applied_in_func_shifted, final_longitudes_masked_shifted, 
                                                                        max_lag, max_corr_lim, delta_long_in_lat, delta_long_step)
            
            error_array, fwhm_array = ZWP_Class.generateErrors(shifted_f1, shifted_f2, delta_long_in_lat, lat_array, lat_delta_indices, margin_applied_in_func_shifted, max_lag, delta_long_step)
            
            return shifted_f1, shifted_f2, shifted_BIM, margin_applied_in_func_shifted, final_longitudes_masked_shifted, delta_long_in_lat, max_corr_lim, error_array, fwhm_array






