#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 13:57:20 2025

@author: ahyder
"""

# using planetmapper to plot some spex data

import planetmapper
from planetmapper.kernel_downloader import download_urls
import pylab as plt
import numpy as np
from astropy.io import fits
import math

# issued once. Already downloaded to ~/spice_kernels
# download_urls('https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/')
# print('lsk done...')
# download_urls('https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/')
# print('pck done...')

# planetary ephemeris - needed once; already downloaded
# Locations of planetary system barycentres:
# download_urls('https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp')
# Locations of Jupiter and its major satellites:
# download_urls('https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/jup365.bsp')


#  HST kernel - downloaded.
# download_urls('https://naif.jpl.nasa.gov/pub/naif/HST/kernels/spk/')

print('All appropriate kernels downloaded to ~/spice_kernels. Jupiter, major satellites, and HST information already stored...')

file_list = ['jc055090.fits.gz', 'jc277312.fits.gz', 'jc451486.fits.gz', 'jc505516.fits.gz']

# jc451486:
# disc_r0 = 198.25 # doesnt really change within a single day so the scale remains the same as well
# img_scale_kmpix     = 360.61538
# img_scale_arcsecpix = 0.11354471

# Number of images
num_images = len(file_list)
ncols = math.ceil(math.sqrt(num_images))
nrows = math.ceil(num_images / ncols)

fig1, axs1 = plt.subplots(
nrows=nrows, ncols=ncols, figsize=(7, 7), dpi=150)
axs1 = axs1.flatten()

fig2, axs2 = plt.subplots(
nrows=nrows, ncols=ncols, figsize=(10, 6), dpi=150)
axs2 = axs2.flatten()

for ax1, f in zip(axs1, file_list):
    print(f, type(f))
    
    hdul = fits.open('spex/2024oct22/' + f)
    hdr = hdul[0].header
    print(hdr['DATE_OBS'], hdr['TIME_OBS'])
    utc_str = hdr['DATE_OBS'] + 'T' + hdr['TIME_OBS']

    disc_x, disc_y = hdr['CX'], hdr['CY'] # obtained automatically so the disc will be around the correct point
    print(disc_x, type(disc_x))
    print()


    # navigational data made using the GUI - doing it here results in the NASA IRTF being set as observer.
    # Any spice available for it though? Would make things easier!
    observation = planetmapper.Observation('spex/2024oct22/' + f, 
                                           target='jupiter', 
                                           utc=utc_str, 
                                           observer='EARTH')

    # body = planetmapper.Body('jupiter', utc_str, observer='EARTH')
    # body = planetmapper.BodyXY('jupiter', utc_str, observer='EARTH')
    # body.plot_wireframe_radec(show=False)
    # body.set_disc_params(x0=disc_x, y0=disc_y)

    # Run the GUI to fit the disc interactively
    #observation.run_gui() # good for fitting the disc directly
    disc_r0 = 198.25


    # plt.figure(dpi=150)
    # observation.plot_backplane_img('LON-GRAPHIC')
    # plt.show()
    
    

    ax1.imshow(observation.data[0, ...], origin='lower', cmap='gist_heat')
    observation.set_disc_params(x0=disc_x, y0=disc_y, r0=disc_r0)
    observation.plot_wireframe_xy(ax1, formatting={
                                      'grid': {'linestyle': '-', 'linewidth': 0.5, 'alpha': 0.3, 'color':'w'},
                                      'prime_meridian': {'linewidth': 1, 'color': 'r'},
                                      'terminator': {'linewidth': 0.5, 'color':'w'},
                                      'equator': {'linewidth': 2, 'color':'w'},
                                      'limb_illuminated': {'linewidth': 0.5, 'color':'w'},},
                                    indicate_equator=True,
                                )
    ax1.set_title(hdr['TIME_OBS']+'UTC')
    # plt.show()
fig1.tight_layout()
    
for ax2, f in zip(axs2, file_list):
    
    hdul = fits.open('spex/2024oct22/' + f)
    hdr = hdul[0].header
    print(hdr['DATE_OBS'], hdr['TIME_OBS'])
    utc_str = hdr['DATE_OBS'] + 'T' + hdr['TIME_OBS']
    disc_x, disc_y = hdr['CX'], hdr['CY']
    
    observation = planetmapper.Observation('spex/2024oct22/' + f, 
                                            target='jupiter', 
                                            utc=utc_str, 
                                            observer='EARTH')
    observation.set_disc_params(x0=disc_x, y0=disc_y, r0=disc_r0)
    
    # Plot a mapped RGB image of the data in the top right
    degree_interval = 1.0  # Plot maps with 4 pixels/degree
    emission_cutoff = 80

    mapped_data = observation.get_mapped_data('cubic')  # get the mapped data

    # Display mapped image and add a useful annotation
    observation.imshow_map(mapped_data[0, ...], ax2, cmap='inferno')
    ax2.set_title(hdr['TIME_OBS']+'UTC')
    # plt.figure(dpi=125)
    # ax2.contourf(np.linspace(360,0,360), np.linspace(-90,90,180), mapped_data[0, ...], 256, cmap='magma')
    # ax2.invert_xaxis()
    # plt.grid(True, 'both', ls=':')
    # ax2.title('Day and Time: ' + '\n' + hdr['DATE_OBS'] + ' ' + hdr['TIME_OBS'] + ' UTC')
    # ax2.set_xlabel('Planetographic longitude (W)'), plt.ylabel('Planetographic latitude')

# fig2.tight_layout()

plt.show()




