#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 19:23:56 2018

@author: Steffen_KJ
"""

from __future__ import division
from __future__ import print_function
from pvextractor import extract_pv_slice, Path
import astropy.io.fits as pyfits
import numpy as np
import matplotlib.pyplot as plt


def load_fits_cube(cube_name):
    """
    This function loads a fits cube.

    Args:
        cube: spectral data cube, assumed dims are
        (im_x, im_y, n_chan, stokes_nr) or
        (im_x, im_y, n_chan).
    Returns:
        data, header
    """

    data_cube, hdr = pyfits.getdata(cube_name, 0, header=True)

    # Remove NaNs
    NaNs = np.isnan(data_cube)
    data_cube[NaNs] = 0.0

    assert len(data_cube.shape) == 4 or len(data_cube.shape) == 3, (
               'Wrong input cube dimensions. Cube must have 3 or 4 dimensions')
    # Account for the cube file format

    if len(data_cube.shape) == 4:
        assert data_cube.shape[0] == 1, 'Stokes nr expected to be 1.'
        stokes_nr, n_chan, ny, nx = data_cube.shape
        data_cube = data_cube[0, :, :, :]

    assert data_cube.shape[0] > 1, ('Must have more than 1 channel in cube. '
                                    'Did you use an image and not a cube? Then'
                                    ' use load_fits_im instead.')
    print('Shape returned from load_fits_cube is {}'.format(data_cube.shape))

    return data_cube, hdr


def load_fits_im(cube_name, is_PV=False):
    """
    This function loads a fits image.

    Args:
        cube: spectral data cube, assumed dims are [im_x, im_y, stokes_nr]
    Returns:
        data, header
    """

    data_cube_temp, hdr = pyfits.getdata(cube_name, 0, header=True)
    if len(data_cube_temp.shape) == 4:
        assert data_cube_temp.shape[0] == 1
        assert data_cube_temp.shape[1] == 1

    elif len(data_cube_temp.shape) == 3:
        assert data_cube_temp.shape[0] == 1, ('First fits file dimension '
                                              'should be the stokes nr.')

    if is_PV:
        assert (hdr['CTYPE2'] == 'VELO-LSR' or
                hdr['CTYPE2'] == 'VRAD' or
                hdr['CTYPE2'] == 'VEL'), ('Second dimension of PV file '
                                          'should be velocity')

    # Remove NaNs
    NaNs = np.isnan(data_cube_temp)
    data_cube_temp[NaNs] = 0.0

    assert len(data_cube_temp.shape) == 4 or len(data_cube_temp.shape) == 3 or len(data_cube_temp.shape) == 2
    if len(data_cube_temp.shape) == 4:
        ni, nz, ny, nx = data_cube_temp.shape
        data_cube = np.zeros((ny, nx), dtype=np.float64)

        # Collapse redundant dimension containing stokes_nr.
        data_cube[:, :] = data_cube_temp[0, 0, :, :]
    if len(data_cube_temp.shape) == 3:
        nz, ny, nx = data_cube_temp.shape
        data_cube = np.zeros((ny, nx), dtype=np.float64)

        # Collapse redundant dimension containing stokes_nr.
        data_cube[:, :] = data_cube_temp[0, :, :]

    elif len(data_cube_temp.shape) == 2:
        ny, nx = data_cube_temp.shape
        data_cube = data_cube_temp

    return data_cube, hdr


PV_file_name = 'pv_test_file.fits'
cube_file = 'im_H13CN_c_cnv.fits'

# Define a path in the image from pv_path_1 to pv_path_2
rotation_axis_len = 50
pv_path_1 = (500, 500 - rotation_axis_len/2.0)  # (x, y) pixel values
pv_path_2 = (500, 500 + rotation_axis_len/2.0)  # (x, y) pixel values

# The path defines the line, whereupon normal
# dimension is collapsed. I.e. the path line defines the offset direction
# vector in the pv diagram.

image_path = Path([pv_path_1, pv_path_2])
data_cube, hdr = load_fits_cube(cube_file)

hdu = extract_pv_slice(data_cube, image_path)

# Then create a HDUList to contain the newly created primary
# HDU, and write to a new file:
hdulist = pyfits.HDUList([hdu])

# Copy original header files onto the new fits file - this is very
# important in order for miriad functions to work properly, i.e. CDELT3
# (channel width) is needed in order to calculate the zeroth moment map
# properly as a channel width of 1 km/s is otherwise assumed.
exclude_keywords_list = ['CTYPE3', 'CRVAL3', 'CDELT3', 'CRPIX3',
                         'CUNIT3', 'CTYPE4', 'CRVAL4', 'CDELT4',
                         'CRPIX4', 'CUNIT4']
prihdr = hdulist[0].header
for x in hdr:
    if x not in exclude_keywords_list:
        try:
            prihdr[x] = hdr[x]
        except ValueError:
            pass

prihdr['NAXIS2'] = hdr['NAXIS3']
prihdr['CTYPE2'] = hdr['CTYPE3']
prihdr['CDELT2'] = hdr['CDELT3']
prihdr['CRVAL2'] = hdr['CRVAL3']
prihdr['CRPIX2'] = hdr['CRPIX3']
hdulist.writeto(PV_file_name, clobber=True)
hdulist.close()

n_x = int(hdr['NAXIS1'])
n_chan = int(hdr['NAXIS2'])

font_size = 12

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

image_data, im_hdr = load_fits_im(PV_file_name)
im = ax.imshow(image_data, aspect='auto', interpolation='nearest',
               cmap=plt.cm.YlOrBr, origin='lower')

ax.set_xlabel('x [pixels]', fontsize=font_size)
ax.set_ylabel(r'y [channels]', fontsize=font_size)

# Colorbar
colorbar_ax = fig.add_axes([0.9, 0.11, 0.05, 0.77])
fig.colorbar(im, cax=colorbar_ax, label=r'Jy Beam$^{-1}$')
colorbar_ax.tick_params(labelsize=10)

plt.show()
