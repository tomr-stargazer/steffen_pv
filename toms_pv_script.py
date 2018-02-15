"""
Tom's shot at simplifying the PV approach.

"""

import os
import numpy as np
import matplotlib.pyplot as plt
import spectral_cube
from astropy.coordinates import SkyCoord
import astropy.units as u
import astropy.stats
import aplpy

# These coords were read by-eye out of CASAVIEWER.
source_A = SkyCoord('03h27m39.109s', '+30d13m03.019s', frame='fk5')
source_B = SkyCoord('03h27m39.119s', '+30d13m02.800s', frame='fk5')

AB_position_angle = source_A.position_angle(source_B)

print("Position angle: {:.2f}".format(AB_position_angle.to(u.deg)) )

A_ra = source_A.fk5.ra
A_dec = source_A.fk5.dec

B_ra = source_B.fk5.ra
B_dec = source_B.fk5.dec

files_list = ["subcube_spw25_20kms_10arcsec.fits",
    "subcube_spw27_20kms_10arcsec.fits",
    "subcube_spw29_20kms_10arcsec.fits",
    "subcube_spw31_20kms_10arcsec.fits",
    "subcube_spw33_20kms_10arcsec.fits",
    "subcube_spw37_20kms_10arcsec.fits",
    "subcube_spw39_20kms_10arcsec.fits",]

labels_list = ["HC15N", "H13CN", "SiO", "H13CO+", "SO", "CH3OH", "CS"]




# prototype land
if False:

    filename = files_list[1]
    label = labels_list[1]

    cube = spectral_cube.SpectralCube.read(filename)

    # let's make the "parallel" PV diagram.

    # It's centered on whatever pixel is at the B source center
    B_ra_px, B_dec_px, zero_vel_channel = cube.wcs.all_world2pix(source_B.fk5.ra, source_B.fk5.dec, 0, 0)


    # its direction is along the binary axis

    # and it is 4'' long.


    # position of the H13CN peak emission, in pixel coordinates
    source_position_px = (302, 360)

    # Define a path in the image from pv_path_1 to pv_path_2
    rotation_axis_len = 50
    pv_path_1 = (source_position_px[0], source_position_px[1] - rotation_axis_len/2.0)  # (x, y) pixel values
    pv_path_2 = (source_position_px[0], source_position_px[1] + rotation_axis_len/2.0)  # (x, y) pixel values

    # The path defines the line, whereupon normal
    # dimension is collapsed. I.e. the path line defines the offset direction
    # vector in the pv diagram.

    image_path = Path([pv_path_1, pv_path_2])
    data_cube, hdr = load_fits_cube(cube_file)

    hdu = extract_pv_slice(data_cube, image_path)



def contpeak_spectra(filename, label):
    """
    Plots the single-pixel spectrum from each of the two continuum peaks.

    """

    cube = spectral_cube.SpectralCube.read(filename)

    A_ra_px, A_dec_px, zero_vel_channel = cube.wcs.all_world2pix(source_A.fk5.ra, source_A.fk5.dec, 0, 0)

    source_A_spectrum = cube[:, int(A_dec_px), int(A_ra_px)]

    B_ra_px, B_dec_px, zero_vel_channel = cube.wcs.all_world2pix(source_B.fk5.ra, source_B.fk5.dec, 0, 0)

    source_B_spectrum = cube[:, int(B_dec_px), int(B_ra_px)]


    fig = plt.figure()
    plt.plot(source_A_spectrum.spectral_axis, source_A_spectrum, label='Continuum peak A')
    plt.plot(source_B_spectrum.spectral_axis, source_B_spectrum, label='Continuum peak B')
    plt.title(label)
    plt.legend()
    plt.savefig("continuum_peak_spectra/{0}_cont_peaks_spectra.png".format(label), bbox_inches='tight')

    return fig


def redblue_moments(filename, label):
    """
    Contour plots for the redshifted/blueshifted halves of emission.

    """

    # moment 0 with red/blue
    central_velocity = 5*u.km/u.s
    red_range = (central_velocity, central_velocity+10*u.km/u.s)
    blue_range = (central_velocity-10*u.km/u.s, central_velocity)

    cube = spectral_cube.SpectralCube.read(filename)

    blue_slab = cube.spectral_slab(*blue_range)
    red_slab = cube.spectral_slab(*red_range)

    cube_mom0 = cube.moment0()
    blue_mom0 = blue_slab.moment0()
    red_mom0 = red_slab.moment0()

    blue_sigma = astropy.stats.mad_std(blue_mom0).value
    red_sigma = astropy.stats.mad_std(red_mom0).value
    mean_sigma = (blue_sigma+red_sigma)/2

    # Plot some contours and also
    fig = aplpy.FITSFigure(cube_mom0.hdu)

    fig.show_grayscale(invert=True, stretch='log', vmax=cube_mom0.max().value/2, vmid=cube_mom0.min().value, interpolation='nearest')
    fig.show_contour(data=blue_mom0.hdu, colors='b', levels=[3*mean_sigma, 9*mean_sigma, 15*mean_sigma])
    fig.show_contour(data=red_mom0.hdu, colors='r', levels=[3*mean_sigma, 9*mean_sigma, 15*mean_sigma])
    fig.set_title(r"$\rm{{ {0} }}$ Moment 0 map of L1455-IRS1 from ALMA".format(label))
    fig.add_beam(color='black', edgecolor='w')

    fig.show_markers([source_A.ra.value, source_B.ra.value], [source_A.dec.value, source_B.dec.value], edgecolor='w')

    # center_in_world_coords = fig.pixel2world(*(x/2 for x in cube_mom0.shape))
    fig.recenter(source_A.ra, source_A.dec, radius=2.5/3600)

    plt.savefig("redblue_moments/{0}_redblue_moments.png".format(label), bbox_inches='tight')

    return fig


for file, label in zip(files_list, labels_list):

    filepath = os.path.expanduser("~/ALMA_subcubes/")+file

    if True:
        contpeak_spectra(filepath, label)

    redblue_moments(filepath, label)


plt.show()


