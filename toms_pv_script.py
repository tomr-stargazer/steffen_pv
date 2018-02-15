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

from pvextractor import extract_pv_slice, Path
from pv_test_script import save_extracted_pv_slice_hdu, load_fits_im

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
def PV_diagram_plot(filename, label, length=4*u.arcsec):

    cube = spectral_cube.SpectralCube.read(filename)

    # let's make the "parallel" PV diagram.

    # It's centered on whatever pixel is at the B source center
    B_ra_px, B_dec_px, zero_vel_channel = cube.wcs.all_world2pix(source_B.fk5.ra, source_B.fk5.dec, 0, 0)
    source_position_px = (B_ra_px, B_dec_px)
    cen_px = source_position_px

    # and it is 4'' long.
    semi_length_in_arcsec = length/2
    semi_length_in_px = semi_length_in_arcsec.to(u.deg)/(cube.header['CDELT2']*u.Unit(cube.header['CUNIT2']))

    # its direction is along the binary axis
    delta_x_px = semi_length_in_px * np.cos(AB_position_angle)
    delta_y_px = semi_length_in_px * np.sin(AB_position_angle)

    # Define a path in the image from pv_path_1 to pv_path_2
    pv_path_1 = (cen_px[0] + delta_x_px, cen_px[1] + delta_y_px)  # (x, y) pixel values
    pv_path_2 = (cen_px[0] - delta_x_px, cen_px[1] - delta_y_px)  # (x, y) pixel values

    # The path defines the line, whereupon normal
    # dimension is collapsed. I.e. the path line defines the offset direction
    # vector in the pv diagram.

    image_path = Path([pv_path_1, pv_path_2], width=3)

    hdu = extract_pv_slice(cube.hdu.data, image_path)
    PV_file_name = 'pv_test_file.fits'
    save_extracted_pv_slice_hdu(hdu, cube.hdu.header, PV_file_name)
    image_data, im_hdr = load_fits_im(PV_file_name)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    im = ax.imshow(image_data, aspect='auto', interpolation='nearest',
                   cmap=plt.cm.YlOrBr, origin='lower')
    ax.set_xlabel('x [pixels]', fontsize=12)
    ax.set_ylabel(r'y [channels]', fontsize=12)
    ax.set_title("{0} PV diagram".format(label))

    # Colorbar
    colorbar_ax = fig.add_axes([0.9, 0.11, 0.05, 0.77])
    fig.colorbar(im, cax=colorbar_ax, label=r'Jy Beam$^{-1}$')
    colorbar_ax.tick_params(labelsize=10)

    return fig

    # plt.show()


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

    if False:
        contpeak_spectra(filepath, label)

        redblue_moments(filepath, label)
    PV_diagram_plot(filepath, label)



plt.show()


