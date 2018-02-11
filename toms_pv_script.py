"""
Tom's shot at simplifying the PV approach.

"""

import os
import numpy as np
import matplotlib.pyplot as plt
import spectral_cube
from astropy.coordinates import SkyCoord
import astropy.units as u

# These coords were read by-eye out of CASAVIEWER.
source_A = SkyCoord('03h27m39.109s', '+30d13m03.019s', frame='fk5')
source_B = SkyCoord('03h27m39.119s', '+30d13m02.800s', frame='fk5')

AB_position_angle = source_A.position_angle(source_B)

print("Position angle: {:.2f}".format(AB_position_angle.to(u.deg)) )

# cube_file = os.path.expanduser("~/ALMA_subcubes/subcube_spw27_20kms_10arcsec.fits")
# cube = spectral_cube.SpectralCube.read(cube_file)

A_ra = source_A.fk5.ra
A_dec = source_A.fk5.dec

# A_ra_px, A_dec_px, zero_vel_channel = cube.wcs.all_world2pix(source_A.fk5.ra, source_A.fk5.dec, 0, 0)

# source_A_spectrum = cube[:, int(A_dec_px), int(A_ra_px)]

B_ra = source_B.fk5.ra
B_dec = source_B.fk5.dec

# B_ra_px, B_dec_px, zero_vel_channel = cube.wcs.all_world2pix(source_B.fk5.ra, source_B.fk5.dec, 0, 0)

# source_B_spectrum = cube[:, int(B_dec_px), int(B_ra_px)]


# plt.figure()
# plt.plot(source_A_spectrum.spectral_axis, source_A_spectrum, label='Continuum peak A')
# plt.plot(source_B_spectrum.spectral_axis, source_B_spectrum, label='Continuum peak B')
# plt.legend()
# plt.show()
# plt.savefig("H13CN_cont_peaks_spectra.png", bbox_inches='tight')


files_list = ["subcube_spw25_20kms_10arcsec.fits",
    "subcube_spw27_20kms_10arcsec.fits",
    "subcube_spw29_20kms_10arcsec.fits",
    "subcube_spw31_20kms_10arcsec.fits",
    "subcube_spw33_20kms_10arcsec.fits",
    "subcube_spw37_20kms_10arcsec.fits",
    "subcube_spw39_20kms_10arcsec.fits",]

labels_list = ["HC15N", "H13CN", "SiO", "H13CO+", "SO", "CH3OH?", "CH3OH", "CS"]

def contpeak_spectra(filename, label):

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
    plt.show()
    plt.savefig("continuum_peak_spectra/{0}_cont_peaks_spectra.png".format(label), bbox_inches='tight')

    return fig


for file, label in zip(files_list, labels_list):

    filepath = os.path.expanduser("~/ALMA_subcubes/")+file

    contpeak_spectra(filepath, label)


