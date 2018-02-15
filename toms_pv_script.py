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


# moment 0 with red/blue
central_velocity = 5*u.km/u.s
red_range = (central_velocity, central_velocity+10*u.km/u.s)
blue_range = (central_velocity-10*u.km/u.s, central_velocity)

filename = files_list[1]
label = labels_list[1]
cube = spectral_cube.SpectralCube.read('/Users/tsrice/ALMA_subcubes/'+filename)

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

# plt.imshow(cube_mom0.data, origin='lower', cmap='Greys')
# plt.contour(blue_mom0.data, color='blue')
fig.show_grayscale(invert=True, stretch='log', vmax=cube_mom0.max().value/2, vmid=cube_mom0.min().value, interpolation='nearest')
fig.show_contour(data=blue_mom0.hdu, colors='b', levels=[3*mean_sigma, 9*mean_sigma, 15*mean_sigma])
fig.show_contour(data=red_mom0.hdu, colors='r', levels=[3*mean_sigma, 9*mean_sigma, 15*mean_sigma])
fig.set_title(r"$\rm{{ {0} }}$ Moment 0 map of L1455-IRS1 from ALMA".format(label))
fig.add_beam(color='black', edgecolor='w')

fig.show_markers([source_A.ra.value, source_B.ra.value], [source_A.dec.value, source_B.dec.value], edgecolor='w')

# center_in_world_coords = fig.pixel2world(*(x/2 for x in cube_mom0.shape))
fig.recenter(source_A.ra, source_A.dec, radius=2.5/3600)

plt.show()


for file, label in zip(files_list, labels_list):

    if True:
        break

    filepath = os.path.expanduser("~/ALMA_subcubes/")+file

    contpeak_spectra(filepath, label)




