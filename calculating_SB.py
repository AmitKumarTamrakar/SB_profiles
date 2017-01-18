"""
Program for calculating the surface brightness profile of a given galaxy
Output:
- A plot of the galaxy with the ALL the elliptical annuli applied to get its surface brightness profile
- The surface brightness profile itself
By Fernando Buitrago (fbuitrago@gmail.com)
"""

from astropy.io import fits
import numpy as np
import math
from photutils import EllipticalAperture, EllipticalAnnulus
from photutils import aperture_photometry
from astropy.table import hstack
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pdb #for debugging purposes only
from astropy.cosmology import FlatLambdaCDM


def finding_coo_max_pix(img):
    return np.unravel_index(np.nanargmax(img),img.shape)
    #np.nanargmax->gets the maximum element of the image without taking into account the NaN values
    #np.unravel_index->gets the two dimensional position from a coordinate in a flatten vector


path_images       = '/mnt/disk2/fb/ellipticals_hudf09/galfiting_quadruple_sersic_natural_star_with_new_mask/'
path_sigma_images = '/mnt/disk2/fb/ellipticals_hudf09/sigma_images_h_band/'
pix_scale = 0.06 #[arcsec/pix]
zz = 0.667
step_1 = 0.5        #[kpc]
transition_1_2 = 2. #[kpc]
step_2 = 2.         #[kpc]
x1 = 15.        #[arcsec]
y0 = 32.        #[mag/arcsec^2]
y1 = 16.        #[mag/arcsec^2]
figsize_x = 10. #[inches]
figsize_y = 10. #[inches]
zeropoint = 25.96

#moving distances from kpc to arcsec
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
kpc_per_arcsec = 1./cosmo.arcsec_per_kpc_proper(zz)
step_1         = step_1        /kpc_per_arcsec; step_1         = step_1.value
transition_1_2 = transition_1_2/kpc_per_arcsec; transition_1_2 = transition_1_2.value
step_2         = step_2        /kpc_per_arcsec; step_2         = step_2.value

#reading the images====
img = fits.open(path_images+'model_3_h_conv.fits')
image  = img[0].data
header = img[0].header

rms = fits.open(path_sigma_images+'3_h_rms.fits' )
rms_image  = rms[0].data
rms_header = rms[0].header
#===

#getting the galaxy center
centroid_y,centroid_x= finding_coo_max_pix(image)
centroid = (centroid_x,centroid_y)

#calculating the semi-major axis
x0_1 = step_1
x1_1 = transition_1_2 
first_part_vect_dist  = np.arange(x0_1/pix_scale,x1_1/pix_scale,step_1/pix_scale)
x0_2 = x1_1
x1_2 = x1
second_part_vect_dist = np.arange(x0_2/pix_scale,x1_2/pix_scale,step_2/pix_scale) 
aa = np.concatenate((first_part_vect_dist,second_part_vect_dist))
resolution = aa.size

#defining the apertures parameters
axis_ratio = 0.90
position_angle = math.radians(75.22)+math.radians(90.) #math.radians(90.) for having the same reference system
bb = aa*axis_ratio #semi-minor axis

#creating the figure with the galaxy and its apertures
fig = plt.figure(figsize=(figsize_x,figsize_y))
ax  = plt.subplot() #nothing inside because it is the only plot
plt.imshow(image, cmap = 'cubehelix', norm = colors.LogNorm(), interpolation = 'nearest')

#calculating the flux and the area in each elliptical annulus
apers = []
area  = []
for ii in range(len(aa)):
    """ THIS WORKED WELL AS WELL
    ellip_aper = EllipticalAperture(centroid, aa[ii], bb[ii], position_angle)
    area = ellip_aper.area()
    flux.append(aperture_photometry(image, ellip_aper, error = rms_image))
    """
    if ii == 0:
        ellip_annulus = EllipticalAnnulus(centroid,       0,aa[ii],bb[ii],position_angle)
    else:
        ellip_annulus = EllipticalAnnulus(centroid,aa[ii-1],aa[ii],bb[ii],position_angle)
    apers.append(aperture_photometry(image, ellip_annulus, error = rms_image))
    area_annulus = ellip_annulus.area()
    area.append(area_annulus)
    #I plot the apertures
    ellip_annulus.plot(ax=ax)

#area's rows are scalars, but flux rows are aperture_photometry tables, and I need to join all of them. I get this by doing
table_apers = hstack(apers)

#printing the image with the apertures
plt.savefig("pyapertures.pdf")
plt.close()
plt.clf()

#create two vectors, one with the fluxes per area and the other with the error in the fluxes per area
flux       = np.zeros(resolution)
sigma_flux = np.zeros(resolution)
for ii in range(resolution):

    keyword = 'aperture_sum_'+str(ii+1)
    ff       = float(table_apers[keyword])
    flux[ii] = ff    

    keyword = 'aperture_sum_err_'+str(ii+1)
    ff_sigma = float(table_apers[keyword])
    sigma_flux[ii] = ff_sigma

#calculating the SB profile
surf_bright       =             -2.5*np.log10(  flux            /area )+zeropoint+5.*np.log10(pix_scale)
surf_bright_error = np.absolute(-2.5*np.log10( (flux+sigma_flux)/area )+zeropoint+5.*np.log10(pix_scale)) - surf_bright
distance_in_arcsec = aa*pix_scale

#creating the figure with the surface brightness profile
fig = plt.figure(figsize=(figsize_x,figsize_y))
ax  = plt.subplot() #nothing inside because it is the only plot
ax.set_xlim(0.,x1)
ax.set_ylim(y0,y1)
ax.scatter  (distance_in_arcsec,surf_bright, marker =".", color = "c")
plt.errorbar(distance_in_arcsec,surf_bright, yerr= surf_bright_error, linestyle = "None", ecolor = "r")
plt.xlabel('Radius [arcsec]')
plt.ylabel('Surface brightness [mag arcsec$^{-2}$]')
plt.grid()
ax2 = ax.twiny() #they share the same y-axis
ax2.plot(distance_in_arcsec*kpc_per_arcsec.value,surf_bright, alpha = 0)
plt.xlabel('Radius [kpc]')
plt.savefig("pyapertures_surfbright.pdf")
plt.close()
plt.clf()
