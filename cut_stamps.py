"""
Program for cutting stamps for individual galaxies from a fits mosaic
Input
- Science image
- Weight image already converted into rms image
Output:
- The galaxy stamp themselves
- A text file containing the positions in pixel of those stamps in the total mosaic
By Fernando Buitrago (fbuitrago@gmail.com)
"""

#import montage_wrapper as montage
import numpy as np
import os
import pdb
from astropy.io import fits
from astropy import wcs
from astropy.io import ascii
import montage_wrapper as montage


def ids_as_str(ids):
    #it returns a numpy array of strings, when receiving the ids as float numbers
    ids = np.array(ids, dtype=np.str)
    ids = np.char.rstrip(ids, '.0')
    return ids


path_img = "/mnt/disk1/fb/" 
path_rms = "/mnt/disk1/fb/"
path_cat = "/mnt/disk1/fb/massive_disks/"
name_img = "hlsp_hlf_hst_wfc3-60mas_goodss_f160w_v1.0_sci.fits" 
name_rms = "hlsp_hlf_hst_wfc3-60mas_goodss_f160w_v1.0_wht.fits" 
name_cat = "massive_high_z_bin_goodss.cat"
name_coo = "massive_high_z_bin_goodss.coo"
flag_multiply_by_exptime  = 0
side = 199 #[pix]

#reading the catalog
gal,ra,dec = np.loadtxt(path_cat+name_cat, dtype = float, unpack = True) #OJO: the names are also floats!
gal = ids_as_str(gal)

#creating folders for the stamps and the rms stamps
if not os.path.exists("./galaxy_images"):
    os.mkdir("./galaxy_images")
if not os.path.exists("./sigma_images"):
    os.mkdir("./sigma_images")

#obtaining the galaxy coordinates in pixels
img = fits.open(path_img+name_img)
ww = wcs.WCS(img[0].header)
coo = []
for ii in range(len(ra)):
    coo.append([ra[ii],dec[ii]])
coo = np.array(coo, dtype = float)
pix_coo = ww.wcs_world2pix(coo,1)
for ii in range(len(gal)):
    print(gal[ii],pix_coo[ii,0],pix_coo[ii,1])
    
#calculating the coordinates in case I need to paste the galaxies back in the image
x0 = pix_coo[:,0] - (side/2.)
x1 = pix_coo[:,0] + (side/2.)
y0 = pix_coo[:,1] - (side/2.)
y1 = pix_coo[:,1] + (side/2.)

#Round elements of the array to the nearest integer
x0 = np.rint(x0); x0 = x0.astype(int)
x1 = np.rint(x1); x1 = x1.astype(int)
y0 = np.rint(y0); y0 = y0.astype(int)
y1 = np.rint(y1); y1 = y1.astype(int)

#creating the galaxy stamps
for ii in range(len(gal)):
    montage.mSubimage_pix(path_img+name_img,"./galaxy_images/"+gal[ii]+".fits"    , x0[ii], y0[ii], side)
    montage.mSubimage_pix(path_rms+name_rms,"./sigma_images/" +gal[ii]+"_rms.fits", x0[ii], y0[ii], side)
    
    #in case we are managing the weight image instead of the rms image
    if name_rms.find("wht"):
        new_img = fits.open("./sigma_images/" +gal[ii]+"_rms.fits")
        fits.writeto("./sigma_images/" +gal[ii]+"_rms.fits",1./np.sqrt(new_img[0].data),new_img[0].header,clobber=True)
        
#saving the coordinates where I cut the galaxies from
ascii.write([gal,x0,x1,y0,y1],names=['id','x0','x1','y0','y1'],output=name_coo,Writer = ascii.CommentedHeader)
