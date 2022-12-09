
import aplpy
import astropy
from astropy.io import fits
import numpy as np
import pyfits
from os import sys
from astropy.wcs import WCS
import matplotlib.pyplot as plt



fits_image1=sys.argv[1]
fits_image2=sys.argv[2]
fits_image3=sys.argv[3]
fits_image4=sys.argv[4]
fits_image5=sys.argv[5]


fits_images=[fits_image1,fits_image2,fits_image3,fits_image4,fits_image5]

def correct_axis(fits_image):
    
    data_hdu = fits.open(fits_image)[0]
    data = data_hdu.data
    header = data_hdu.header
    data = data[0,0,:,:]
    header = header.copy()
    header['NAXIS'] = 	2
    del header['*3']
    del header['*4']
    new_fits_image=str(fits_image).replace('image.fits','deleted_axis.fits') 
        
    fits.writeto(new_fits_image, header=header, data=data, overwrite=True)
    
    plot_image(new_fits_image)
    




def plot_image(fits_image):

    f= aplpy.FITSFigure(fits_image) 

    #f.show_colorscale(cmap='inferno',vmin =0,vmax=0.00005) 
    f.show_colorscale(cmap='gist_heat',vmin =0,vmax=0.00006) # cutouts
    f.recenter(cluster_centre_RA,cluster_centre_DEC, width =1.2,height = 1.2) ## Specify size and co-ordiante of specific region	for cutouts
    f.add_colorbar()
    f.add_beam()
    f.beam.set_major(header['BMAJ'])  # degrees
    f.beam.set_minor(header['BMIN'])
    #f.show_regions(regions)
    f.beam.set_color("white")
    f.beam.set_hatch('+')
    f.beam.set_edgecolor('yellow')
    f.colorbar.set_font(size=20)
    f.tick_labels.set_font(size=20)
    f.axis_labels.set_font(size=20)
    f.axis_labels.set_xtext("RA (J2000)")
    f.axis_labels.set_ytext("Dec (J2000)")
    f.save(str(fits_image).replace('.fits','.pdf'))
        

if __name__ == "__main__":

    for i in range(len(fits_images)):

        radiohdu = pyfits.open(fits_images[i])[0]
        data = radiohdu.data
        header = radiohdu.header
        cluster_centre_RA=header['CRVAL1']-0.2
        cluster_centre_DEC=header['CRVAL2']+0.05
        if header['NAXIS'] >2:

     
            correct_axis(fits_images[i])
            
        
        else:

            
            plot_image(fits_image)	
     	
