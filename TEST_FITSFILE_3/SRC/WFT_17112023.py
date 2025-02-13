#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  3 08:30:43 2021

@author: alex

"""
import astropy.io.fits as pyfits
import astropy.units as _u
from os import system
from astropy import coordinates
import TicToc
from os import listdir, getcwd

#cwd = getcwd()
#files = [f for f in listdir(cwd) if f.startswith('cdobj') and f.endswith('fit')]
#files.sort()

#for fname_input in files:

fname_input = 'obj_20231021201704.fit'
fname_output = 'TMP_'+fname_input
    
header = pyfits.getheader(fname_input, ext=0)
image_data = pyfits.getdata(fname_input, ext=0)
    
PIXSIZE = header['XPIXSZ']
WIDTH = header['NAXIS1']
HEIGHT = header['NAXIS2']
FOCALLEN = 551 #header['FOCALLEN']
BINNING = header['XBINNING']
       
SCALE = BINNING * PIXSIZE/FOCALLEN*_u.rad.to('arcsec')*_u.um.to('mm')
    
RA = header['OBJCTRA']
DEC = header['OBJCTDEC']
    
COORDS = coordinates.SkyCoord(RA,DEC,unit=('hour','deg'), frame='icrs', equinox='J2000')
    
print(f'RA = {RA}, DEC = {DEC}')
print(f'RA(deg.) = {COORDS.ra.deg:.5f},  DEC(deg.) = {COORDS.dec.deg:.5f}')
    
str1 = 'solve-field --config astrometry.cfg --overwrite --downsample 4 --cpulimit 3600 --no-plots '
str2 = f'--scale-units app --scale-low {SCALE*0.95:.2f} --scale-high {SCALE*1.05:.2f} '
str3 = f'--ra {COORDS.ra.deg:.5f} --dec {COORDS.dec.deg:.5f} --radius 5.0 '
str4 = f'--x-column X_CENTER --y-column Y_CENTER --sort-column AREA --width {WIDTH} --height {HEIGHT} '
    
cmd = str1 + str2 + str3 + str4 + 'XY.fits'
    
print(cmd)
    
try:
   TicToc.tic()
   system(cmd)
   TicToc.toc()
        
   cmd = '/usr/local/astrometry/bin/new-wcs -i ' + fname_input + ' -w ' + 'XY.wcs -o WCS' + fname_output
   system(cmd)
        
   header = pyfits.getheader('WCS'+fname_output, ext=0)
        
   fits = pyfits.open(fname_input)
   fits[0].header = header
   fits[0].data = image_data 
   fits.writeto(fname_output,output_verify='silentfix',overwrite=True)
   fits.close()
except:
   print(f"No WCS solution for {fname_input}")    
