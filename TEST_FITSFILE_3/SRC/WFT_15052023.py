#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  3 08:30:43 2021

@author: alex

"""
import astropy.io.fits as pyfits
import astropy.units as _u
from os import system, listdir
import numpy as np
from astropy import coordinates
import TicToc

DIR_IN = 'LOAD_FILE'
DIR_OUT = 'PROCESS_FILE'
fname = listdir(DIR_IN)[0]

fname_input = f'{DIR_IN}/{fname}'
fname_output = f'{DIR_OUT}/TMP_{fname}'

header = pyfits.getheader(fname_input, ext=0)
image_data = pyfits.getdata(fname_input, ext=0)

PIXSIZE = header['XPIXSZ']
FOCALLEN = header['FOCALLEN']
NX = header['NAXIS2']
NY = header['NAXIS1']

SCALE = PIXSIZE/FOCALLEN*_u.rad.to('arcsec')*_u.um.to('mm')

fits = pyfits.open(fname_input)
fits[0].header = header   
fits[0].data = image_data 
fits.writeto(fname_output,output_verify='silentfix',overwrite=True)
fits.close()

RA = header['OBJCTRA']
DEC = header['OBJCTDEC']

COORDS = coordinates.SkyCoord(RA,DEC,unit=('hour','deg'), frame='icrs', equinox='J2000')

print(f'RA = {RA}, DEC = {DEC}')
print(f'RA(deg.) = {COORDS.ra.deg},  DEC(deg.) = {COORDS.dec.deg}')

ASTROMETRY_CONFIG = 'astrometry.cfg'
DOWNSAMPLE = 4
SCALE_LOW = SCALE*_u.arcsec.to('deg')*np.min([NX,NY])/1.2
SCALE_HIGH = SCALE*_u.arcsec.to('deg')*np.max([NX,NY])*1.2
SEARCH_RAD = SCALE_HIGH 

####  this is based on field of view
#str1 = f'solve-field --config {ASTROMETRY_CONFIG} --overwrite --downsample {DOWNSAMPLE} --cpulimit 3600 --no-plots '
#str2 = f'--scale-units degwidth --scale-low {SCALE_LOW} --scale-high {SCALE_HIGH} ' 
#str3 = f'--ra {COORDS.ra.deg} --dec {COORDS.dec.deg} --radius {SEARCH_RAD} {fname_output}'

####  this is based on scale (arcsec per pixel)
str1 = f'solve-field --config {ASTROMETRY_CONFIG} --overwrite --downsample {DOWNSAMPLE} --cpulimit 3600 --no-plots '
str2 = f'--scale-units arcsecperpix --scale-low {SCALE/1.1} --scale-high {SCALE*1.1} '
str3 = f'--ra {COORDS.ra.deg} --dec {COORDS.dec.deg} --radius {SEARCH_RAD} {fname_output}'

cmd = str1 + str2 + str3

print(cmd)

TicToc.tic()
system(cmd)
TicToc.toc()

#fn1 = fname_output.split('.')[0]+'.new'
#fn2 = 'WCS' + fname_output[3:]

#print('Rewrite {0} to {1}'.format(fn1, fn2))

#cmd = ('mv {0} {1}'.format(fn1, fn2))
#system(cmd)
