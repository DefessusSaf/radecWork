#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 14:37:37 2024

@author: alex

Rewrite results after SExtractor for region plot on DS9, etc.

Input information format:
    X_IMAGE         Object position along x                     [pixel]
    Y_IMAGE         Object position along y                     [pixel]
    ERRCXX_IMAGE    Cxx error ellipse parameter                 [pixel**(-2)]
    ERRCYY_IMAGE    Cyy error ellipse parameter                 [pixel**(-2)]
    A_IMAGE         Profile RMS along major axis                [pixel]
    B_IMAGE         Profile RMS along minor axis                [pixel]
    XMIN_IMAGE      Minimum x-coordinate among detected pixels  [pixel]
    YMIN_IMAGE      Minimum y-coordinate among detected pixels  [pixel]
    XMAX_IMAGE      Maximum x-coordinate among detected pixels  [pixel]
    YMAX_IMAGE      Maximum y-coordinate among detected pixels  [pixel]
    THETA_IMAGE     Position angle (CCW/x)                      [deg]
    FLAGS           Extraction flags                                         
    FLUX_ISOCOR     Corrected isophotal flux                   [count]

output format:
    ellipse(   269.4581,   397.0928,   46.4340,   1.7650,  -69.43)
    
Запускать под оболочкой astroconda (есть конфликт в numpy версиях  
1.21.5 в astropy vs. 1.26.4 в seiferts etc.

"""

import numpy as np
from os import system

#%%
DIR = 'TMP/'
f_res = f'{DIR}k1-impTEST.fts.sx'
f_out = f'{DIR}elik1-imp-field-TEST.reg'

X,Y,ERRX,ERRY,A,B,XMIN,YMIN,XMAX,YMAX,TH,FLAG,FLUX= np.genfromtxt(f_res, unpack=True)

output = []

for x,y,a,b,th in zip(X,Y,A,B,TH):
    output.append(f'ellipse( {x} , {y} , {a:.4f} , {b:.4f},{th:.2f}')

np.savetxt(f_out, output, fmt='%s')

#%%

fn = 'XY.txt'
arr = np.vstack((X, Y, FLUX))
np.savetxt(f'{DIR}{fn}', arr.T, fmt=('%s %s %s'), header='X_CENTER  Y_CENTER  AREA')

## use  text2fits  utility file from astrometry.net

cmd = f'~/anaconda3/envs/radec/bin/text2fits {DIR}{fn}  {DIR}{fn.split(".")[0]}.fits' #for VM PC on work
# cmd = f'~/miniconda3/envs/radec/bin/text2fits {DIR}{fn}  {DIR}{fn.split(".")[0]}.fits' #for server assy
system(cmd) 