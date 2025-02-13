import astropy.io.fits as pyfits
import astropy.units as _u
from os import system
from astropy import coordinates
import TicToc
from os import listdir
import numpy as np

DIR_IN = 'LOAD_FILE'
DIR_OUT = 'PROCESSED'
fname = listdir(DIR_IN)[0]

fname_input = f'{DIR_IN}/{fname}'
fname_output = f'{DIR_OUT}/TMP_{fname}'

# Чтение заголовка и данных из FITS-файла
header = pyfits.getheader(fname_input, ext=0)
image_data = pyfits.getdata(fname_input, ext=0)

PIXSIZE = header['XPIXSZ']
WIDTH = header['NAXIS1']
HEIGHT = header['NAXIS2']
# FOCALLEN = 551 #header['FOCALLEN']
BINNING = header['XBINNING']

RA = header['OBJCTRA']
DEC = header['OBJCTDEC']

COORDS = coordinates.SkyCoord(RA, DEC, unit=('hour', 'deg'), frame='icrs', equinox='J2000')

print(f'RA = {RA}, DEC = {DEC}')
print(f'RA(deg.) = {COORDS.ra.deg:.5f},  DEC(deg.) = {COORDS.dec.deg:.5f}')

ASTROMETRY_CONFIG = 'CONFIGS/astrometry.cfg'
DOWNSAMPLE = 4

# Выполнение новой команды solve-field
cmd = (
    "solve-field --overwrite --scale-units arcsecperpix --scale-low 0.1 --scale-high 10 --no-plots "
    f"--x-column X_CENTER --y-column Y_CENTER --sort-column AREA "
    f"--width {WIDTH} --height {HEIGHT} TMP/XY.fits"
)

print(f"Executing command: {cmd}")

try:
    TicToc.tic()
    system(cmd)
    TicToc.toc()

    # Обновление WCS в исходном файле
    cmd = f'/usr/local/astrometry/bin/new-wcs -i {fname_input} -w TMP/XY.wcs -o TMP/WCS_{fname}'
    system(cmd)

    # Сохранение обновленного FITS-файла с новым WCS
    header = pyfits.getheader(f'TMP/WCS_{fname}', ext=0)

    fits = pyfits.open(fname_input)
    fits[0].header = header
    fits[0].data = image_data  # сохраним данные изображения
    fits.writeto(fname_output, output_verify='silentfix', overwrite=True)
    fits.close()

except Exception as e:
    print(f'No WCS solution for {fname_input}: {e}')