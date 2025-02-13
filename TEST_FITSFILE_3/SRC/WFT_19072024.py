import sys
import os
import astropy.io.fits as pyfits
import astropy.units as _u
from astropy import coordinates
import numpy as np
import TicToc

# Получаем имя файла из аргументов командной строки
if len(sys.argv) != 2:
    print("Usage: python WFT_19072024.py <filename>")
    sys.exit(1)

filename = sys.argv[1]

# Убедитесь, что файл существует
if not os.path.exists(filename):
    print(f"Error: {filename} does not exist")
    sys.exit(1)

DIR_IN = 'LOAD_FILE'
DIR_OUT = 'PROCESSED'
fname_input = filename
fname_output = f'{DIR_OUT}/TMP_{os.path.basename(filename)}'

# Чтение заголовка и данных из FITS-файла
header = pyfits.getheader(fname_input, ext=0)
image_data = pyfits.getdata(fname_input, ext=0)

PIXSIZE = header['XPIXSZ']
# FOCALLEN = header['FOCALLEN']
NX = header['NAXIS2']
NY = header['NAXIS1']
WIDTH = NY
HEIGHT = NX

FOCALLEN = 551

# Расчет масштаба и радиуса поиска
SCALE = PIXSIZE / FOCALLEN * _u.rad.to('arcsec') * _u.um.to('mm')
SCALE_LOW = SCALE * _u.arcsec.to('deg') * np.min([NX, NY]) / 1.2
SCALE_HIGH = SCALE * _u.arcsec.to('deg') * np.max([NX, NY]) * 1.2
SEARCH_RAD = SCALE_HIGH

RA = header['OBJCTRA']
DEC = header['OBJCTDEC']
COORDS = coordinates.SkyCoord(RA, DEC, unit=('hour', 'deg'), frame='icrs', equinox='J2000')

print(f'RA = {RA}, DEC = {DEC}')
print(f'RA(deg.) = {COORDS.ra.deg:.5f},  DEC(deg.) = {COORDS.dec.deg:.5f}')

ASTROMETRY_CONFIG = 'CONFIGS/astrometry.cfg'
DOWNSAMPLE = 4

# Формирование команды solve-field с обновленными параметрами
cmd = (
    f"solve-field --config {ASTROMETRY_CONFIG} --overwrite --downsample {DOWNSAMPLE} --cpulimit 3600 --no-plots "
    f"--scale-units arcsecperpix --scale-low {SCALE * 0.9} --scale-high {SCALE * 1.1} "
    f"--x-column X_CENTER --y-column Y_CENTER --sort-column AREA "
    f"--ra {COORDS.ra.deg} --dec {COORDS.dec.deg} --radius {SEARCH_RAD} "
    f"--width {WIDTH} --height {HEIGHT} TMP/XY.fits"
    

)

print(f"Executing command: {cmd}")

try:
    TicToc.tic()
    os.system(cmd)
    TicToc.toc()


    cmd = f'~/anaconda3/envs/radec/bin/new-wcs -i {fname_input} -w TMP/XY.wcs -o TMP/WCS_{os.path.basename(filename)}' #for server assy

    # Обновление WCS в исходном файле
    # cmd = f'~/miniconda3/envs/radec/bin/new-wcs -i {fname_input} -w TMP/XY.wcs -o TMP/WCS_{os.path.basename(filename)}' #for server assy
    os.system(cmd)

    # Сохранение обновленного FITS-файла с новым WCS
    header = pyfits.getheader(f'TMP/WCS_{os.path.basename(filename)}', ext=0)

    fits = pyfits.open(fname_input)
    fits[0].header = header
    fits[0].data = image_data  # сохраним данные изображения
    fits.writeto(fname_output, output_verify='silentfix', overwrite=True)
    fits.close()

except Exception as e:
    print(f'No WCS solution for {fname_input}: {e}')