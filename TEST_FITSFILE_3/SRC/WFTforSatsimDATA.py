#!/usr/bin/env python3
import sys
import os
import re
import astropy.io.fits as pyfits
import TicToc

# Директории
TMP_DIR = "TMP"
OUTPUT_DIR = "PROCESSED"

# Проверка аргументов
if len(sys.argv) != 2:
    print("Usage: python astrometry_script.py <filename>")
    sys.exit(1)

fname_input = sys.argv[1]

if not os.path.exists(fname_input):
    print(f"Error: {fname_input} does not exist")
    sys.exit(1)

# Определяем базовое имя без расширения
base_name = os.path.splitext(os.path.basename(fname_input))[0]

# Итоговый FITS-файл будет сохранен в OUTPUT_DIR с именем, равным базовому имени
fname_output = os.path.join(OUTPUT_DIR, f"{base_name}.fits")

# Создаем необходимые директории, если их нет
for d in [TMP_DIR, OUTPUT_DIR]:
    os.makedirs(d, exist_ok=True)

# Опции для astrometry.net
ASTROMETRY_CONFIG = "/home/hellnim/radecWork/TEST_FITSFILE_3/CONFIGS/astrometry.cfg"
DOWNSAMPLE = 4

# Формирование команды astrometry.net.
# Здесь мы используем опцию --dir TMP_DIR, чтобы все выходные файлы сохранялись в TMP,
# и опцию -o, задающую базовое имя, равное base_name.
cmd = (
    f"solve-field --config {ASTROMETRY_CONFIG} --continue --downsample {DOWNSAMPLE} --no-plots "
    f"--dir {TMP_DIR} -o {base_name} {fname_input}"
)

print("Executing command:")
print(cmd)

try:
    TicToc.tic()
    os.system(cmd)
    TicToc.toc()

    # Проверяем, создался ли WCS-файл (ожидается TMP/<base_name>.wcs)
    wcs_file = os.path.join(TMP_DIR, f"{base_name}.wcs")
    if not os.path.exists(wcs_file):
        raise FileNotFoundError(f"WCS file not found: {wcs_file}")

    # Команда new-wcs: используем входной файл, созданный WCS-файл, и записываем новый файл с WCS
    cmd_newwcs = (
        f"~/anaconda3/envs/radec/bin/new-wcs -i {fname_input} -w {wcs_file} -o {TMP_DIR}/WCS_{base_name}.fits"
    )
    print("Executing new-wcs command:")
    print(cmd_newwcs)
    os.system(cmd_newwcs)

    # Проверяем, создался ли обновленный FITS-файл с WCS
    updated_fits = os.path.join(TMP_DIR, f"WCS_{base_name}.fits")
    if not os.path.exists(updated_fits):
        raise FileNotFoundError(f"Updated FITS file not found: {updated_fits}")

    # Считываем новый header и данные из исходного FITS-файла
    header_new = pyfits.getheader(updated_fits, ext=0)
    image_data = pyfits.getdata(fname_input, ext=0)

    # Корректируем DATE-OBS в новом заголовке, удаляя часовой пояс, если он присутствует
    if 'DATE-OBS' in header_new:
        date_obs = header_new['DATE-OBS']
        m = re.match(r"(.+T\d{2}:\d{2}:\d{2})", date_obs)
        if m:
            header_new['DATE-OBS'] = m.group(1)
            print("Updated DATE-OBS in new header:", header_new['DATE-OBS'])

    # Обновляем исходный FITS-файл новым header и сохраняем в OUTPUT_DIR под именем base_name.fits
    fits_obj = pyfits.open(fname_input)
    fits_obj[0].header = header_new
    fits_obj[0].data = image_data
    fits_obj.writeto(fname_output, output_verify='silentfix', overwrite=True)
    fits_obj.close()

    print(f"Updated FITS file saved as {fname_output}")

except Exception as e:
    print(f'No WCS solution for {fname_input}: {e}')