import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib import pyplot as plt
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
from astropy.coordinates import Angle
from astropy.time import Time
from datetime import timedelta
import astropy.units as u
import os
import glob

def is_outlier(points, thresh=3.5):
    if len(points.shape) == 1:
        points = points[:, None]
    median = np.median(points, axis=0)
    diff = np.sqrt(np.sum((points - median) ** 2, axis=-1))
    med_abs_deviation = np.median(diff)
    modified_z_score = 0.6745 * diff / med_abs_deviation

    return np.where(modified_z_score > thresh)[0]


def load_data(file_path):
    data =  np.genfromtxt(file_path, unpack=True)
    return data


def preprocess_data(X,Y,ERRX,ERRY,A,B,XMIN,YMIN,XMAX,YMAX,TH,FLAG,FLUX , x_min, x_max, y_min, y_max):
    # y_threshold = np.percentile(Y, percentile)
    # y_mask = Y < y_threshold

    y_mask = (Y >= y_min) & (Y <= y_max)

    x_mask = (X >= x_min) & (X <= x_max)
    mask = x_mask & y_mask

    X,Y,ERRX,ERRY,A,B,XMIN,YMIN,XMAX,YMAX,TH,FLAG,FLUX  = X[mask], Y[mask], ERRX[mask], ERRY[mask], A[mask], B[mask], XMIN[mask], YMIN[mask], XMAX[mask], YMAX[mask], TH[mask], FLAG[mask], FLUX[mask]

    return X,Y,ERRX,ERRY,A,B,XMIN,YMIN,XMAX,YMAX,TH,FLAG,FLUX


def compute_elongation(A, B):
    return A / B

def scale_features(features):
    scaler = StandardScaler()
    return scaler.fit_transform(features)


def ab_ratio(ELONG, threshold):
    return ELONG < threshold


def cluster_data(features_scaled, eps, min_samples):
    dbscan = DBSCAN(eps=eps, min_samples=min_samples)
    return dbscan.fit_predict(features_scaled)


def compute_errores(ERRX, ERRY):
    erroreX = np.sqrt(1/ERRX)
    erroreY = np.sqrt(1/ERRY)
    return erroreX, erroreY


def choose_fits_file():
    log_file = 'TMP/processing_log.txt'
    if not os.path.exists(log_file):
        raise FileNotFoundError("Log file does not exist.")

    with open(log_file, 'r') as f:
        last_file = f.readline().strip()

    if not last_file:
        raise ValueError("No file name found in the log file.")

    last_file_base = os.path.basename(last_file)
    # Извлекаем суффикс (все символы до расширения)
    last_file_suffix = last_file_base.split('.')[0]
    print(f"Extracted suffix: {last_file_suffix}")

    fits_files = glob.glob('TMP/*.fits') + glob.glob('TMP/*.fit')

    # Печать всех найденных файлов
    print("Found FITS files:")
    for f in fits_files:
        print(f)

    matching_file = None
    for fits_file in fits_files:
        file_name = os.path.basename(fits_file)
        # Печать для отладки
        print(f"Checking file: {file_name}")
        # Проверка наличия суффикса в имени файла
        if last_file_suffix in file_name:
            matching_file = fits_file
            break

    if not matching_file:
        raise FileNotFoundError(f"No matching FITS file found for suffix {last_file_suffix}")

    print(f"Selected FITS file: {matching_file}")

    return matching_file, last_file_base


def convert(ra_deg, dec_deg):
    ra_angle = Angle(ra_deg, unit=u.deg)
    ra_hms = ra_angle.to_string(unit=u.hour, sep=':', precision=2)

    dec_angle = Angle(dec_deg, unit=u.deg)
    dec_dms = dec_angle.to_string(unit=u.deg, sep=':', precision=2, alwayssign=True)  # Добавляем alwayssign=True

    return ra_hms, dec_dms


def save_results(coords_first, coords_second, base_filename, fits_filename, x_y_a_b_values, errors):
    """
    Сохраняет результаты кластеризаций вместе с данными X, Y, A, B, XMIN, YMIN, XMAX, YMAX

    Параметры:
        coords_first (list): Координаты RA и DEC первой кластеризации.
        coords_second (list): Координаты RA и DEC второй кластеризации.
        base_filename (str): Базовое имя файла.
        fits_filename (str): Имя FITS файла для чтения заголовка.
        x_y_a_b_values (dict): Словарь со значениями X, Y, A, B, XMIN, YMIN, XMAX, YMAX.
        errors (dict): Словарь с ошибками по X и Y.
    """
    output_dir = 'PROCESS_FILE'
    os.makedirs(output_dir, exist_ok=True)

    txt_filename = os.path.join(output_dir, f'{base_filename}.txt')
    with open(txt_filename, 'w') as f:
        f.write(f"File: {base_filename}\n")

        with fits.open(fits_filename) as hdul:
            header = hdul[0].header
            date_obs = header.get('DATE-OBS', '00000')
            exptime = header.get('EXPTIME', 0)
            # Вычисление среднего времени экспозиции, если EXPTIME задано
            if date_obs != '00000' and exptime > 0:
                # Извлекаем дату и время
                date_str, time_str = date_obs.split('T')
                # Преобразуем строку времени в объект времени
                time_obs = Time(f'{date_str} {time_str}', format='iso')
                # Добавляем половину EXPTIME (в секундах)
                avg_exposure_time = time_obs + timedelta(seconds=(exptime / 2.0))
                # Формируем новую строку DATE-OBS с сохранением даты
                new_time_str = avg_exposure_time.iso.replace(" ", "T")
                f.write(f"{new_time_str}\n")
            else:
                f.write(f"{date_obs}\n")

        # Сохранение первой кластеризации
        for (ra_hms, dec_dms), x, y, a, b, xmin, ymin, xmax, ymax, erroreX, erroreY in zip(
            coords_first,
            x_y_a_b_values['X_first'],
            x_y_a_b_values['Y_first'],
            x_y_a_b_values['A_first'],
            x_y_a_b_values['B_first'],
            x_y_a_b_values['XMIN_first'],
            x_y_a_b_values['YMIN_first'],
            x_y_a_b_values['XMAX_first'],
            x_y_a_b_values['YMAX_first'],
            errors["err_x_first"],
            errors["err_y_first"]
        ):
            f.write(f"{ra_hms} {dec_dms} {x} {y} {erroreX} {erroreY} {a} {b} {xmin} {ymin} {xmax} {ymax}\n")

        # Если есть вторая кластеризация, сохраняем её
        if coords_second:
            f.write(f"#Second cluster:\n")
            for (ra_hms, dec_dms), x, y, a, b, xmin, ymin, xmax, ymax, erroreX, erroreY in zip(
                coords_second,
                x_y_a_b_values['X_second'],
                x_y_a_b_values['Y_second'],
                x_y_a_b_values['A_second'],
                x_y_a_b_values['B_second'],
                x_y_a_b_values['XMIN_second'],
                x_y_a_b_values['YMIN_second'],
                x_y_a_b_values['XMAX_second'],
                x_y_a_b_values['YMAX_second'],
                errors["err_x_second"],
                errors["err_y_second"]
            ):
                f.write(f"{ra_hms} {dec_dms} {x} {y} {erroreX} {erroreY} {a} {b} {xmin} {ymin} {xmax} {ymax}\n")

def main():
    DIR = 'TMP/'
    fn = 'k1-impTEST.fts.sx'

    try:
        fits_filename, base_filename = choose_fits_file()
    except (FileNotFoundError, ValueError) as e:
        print(e)
        return

    X,Y,ERRX,ERRY,A,B,XMIN,YMIN,XMAX,YMAX,TH,FLAG,FLUX= load_data(f'{DIR}{fn}')

    X,Y,ERRX,ERRY,A,B,XMIN,YMIN,XMAX,YMAX,TH,FLAG,FLUX = preprocess_data(X,Y,ERRX,ERRY,A,B,XMIN,YMIN,XMAX,YMAX,TH,FLAG,FLUX, x_min=100, x_max=3100, y_min=50, y_max=2105)

    ELONG = compute_elongation(A, B)

    likely_satelite = ab_ratio(ELONG, threshold=5)

    erroreX, erroreY = compute_errores(ERRX, ERRY)

    hight_flux = FLUX >= 1000

    outlier_indices = is_outlier(ELONG)
    satellites = np.zeros(len(ELONG), dtype=bool)
    satellites[outlier_indices] = True
    # print(f'ELONG: {likely_satelite}, FLUX: {hight_flux}')

    features = np.column_stack((TH, likely_satelite, hight_flux.astype(int)))
    # print(f'features: {features}')
    features_scaled = scale_features(features)

    # Первая кластеризация для обнаружения потенциального спутника
    labels_first_cluster = cluster_data(features_scaled, eps=5, min_samples=5)
    satellite_mask = labels_first_cluster == -1
    satellite_x_coords = X[satellite_mask]
    satellite_y_coords = Y[satellite_mask]
    print(f'First cluster - X: {satellite_x_coords} Y: {satellite_y_coords}')

    with fits.open(fits_filename) as hdul:
        wcs = WCS(hdul[0].header)
    sky_coords_first = wcs.pixel_to_world(satellite_x_coords, satellite_y_coords)

    print("RA and DEC: ")
    coords_first = []
    for coord in sky_coords_first:
        ra_deg = coord.ra.deg
        dec_deg = coord.dec.deg

        ra_hms, dec_dms = convert(ra_deg, dec_deg)
        print(f"RA: {ra_hms}, DEC: {dec_dms}")
        coords_first.append((ra_hms, dec_dms))


    # Вторая кластеризация для обнаружения аномалий (треков)
    # Исключаем объекты, найденные в первой кластеризации
    non_satellite_mask = ~satellite_mask
    X_non_satellite = X[non_satellite_mask]
    Y_non_satellite = Y[non_satellite_mask]
    TH_non_satellite = TH[non_satellite_mask]
    ELONG_non_satellite = ELONG[non_satellite_mask]

    # Используем ELONG и TH в качестве признаков для обнаружения аномалий
    features_anomalies = np.column_stack((ELONG_non_satellite, TH_non_satellite))
    features_anomalies_scaled = scale_features(features_anomalies)

    labels_anomalies = cluster_data(features_anomalies_scaled, eps=0.5, min_samples=3)
    anomaly_mask = labels_anomalies == -1
    anomaly_x_coords = X_non_satellite[anomaly_mask]
    anomaly_y_coords = Y_non_satellite[anomaly_mask]
    print(f'Second cluster - X: {anomaly_x_coords} Y: {anomaly_y_coords}')

    # Вывод координат для второй кластеризации
    sky_coords_second = wcs.pixel_to_world(anomaly_x_coords, anomaly_y_coords)

    print("Second cluster: ")
    coords_second = []
    for coord in sky_coords_second:
        ra_deg = coord.ra.deg
        dec_deg = coord.dec.deg

        ra_hms, dec_dms = convert(ra_deg, dec_deg)
        print(f"RA: {ra_hms}, DEC: {dec_dms}")
        coords_second.append((ra_hms, dec_dms))

    # Сохраняем результаты обеих кластеризаций
    x_y_a_b_values = {
        'X_first': satellite_x_coords,
        'Y_first': satellite_y_coords,
        'A_first': A[satellite_mask],
        'B_first': B[satellite_mask],
        'XMIN_first': XMIN[satellite_mask],
        'YMIN_first': YMIN[satellite_mask],
        'XMAX_first': XMAX[satellite_mask],
        'YMAX_first': YMAX[satellite_mask],

        'X_second': anomaly_x_coords,
        'Y_second': anomaly_y_coords,
        'A_second': A[non_satellite_mask][anomaly_mask],
        'B_second': B[non_satellite_mask][anomaly_mask],
        'XMIN_second': XMIN[non_satellite_mask][anomaly_mask],
        'YMIN_second': YMIN[non_satellite_mask][anomaly_mask],
        'XMAX_second': XMAX[non_satellite_mask][anomaly_mask],
        'YMAX_second': YMAX[non_satellite_mask][anomaly_mask],
    }

    errors = {
        "err_x_first": erroreX[satellite_mask],
        "err_y_first": erroreY[satellite_mask],
        "err_x_second": erroreX[non_satellite_mask][anomaly_mask],
        "err_y_second": erroreY[non_satellite_mask][anomaly_mask]

    }

    # Сохранение результатов
    save_results(coords_first, coords_second, base_filename, fits_filename, x_y_a_b_values, errors)

    # Визуализация 1: TH vs ELONG
    fig1, ax1 = plt.subplots(figsize=(12, 9))

    # Все объекты
    ax1.scatter(TH, ELONG, c='blue', s=10, label='All data')

    # «Спутники» по первой кластеризации
    if np.any(satellite_mask):
        ax1.scatter(TH[satellite_mask], ELONG[satellite_mask],
                    facecolors='none', edgecolors='red', s=200, label='First cluster (potential satellite)')

    # «Аномалии» по второй кластеризации
    if np.any(anomaly_mask):
        ax1.scatter(TH_non_satellite[anomaly_mask], ELONG_non_satellite[anomaly_mask],
                    facecolors='none', edgecolors='green', s=200,
                    label='Second cluster (anomalies)')

    ax1.set_xlabel('Angle (TH)')
    ax1.set_ylabel('A/B')
    ax1.set_title('DBSCAN in (TH, A/B) feature space')
    ax1.legend()
    ax1.grid(True)
    plt.tight_layout()
    plt.show()

    # Визуализация 2: ELONG vs TH
    fig2, ax2 = plt.subplots(figsize=(12, 9))

    # Все объекты
    ax2.scatter(ELONG, TH, c='blue', s=10, label='All data')

    # «Спутники» по первой кластеризации
    if np.any(satellite_mask):
        ax2.scatter(ELONG[satellite_mask], TH[satellite_mask],
                    facecolors='none', edgecolors='red', s=200, label='First cluster (potential satellite)')

    # «Аномалии» по второй кластеризации
    if np.any(anomaly_mask):
        ax2.scatter(ELONG_non_satellite[anomaly_mask], TH_non_satellite[anomaly_mask],
                    facecolors='none', edgecolors='green', s=200,
                    label='Second cluster (anomalies)')

    ax2.set_xlabel('A/B')
    ax2.set_ylabel('Angle (TH)')
    ax2.set_title('DBSCAN in (A/B, TH) feature space')
    ax2.legend()
    ax2.grid(True)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()