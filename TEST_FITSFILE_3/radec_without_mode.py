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
from astropy.table import QTable

def is_outlier(points, thresh=3.5):
    if len(points.shape) == 1:
        points = points[:, None]
    median = np.median(points, axis=0)
    diff = np.sqrt(np.sum((points - median) ** 2, axis=-1))
    med_abs_deviation = np.median(diff)
    modified_z_score = 0.6745 * diff / med_abs_deviation

    return np.where(modified_z_score > thresh)[0]


def load_data(file_path):
    """Загружает данные из файла в QTable."""
    try:
        # Попробуем прочитать как ASCII-таблицу с разделителями
        data = QTable.read(file_path, format='ascii.fast_no_header', delimiter=' ', names=('X', 'Y', 'ERRX', 'ERRY', 'A', 'B', 'XMIN', 'YMIN', 'XMAX', 'YMAX', 'TH', 'FLAG', 'FLUX'))
        return data
    except Exception as e:
        print(f"Ошибка при чтении файла как ASCII: {e}")
        raise


def preprocess_data(table, x_min, x_max, y_min, y_max):
    """Предварительная обработка данных с использованием QTable."""
    y_mask = (table['Y'] >= y_min) & (table['Y'] <= y_max)
    x_mask = (table['X'] >= x_min) & (table['X'] <= x_max)
    mask = x_mask & y_mask
    return table[mask]


def compute_elongation(A, B):
    """Вычисляет элонгацию (соотношение A/B)."""
    return A / B

def scale_features(features):
    """Масштабирует признаки."""
    scaler = StandardScaler()
    return scaler.fit_transform(features)


def ab_ratio(ELONG, threshold):
    """Определяет объекты с низким соотношением A/B."""
    return ELONG < threshold


def cluster_data(features_scaled, eps, min_samples):
    """Выполняет кластеризацию DBSCAN."""
    dbscan = DBSCAN(eps=eps, min_samples=min_samples)
    return dbscan.fit_predict(features_scaled)


def compute_errores(ERRX, ERRY):
    """Вычисляет ошибки на основе ERRX и ERRY."""
    erroreX = np.sqrt(1/ERRX)
    erroreY = np.sqrt(1/ERRY)
    return erroreX, erroreY


def choose_fits_file():
    """Выбирает FITS-файл на основе содержимого лог-файла."""
    log_file = 'TMP/processing_log.txt'
    if not os.path.exists(log_file):
        raise FileNotFoundError("Log file does not exist.")

    with open(log_file, 'r') as f:
        last_file = f.readline().strip()

    if not last_file:
        raise ValueError("No file name found in the log file.")

    last_file_base = os.path.basename(last_file)
    last_file_suffix = last_file_base.split('.')[0]
    print(f"Extracted suffix: {last_file_suffix}")

    fits_files = glob.glob('TMP/*.fits') + glob.glob('TMP/*.fit')

    print("Found FITS files:")
    for f in fits_files:
        print(f)

    matching_file = None
    for fits_file in fits_files:
        file_name = os.path.basename(fits_file)
        print(f"Checking file: {file_name}")
        if last_file_suffix in file_name:
            matching_file = fits_file
            break

    if not matching_file:
        raise FileNotFoundError(f"No matching FITS file found for suffix {last_file_suffix}")

    print(f"Selected FITS file: {matching_file}")

    return matching_file, last_file_base


def convert(ra_deg, dec_deg):
    """Преобразует градусы в формат ЧЧ:ММ:СС и ГГ:ММ:СС."""
    ra_angle = Angle(ra_deg, unit=u.deg)
    ra_hms = ra_angle.to_string(unit=u.hour, sep=':', precision=2)

    dec_angle = Angle(dec_deg, unit=u.deg)
    dec_dms = dec_angle.to_string(unit=u.deg, sep=':', precision=2, alwayssign=True)

    return ra_hms, dec_dms


def save_results(coords_first, coords_second, base_filename, fits_filename, x_y_a_b_values, errors):
    """Сохраняет результаты кластеризаций."""
    output_dir = 'PROCESS_FILE'
    os.makedirs(output_dir, exist_ok=True)

    txt_filename = os.path.join(output_dir, f'{base_filename}.txt')
    with open(txt_filename, 'w') as f:
        f.write(f"File: {base_filename}\n")

        with fits.open(fits_filename) as hdul:
            header = hdul[0].header
            date_obs = header.get('DATE-OBS', '00000')
            exptime = header.get('EXPTIME', 0)
            if date_obs != '00000' and exptime > 0:
                date_str, time_str = date_obs.split('T')
                time_obs = Time(f'{date_str} {time_str}', format='iso')
                avg_exposure_time = time_obs + timedelta(seconds=(exptime / 2.0))
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

    try:
        data_table = load_data(f'{DIR}{fn}')
    except Exception as e:
        print(f"Ошибка при загрузке данных: {e}")
        return

    processed_table = preprocess_data(data_table, x_min=100, x_max=3100, y_min=50, y_max=2105)

    X = processed_table['X']
    Y = processed_table['Y']
    ERRX = processed_table['ERRX']
    ERRY = processed_table['ERRY']
    A = processed_table['A']
    B = processed_table['B']
    XMIN = processed_table['XMIN']
    YMIN = processed_table['YMIN']
    XMAX = processed_table['XMAX']
    YMAX = processed_table['YMAX']
    TH = processed_table['TH']
    FLAG = processed_table['FLAG']
    FLUX = processed_table['FLUX']

    ELONG = compute_elongation(A, B)

    likely_satelite = ab_ratio(ELONG, threshold=5)

    erroreX, erroreY = compute_errores(ERRX, ERRY)

    hight_flux = FLUX >= 1000

    outlier_indices = is_outlier(ELONG)
    satellites = np.zeros(len(ELONG), dtype=bool)
    satellites[outlier_indices] = True

    features = np.column_stack((TH, likely_satelite, hight_flux.astype(int)))
    features_scaled = scale_features(features)

    # Первая кластеризация для обнаружения потенциального спутника
    labels_first_cluster = cluster_data(features_scaled, eps=5, min_samples=5)
    satellite_mask = labels_first_cluster == -1
    satellite_x_coords = X[satellite_mask]
    satellite_y_coords = Y[satellite_mask]
    # print(f'First cluster:\n{satellite_x_coords}{satellite_y_coords}\n')
    print(f"First cluster:")
    for x, y in zip(satellite_x_coords, satellite_y_coords):
        print(f" X: {x}, Y: {y}")


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
    # print(f'Second cluster:\n{anomaly_x_coords} {anomaly_y_coords}\n')
    print(f"Second cluster:")
    for x, y in zip(anomaly_x_coords, anomaly_y_coords):
        print(f" X: {x}, Y: {y}")

    # Вывод координат для второй кластеризации
    sky_coords_second = wcs.pixel_to_world(anomaly_x_coords, anomaly_y_coords)

    print("RA and DEC: ")
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
                    facecolors='none', edgecolors='red', s=200, label='First cluster')

    # «Аномалии» по второй кластеризации
    if np.any(anomaly_mask):
        ax1.scatter(TH_non_satellite[anomaly_mask], ELONG_non_satellite[anomaly_mask],
                    facecolors='none', edgecolors='green', s=200,
                    label='Second cluster')

    ax1.set_xlabel('Angle (TH)')
    ax1.set_ylabel('A/B')
    ax1.set_title('DBSCAN in (TH, A/B) feature space')
    ax1.legend()
    ax1.grid(True)
    plt.tight_layout()
    plt.show()




if __name__ == "__main__":
    main()