import subprocess
import sys
import os
import shutil
import timeit
import numpy as np
from astropy.io import fits
import astropy.units as u
from astropy.table import Table # Оставляем импорт на всякий случай, но не используем
import logging

# Настройка логирования
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Константа фокусного расстояния с единицами измерения
FOCALLEN = 551 * u.mm



# Извлечение данных из заголовка FITS-файла и расчет масштаба (ИСПРАВЛЕННАЯ ВЕРСИЯ)
def extract_header_and_scale(fits_file):
    """
    It extracts the necessary data from the tits file header and calculates the image scale.    """
    start = timeit.default_timer()
    logging.info(f"Открытие FITS файла: {fits_file}")
    try:
        with fits.open(fits_file) as hdul:
            if not hdul:
                raise ValueError("Не удалось открыть FITS файл или он пуст.")
            header = hdul[0].header

            # Извлечение параметров с проверкой наличия
            try:
                xpixsz_um = header["XPIXSZ"]
                exptime = header["EXPTIME"]
            except KeyError as e:
                logging.error(f"Отсутствует обязательный ключ в заголовке: {e}")
                raise

            # Обработка GAIN/EGAIN
            if "EGAIN" in header:
                gain = header["EGAIN"]
                logging.info(f"Используется EGAIN: {gain}")
            elif "GAIN" in header:
                gain = header["GAIN"]
                logging.info(f"Используется GAIN: {gain}")
            else:
                raise KeyError("Отсутствует GAIN или EGAIN в заголовке FITS.")

            # Присваиваем единицы измерения
            pixel_size = xpixsz_um * u.micron

            # ---- ИСПРАВЛЕННЫЙ РАСЧЕТ МАСШТАБА ----
            try:
                # 1. Вычисляем безразмерное отношение (размер пикселя / фокусное расстояние)
                dimensionless_ratio = (pixel_size / FOCALLEN).to(u.dimensionless_unscaled)

                # 2. Физически это отношение соответствует углу в радианах на пиксель
                #    Присваиваем правильные единицы
                scale_rad_per_pix = dimensionless_ratio * u.rad / u.pixel

                # 3. Конвертируем в угловые секунды на пиксель
                scale_arcsec_per_pix_quantity = scale_rad_per_pix.to(u.arcsec / u.pixel)
                scale_arcsec_per_pix = scale_arcsec_per_pix_quantity.value # Получаем числовое значение

                logging.info(f"Размер пикселя (XPIXSZ): {xpixsz_um} um")
                logging.info(f"Фокусное расстояние: {FOCALLEN}")
                # Используем quantity для вывода с единицами для наглядности
                logging.info(f"Рассчитанный масштаб: {scale_arcsec_per_pix_quantity:.4f}")

            except u.UnitConversionError as e:
                 # Эта ошибка теперь маловероятна здесь, но оставим на всякий случай
                 logging.error(f"Ошибка конвертации единиц при расчете масштаба: {e}")
                 raise ValueError("Несовместимые единицы для расчета масштаба.") from e

    except FileNotFoundError:
        logging.error(f"Файл не найден: {fits_file}")
        raise
    except Exception as e:
        logging.error(f"Произошла ошибка при чтении заголовка FITS: {e}")
        raise

    time_taken = timeit.default_timer() - start
    logging.info(f"Заголовок извлечен и масштаб рассчитан за {time_taken:.4f} секунд.")

    return xpixsz_um, gain, exptime, scale_arcsec_per_pix


# Вычисление начальных параметров DETECT_MINAREA и DETECT_THRESH
def calculate_initial_sextractor_params(scale_arcsec_per_pix):
    """
    Calculates the initial parameters for SEXTRACTOR.

    Args:
        Scale_arcsec_per_pix (Float): Image scale in arcsec/pixel (for information).

    Returns:
        Tuple: Cortege (Detect_Minarea, Detect_thresh)
        - Detect_minarea (int): the initial minimum area of ​​the object in pixels.
        - Detect_thresh (Float): the initial detection threshold (in the Sigma above the background).
    """
    start = timeit.default_timer()

    DETECT_MINAREA = 9
    logging.info(f"Начальное значение DETECT_MINAREA установлено в {DETECT_MINAREA} пикселей.")

    DETECT_THRESH = 1.5
    logging.info(f"Начальное значение DETECT_THRESH установлено в {DETECT_THRESH} сигма.")

    time_taken = timeit.default_timer() - start
    logging.info(f"Начальные параметры SExtractor рассчитаны за {time_taken:.4f} секунд.")
    return DETECT_MINAREA, DETECT_THRESH

# Выполнение SExtractor и подсчёт количества объектов (ВОССТАНОВЛЕННАЯ ЛОГИКА КАТАЛОГА)
def run_sextractor_and_count(fits_file, config_file, temp_dir, detect_minarea, detect_thresh):
    """
    Launches SEXTRACTOR with specified parameters and calculates the number of objects discovered
    From *a rigidly set catalog file *.

    Args:
        Fits_file (str): the path to the input fits file.
        Config_file (str): the path to the configuration file SEXTRACTOR (.SEX).
        TEMP_DIR (str): the path to the temporary directory for some output files (for example, Checkimage).
        Detect_Minarea (Float): The value of the Detect_Minarea parameter.
        DETECT_THRESH (FLOAT): Detect_thresh parameter value.

    Returns:
        Int: The number of objects detected.
    RAISES:
        Filenotfounderror: If the SEXTRACTOR configuration file is not found or the expected catalog is not created.
        Runtimeerror: if SEXTRACTOR ends with an error.
    """
    start = timeit.default_timer()

    if not os.path.exists(config_file):
         logging.error(f"Файл конфигурации SExtractor не найден: {config_file}")
         raise FileNotFoundError(f"Файл конфигурации SExtractor не найден: {config_file}")

    # Имя файла каталога
    catalog_file_path = "TMP/k1-impTEST.fts.sx"
    catalog_dir = os.path.dirname(catalog_file_path)

    # Убедимся, что директория TMP существует
    if catalog_dir and not os.path.exists(catalog_dir):
        try:
            os.makedirs(catalog_dir, exist_ok=True)
            logging.info(f"Создана директория для каталога: {catalog_dir}")
        except OSError as e:
            logging.error(f"Не удалось создать директорию {catalog_dir}: {e}")

    # Генерируем имя для check image во временной директории temp_dir
    base_name = os.path.splitext(os.path.basename(fits_file))[0]
    check_image_path = os.path.join(temp_dir, f"{base_name}_check.fits")

    sextractor_command = [
        "sex", fits_file,  # Или "sextractor"
        "-c", config_file,
        "-DETECT_MINAREA", f"{detect_minarea:.2f}",
        "-DETECT_THRESH", f"{detect_thresh:.2f}",
        "-CHECKIMAGE_NAME", check_image_path # Check image сохраняется во временную папку
    ]

    logging.info(f"Запуск SExtractor с DETECT_MINAREA={detect_minarea:.2f}, DETECT_THRESH={detect_thresh:.2f}")
    logging.info(f"Ожидаемый файл каталога: {catalog_file_path}")
    logging.debug(f"Команда: {' '.join(sextractor_command)}")

    try:
        result = subprocess.run(sextractor_command, capture_output=True, text=True, check=False)

        if result.returncode != 0:
            logging.error(f"SExtractor завершился с ошибкой (код {result.returncode}).")
            logging.error(f"Stderr:\n{result.stderr}")
            # Важно: даже если SExtractor упал, он мог создать пустой или неполный каталог
            # Проверка на существование каталога будет ниже
            raise RuntimeError(f"Ошибка выполнения SExtractor. Подробности в логе.")
        else:
            logging.info("SExtractor успешно выполнен.")
            logging.debug(f"Stdout:\n{result.stdout}")
            if result.stderr:
                 logging.warning(f"Stderr (возможны предупреждения):\n{result.stderr}")

        # Чтение жестко заданного каталога объектов
        logging.info(f"Чтение файла каталога: {catalog_file_path}")
        if not os.path.exists(catalog_file_path):
             # Эта ошибка теперь может возникнуть И если SExtractor упал, И если он отработал, но не создал файл
             logging.error(f"Ожидаемый файл каталога НЕ НАЙДЕН: {catalog_file_path}")
             raise FileNotFoundError(f"Файл каталога {catalog_file_path} не найден после выполнения SExtractor.")

        # Подсчет строк, как в оригинале
        with open(catalog_file_path, "r") as catalog:
            lines = catalog.readlines()
        object_count = len(lines) # Оригинальный код считал все строки, включая комментарии

    except FileNotFoundError as e:
        raise
    except RuntimeError as e:
        raise
    except Exception as e:
        logging.error(f"Произошла ошибка во время выполнения или обработки результатов SExtractor: {e}")
        raise

    time_taken = timeit.default_timer() - start
    logging.info(f"SExtractor обработал файл и найдено {object_count} объектов (из {catalog_file_path}) за {time_taken:.4f} секунд.")
    return object_count

# Очистка временной директории после выполнения
def cleanup_temp_dir(temp_dir):
    """Удаляет временную директорию и ее содержимое."""
    start = timeit.default_timer()
    if os.path.exists(temp_dir):
        try:
            shutil.rmtree(temp_dir)
            logging.info(f"Временная директория {temp_dir} успешно удалена.")
        except OSError as e:
            logging.error(f"Не удалось удалить временную директорию {temp_dir}: {e}")
    else:
        logging.info(f"Временная директория {temp_dir} не найдена, очистка не требуется.")
    time_taken = timeit.default_timer() - start
    logging.info(f"Очистка выполнена за {time_taken:.4f} секунд.")

# Основная функция
def main():
    """
    The main function of the script.Processing fits file, it is launched by SEXTRACTOR
    To achieve the target number of objects and cleans temporary files.
    """
    if len(sys.argv) < 2:
        print(f"Использование: python {sys.argv[0]} <путь_к_fits_файлу>")
        sys.exit(1)

    fits_file = sys.argv[1]
    temp_dir = "./temp_sextractor_files"
    config_file = "CONFIGS/defaultTEST.sex"

    if not os.path.exists(fits_file):
        logging.error(f"Входной FITS файл не найден: {fits_file}")
        sys.exit(1)
    if not os.path.exists(config_file):
         logging.error(f"Файл конфигурации SExtractor не найден: {config_file}")
         sys.exit(1)

    try:
        os.makedirs(temp_dir, exist_ok=True)
        logging.info(f"Временная директория создана/существует: {temp_dir}")
    except OSError as e:
        logging.error(f"Не удалось создать временную директорию {temp_dir}: {e}")
        sys.exit(1)

    try:
        xpixsz_um, gain, exptime, scale_arcsec_per_pix = extract_header_and_scale(fits_file) # Используем ИСПРАВЛЕННУЮ функцию
        logging.info(f"Параметры из заголовка: XPIXSZ={xpixsz_um} um, GAIN={gain}, EXPTIME={exptime} s")

        detect_minarea, detect_thresh = calculate_initial_sextractor_params(scale_arcsec_per_pix)

        target_objects = 80
        max_iter = 50
        min_minarea = 1
        adjustment_factor = 0.9

        current_minarea = float(detect_minarea)
        final_minarea = current_minarea
        found_target = False

        logging.info(f"Начало итеративного поиска {target_objects} объектов (макс. итераций: {max_iter}).")

        for iteration in range(max_iter):
            logging.info(f"--- Итерация {iteration + 1}/{max_iter} ---")
            try:
                detect_obj_count = run_sextractor_and_count(
                    fits_file, config_file, temp_dir, current_minarea, detect_thresh
                )

                if detect_obj_count >= target_objects:
                    logging.info(f"Найдено {detect_obj_count} объектов (>= {target_objects}) с DETECT_MINAREA = {current_minarea:.2f}")
                    final_minarea = current_minarea
                    found_target = True
                    break

                else:
                    logging.info(f"Найдено {detect_obj_count} объектов (< {target_objects}). Уменьшаем DETECT_MINAREA.")
                    current_minarea *= adjustment_factor
                    if current_minarea < min_minarea:
                        logging.warning(f"DETECT_MINAREA ({current_minarea:.2f}) достигло минимального предела ({min_minarea}). Остановка итераций.")
                        final_minarea = min_minarea
                        try:
                            detect_obj_count = run_sextractor_and_count(
                                fits_file, config_file, temp_dir, final_minarea, detect_thresh
                            )
                            logging.info(f"Количество объектов при минимальном DETECT_MINAREA ({final_minarea:.2f}): {detect_obj_count}")
                        except Exception as final_run_e:
                            logging.error(f"Ошибка при последнем запуске SExtractor с min_minarea: {final_run_e}")
                        break

            except (FileNotFoundError, RuntimeError) as e:
                 logging.error(f"Ошибка на итерации {iteration + 1}: {e}. Прерывание цикла.")
                 cleanup_temp_dir(temp_dir)
                 sys.exit(1)
            except Exception as e:
                 logging.error(f"Неожиданная ошибка на итерации {iteration + 1}: {e}. Прерывание цикла.")
                 cleanup_temp_dir(temp_dir)
                 sys.exit(1)

        if found_target:
            print(f"\nУспех! Найдено достаточное количество объектов.")
        elif iteration == max_iter - 1 and not found_target:
             logging.warning(f"Достигнут лимит итераций ({max_iter}). Не удалось найти {target_objects} объектов.")
             print(f"\nНе удалось найти {target_objects} объектов за {max_iter} итераций.")
        else:
             print(f"\nПоиск остановлен из-за достижения минимального DETECT_MINAREA ({final_minarea:.2f}).")

        print(f"Итоговое значение DETECT_MINAREA: {final_minarea:.2f}")
        print(f"Итоговое значение DETECT_THRESH: {detect_thresh:.2f}")

    except (FileNotFoundError, KeyError, ValueError, RuntimeError) as e:
        logging.error(f"Произошла критическая ошибка: {e}")
        print(f"Ошибка выполнения: {e}")
        cleanup_temp_dir(temp_dir)
        sys.exit(1)
    except Exception as e:
        logging.exception("Произошла непредвиденная ошибка:")
        print(f"Непредвиденная ошибка: {e}")
        cleanup_temp_dir(temp_dir)
        sys.exit(1)
    finally:
        cleanup_temp_dir(temp_dir)
        logging.info("Работа скрипта завершена.")

if __name__ == "__main__":
    main()