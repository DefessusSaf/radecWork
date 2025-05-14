import subprocess
import sys
import os
import shutil
import timeit
import numpy as np
from astropy.io import fits
import astropy.units as u
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

FOCALLEN = 551 * u.mm


def extract_header_and_scale(fits_file):
    """
    It extracts the necessary data from the tits file header and calculates the image scale.    """
    start = timeit.default_timer()
    logging.info(f"Opening FITS file: {fits_file}")
    try:
        with fits.open(fits_file) as hdul:
            if not hdul:
                raise ValueError("Failed to open a file file or it is empty.")
            header = hdul[0].header

            try:
                xpixsz_um = header["XPIXSZ"]
                exptime = header["EXPTIME"]
            except KeyError as e:
                logging.error(f"Отсутствует обязательный ключ в заголовке: {e}")
                raise

            if "EGAIN" in header:
                gain = header["EGAIN"]
                logging.info(f"Used EGAIN: {gain}")
            elif "GAIN" in header:
                gain = header["GAIN"]
                logging.info(f"Used GAIN: {gain}")
            else:
                raise KeyError("Absent GAIN or EGAIN in FITS heading.")

            pixel_size = xpixsz_um * u.micron

            try:
                # 1.We calculate the dimensionless attitude (pixel size / focal distance)
                dimensionless_ratio = (pixel_size / FOCALLEN).to(u.dimensionless_unscaled)

                # 2.Physically, this attitude corresponds to the corner in radiana for pixel
                #   We assign the right units
                scale_rad_per_pix = dimensionless_ratio * u.rad / u.pixel

                # 3. We convert in the corner seconds to the pixel
                scale_arcsec_per_pix_quantity = scale_rad_per_pix.to(u.arcsec / u.pixel)
                scale_arcsec_per_pix = scale_arcsec_per_pix_quantity.value

                logging.info(f"Pixel size (XPIXSZ): {xpixsz_um} um")
                logging.info(f"Focal length: {FOCALLEN}")
                logging.info(f"Calculated scale: {scale_arcsec_per_pix_quantity:.4f}")

            except u.UnitConversionError as e:
                 logging.error(f"Error in converting units when calculating the scale: {e}")
                 raise ValueError("Incompatible units for calculating scale.") from e

    except FileNotFoundError:
        logging.error(f"File not found: {fits_file}")
        raise
    except Exception as e:
        logging.error(f"There was an error when reading the header FITS: {e}")
        raise

    time_taken = timeit.default_timer() - start
    logging.info(f"The title is extracted and the scale is designed for {time_taken:.4f} seconds.")

    return xpixsz_um, gain, exptime, scale_arcsec_per_pix


def calculate_initial_sextractor_params(scale_arcsec_per_pix):
    """
    Calculates the initial parameters for SExtractor.

    Args:
        Scale_arcsec_per_pix (Float): Image scale in arcsec/pixel (for information).

    Returns:
        Tuple: Cortege (Detect_Minarea, Detect_thresh)
        - Detect_minarea (int): the initial minimum area of ​​the object in pixels.
        - Detect_thresh (Float): the initial detection threshold (in the Sigma above the background).
    """
    start = timeit.default_timer()

    DETECT_MINAREA = 9
    logging.info(f"The initial value of Detect_Minarea is set in {DETECT_MINAREA} pixels.")

    DETECT_THRESH = 1.5
    logging.info(f"The initial value of Detect_thresh is installed in {DETECT_THRESH} sigma.")

    time_taken = timeit.default_timer() - start
    logging.info(f"The initial SExtractor parameters are designed for {time_taken} seconds.")
    return DETECT_MINAREA, DETECT_THRESH


def run_sextractor_and_count(fits_file, config_file, temp_dir, detect_minarea, detect_thresh):
    """
    Launches SExtractor with specified parameters and calculates the number of objects discovered
    From *a rigidly set catalog file *.

    Args:
        Fits_file (str): the path to the input fits file.
        Config_file (str): the path to the configuration file SExtractor (.SEX).
        TEMP_DIR (str): the path to the temporary directory for some output files (for example, Checkimage).
        Detect_Minarea (Float): The value of the Detect_Minarea parameter.
        DETECT_THRESH (FLOAT): Detect_thresh parameter value.

    Returns:
        Int: The number of objects detected.
    RAISES:
        Filenotfounderror: If the SExtractor configuration file is not found or the expected catalog is not created.
        Runtimeerror: if SExtractor ends with an error.
    """
    start = timeit.default_timer()

    if not os.path.exists(config_file):
         logging.error(f"SExtractor configuration file was not found: {config_file}")
         raise FileNotFoundError(f"SExtractor configuration file was not found: {config_file}")

    catalog_file_path = "TMP/k1-impTEST.fts.sx"
    catalog_dir = os.path.dirname(catalog_file_path)

    if catalog_dir and not os.path.exists(catalog_dir):
        try:
            os.makedirs(catalog_dir, exist_ok=True)
            logging.info(f"Created a directory for the catalog: {catalog_dir}")
        except OSError as e:
            logging.error(f"Failed to create a directory {catalog_dir}: {e}")

    base_name = os.path.splitext(os.path.basename(fits_file))[0]
    check_image_path = os.path.join(temp_dir, f"{base_name}_check.fits")

    sextractor_command = [
        "sex", fits_file, 
        "-c", config_file,
        "-DETECT_MINAREA", f"{detect_minarea:.2f}",
        "-DETECT_THRESH", f"{detect_thresh:.2f}",
        "-CHECKIMAGE_NAME", check_image_path # Check Image is stored in a temporary folder
    ]

    logging.info(f"Запуск SExtractor с DETECT_MINAREA={detect_minarea:.2f}, DETECT_THRESH={detect_thresh:.2f}")
    logging.info(f"Ожидаемый файл каталога: {catalog_file_path}")
    logging.debug(f"Команда: {' '.join(sextractor_command)}")

    try:
        result = subprocess.run(sextractor_command, capture_output=True, text=True, check=False)

        if result.returncode != 0:
            logging.error(f"SExtractor завершился с ошибкой (код {result.returncode}).")
            logging.error(f"Stderr:\n{result.stderr}")
            # IMPORTANT: Even if SEXTRACTOR has fallen, it could create an empty or incomplete directory
            # Catalog's existence will be lower than
            raise RuntimeError(f"Error execution SExtractor.")
        else:
            logging.info("SExtractor Successfully completed.")
            logging.debug(f"Stdout:\n{result.stdout}")
            if result.stderr:
                 logging.warning(f"Stderr (Possible warnings):\n{result.stderr}")

        logging.info(f"Reading the catalog file: {catalog_file_path}")
        if not os.path.exists(catalog_file_path):
             # This error can now occur if Sextractor has fallen, and if it worked, but did not create a file
             logging.error(f"The expected catalog file was not found: {catalog_file_path}")
             raise FileNotFoundError(f"Catalog file {catalog_file_path} not found after execution SExtractor.")

        with open(catalog_file_path, "r") as catalog:
            lines = catalog.readlines()
        object_count = len(lines) 

    except FileNotFoundError as e:
        raise
    except RuntimeError as e:
        raise
    except Exception as e:
        logging.error(f"An error occurred during the performance or processing of the results of SExtractor: {e}")
        raise

    time_taken = timeit.default_timer() - start
    logging.info(f"SExtractor I processed the file and found {object_count} objects (from {catalog_file_path}) for {time_taken:.4f} seconds.")

    return object_count


def cleanup_temp_dir(temp_dir):
    """Removes the temporary directory and its contents."""
    start = timeit.default_timer()
    if os.path.exists(temp_dir):
        try:
            shutil.rmtree(temp_dir)
            logging.info(f"Temporary directory {temp_dir} successfully removed.")
        except OSError as e:
            logging.error(f"Failed to remove the temporary directory {temp_dir}: {e}")
    else:
        logging.info(f"Temporary directory {temp_dir} not found, cleaning is not required")
    time_taken = timeit.default_timer() - start
    logging.info(f"Cleaning is performed for {time_taken:.4f} seconds.")


def main():
    """
    The main function of the script.Processing fits file, it is launched by SExtractor
    To achieve the target number of objects and cleans temporary files.
    """
    if len(sys.argv) < 2:
        print(f"Usage: python {sys.argv[0]} <Way_fits_Fileу>")
        sys.exit(1)

    fits_file = sys.argv[1]
    temp_dir = "./temp_sextractor_files"
    config_file = "CONFIGS/defaultTEST.sex"

    if not os.path.exists(fits_file):
        logging.error(f"The entrance FITS file was not found: {fits_file}")
        sys.exit(1)
    if not os.path.exists(config_file):
         logging.error(f"SExtractor configuration file was not found: {config_file}")
         sys.exit(1)

    try:
        os.makedirs(temp_dir, exist_ok=True)
        logging.info(f"The temporary directory has been created/exists: {temp_dir}")
    except OSError as e:
        logging.error(f"Failed to create a temporary directory {temp_dir}: {e}")
        sys.exit(1)

    try:
        xpixsz_um, gain, exptime, scale_arcsec_per_pix = extract_header_and_scale(fits_file)
        logging.info(f"Parameters from heading: XPIXSZ={xpixsz_um} um, GAIN={gain}, EXPTIME={exptime} s")

        detect_minarea, detect_thresh = calculate_initial_sextractor_params(scale_arcsec_per_pix)

        target_objects = 80
        max_iter = 50
        min_minarea = 1
        adjustment_factor = 0.9

        current_minarea = float(detect_minarea)
        final_minarea = current_minarea
        found_target = False

        logging.info(f"The beginning of an iterative search {target_objects} objects (max Iterations: {max_iter}).")

        for iteration in range(max_iter):
            logging.info(f"--- Iteration {iteration + 1}/{max_iter} ---")
            try:
                detect_obj_count = run_sextractor_and_count(
                    fits_file, config_file, temp_dir, current_minarea, detect_thresh
                )

                if detect_obj_count >= target_objects:
                    logging.info(f"Found {detect_obj_count} objects (>= {target_objects}) DETECT_MINAREA = {current_minarea:.2f}")
                    final_minarea = current_minarea
                    found_target = True
                    break

                else:
                    logging.info(f"Found {detect_obj_count} objects (< {target_objects}). Reduce DETECT_MINAREA.")
                    current_minarea *= adjustment_factor
                    if current_minarea < min_minarea:
                        logging.warning(f"DETECT_MINAREA ({current_minarea:.2f}) reached a minimum limit ({min_minarea}). Stopping iterations.")
                        final_minarea = min_minarea
                        try:
                            detect_obj_count = run_sextractor_and_count(
                                fits_file, config_file, temp_dir, final_minarea, detect_thresh
                            )
                            logging.info(f"The number of objects with minimal DETECT_MINAREA ({final_minarea:.2f}): {detect_obj_count}")
                        except Exception as final_run_e:
                            logging.error(f"Error at the last launch SExtractor with min_minarea: {final_run_e}")
                        break

            except (FileNotFoundError, RuntimeError) as e:
                 logging.error(f"An error for iteration {iteration + 1}: {e}. Interruption of the cycle.")
                 cleanup_temp_dir(temp_dir)
                 sys.exit(1)
            except Exception as e:
                 logging.error(f"Неожиданная ошибка на итерации {iteration + 1}: {e}. Прерывание цикла.")
                 cleanup_temp_dir(temp_dir)
                 sys.exit(1)

        if found_target:
            print(f"\nSuccess!A sufficient number of objects found.")
        elif iteration == max_iter - 1 and not found_target:
             logging.warning(f"The iteration limit has been reached ({max_iter}). It was not possible to find {target_objects} objects.")
             print(f"\nIt was not possible to find {target_objects} objects for {max_iter} iterations.")
        else:
             print(f"\nThe search is stopped due to the achievement of the minimum DETECT_MINAREA ({final_minarea:.2f}).")

        print(f"Final value DETECT_MINAREA: {final_minarea:.2f}")
        print(f"Final value DETECT_THRESH: {detect_thresh:.2f}")

    except (FileNotFoundError, KeyError, ValueError, RuntimeError) as e:
        logging.error(f"A critical error occurred: {e}")
        print(f"Error execution: {e}")
        cleanup_temp_dir(temp_dir)
        sys.exit(1)
    except Exception as e:
        logging.exception("An unforeseen mistake occurred:")
        print(f"An unforeseen mistake: {e}")
        cleanup_temp_dir(temp_dir)
        sys.exit(1)
    finally:
        cleanup_temp_dir(temp_dir)
        logging.info("The work of the script is completed.")


if __name__ == "__main__":
    main()