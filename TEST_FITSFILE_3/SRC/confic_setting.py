import subprocess
import re 
import sys
import numpy as np
from astropy.io import fits 


def extract_header(fits_file):
    with fits.open(fits_file) as hdul:
        header = hdul[0].header
        xpixsz = header["XPIXSZ"]
        gain = header["GAIN"]
        exptime = header["EXPTIME"]    
        
    return xpixsz, gain, exptime


def calcul_scale(fits_file):
    result = subprocess.run(["solve-field", "--scale-units", "arcsecperpix", "--no-plots", "--overwrite", "--silent", fits_file], capture_output=True, text=True)
    
    scale_match = re.search(r"Scale: (\d+\,\d+)", result.stdout)
    if not scale_match:
        raise RuntimeError("Failed to calculate")
    
    scale = float(scale_match.group(1))
    
    return scale 


def calculate_param(xpixsz, gain, exptime, scale):
    D = xpixsz * scale
    DETECT_MINAPEA = (D/ scale) ** 2
    
    signal = gain * exptime
    noise = np.sqrt(signal)
    DETECT_THRESH = signal / noise
    
    return DETECT_MINAPEA, DETECT_THRESH


def main():
    fits_file = "Meteor1_19R_00001.fits"
    # fits_file = sys.argv[1]
    config_file = "CONFIGS/defaultTEST.sex"
    
    xpixsz, gain, exptime = extract_header(fits_file)
    scale = calcul_scale(fits_file)
    detect_minarea, detect_thresh = calculate_param(xpixsz, gain, exptime, scale)
    
    sextrator_command = ["sex", fits_file,
                         "-c", config_file,
                         "-DETECT_MINAREA", f"{detect_minarea:.2f},",
                         "-DETECT_THRESH", f"{detect_thresh:.2f}"]
    
    subprocess.run(sextrator_command)
    
    print(f"DETECT_MINAREA: {detect_minarea:2f}")
    print(f"DETECT_THRESH: {detect_thresh:2f}")
    
    
    