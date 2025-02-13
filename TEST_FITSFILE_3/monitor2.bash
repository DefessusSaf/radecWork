#!/bin/bash
MONITORDIR="LOAD_FILE/"
LOCKFILE="/tmp/process.lock"
TMPDIR="TMP"
LOGFILE="${TMPDIR}/accuracy_log.txt"

# Очистка TMP перед началом работы
clear_tmp_dir() {
    echo "cleaning ${TMPDIR}..."
    rm -rf ${TMPDIR}/*
    echo "Done."
}

process_file() {
    NEWFILE="$1"
    echo "${NEWFILE} uploaded"

    echo "${NEWFILE}" >> ${TMPDIR}/processing_log.txt

    sleep 1
    echo '===========  Sources extractions ======================='
    python confic_setting.py "${NEWFILE}"
    # python config_setting_without_focallen.py "${NEWFILE}"
    

    sleep 1
    # echo '===========  Sources extractions ======================='
    # sex "${NEWFILE}" -c CONFIGS/default.sex
    # sex "${NEWFILE}" -c CONFIGS/defaultTEST.sex

    echo '========================================================'
    python SRC/WriteOutRegionFile_18072024.py

    echo '===========  astrometry up TMP directory  =============='
    python SRC/WFT_19072024.py "${NEWFILE}"

    echo '===========  astrometry up TMP directory  =============='
    # python SRC/xy2radec.py "${NEWFILE}"

    echo '===========  astrometry up TMP directory  =============='
    python radec_without_mode.py

    if [ ! -s TMP/XY.fits ]; then
        echo "Error: TMP/XY.fits does not exist or is empty"
        exit 1
    fi

    echo "${NEWFILE}" > ${TMPDIR}/processing_log.txt

    rm -f "${NEWFILE}"
    : > ${TMPDIR}/processing_log.txt

    echo '==================  DONE  =============================='
}

# Очистка TMP перед запуском мониторинга
clear_tmp_dir

inotifywait -m -e create --format '%w%f' "${MONITORDIR}" | while read NEWFILE
do
    {
        exec 200>$LOCKFILE
        flock -x 200
        process_file "${NEWFILE}"
        # Очистка TMP после обработки каждого файла (если нужно)
        clear_tmp_dir
    }
done