#!/bin/sh
# FFTW Verification Script for Linux/Mac

echo "int main() { return 0; }" > test_fftw_link.c

gfortran -o test_fftw_link test_fftw_link.c -lfftw3 -lfftw3_threads >/dev/null 2>&1

if [ -f test_fftw_link ]; then
    echo "[OK] FFTW is installed."
    rm test_fftw_link.c test_fftw_link
    exit 0
else
    echo "[ERREUR] FFTW is not installed or not found !"
    rm test_fftw_link.c
    exit 1
fi