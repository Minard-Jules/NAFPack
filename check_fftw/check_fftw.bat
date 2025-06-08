@echo off
REM FFTW Verification Script for Windows (cmd.exe)
setlocal

REM Create a minimal C test file
echo int main() { return 0; } > test_fftw_link.c

REM Attempt to link with FFTW
gfortran -o test_fftw_link.exe test_fftw_link.c -lfftw3 -lfftw3_threads >nul 2>&1

REM Checks if the executable has been created
if exist test_fftw_link.exe (
    echo [OK] FFTW is installed.
    del test_fftw_link.c
    del test_fftw_link.exe
    exit /b 0
) else (
    echo [ERREUR] FFTW is not installed or not found !
    del test_fftw_link.c
    exit /b 1
)