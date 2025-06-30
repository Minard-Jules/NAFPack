# NAFPack - Numerical Analysis in Fortran Package

NAFPack is a Fortran-based numerical analysis package, offering a comprehensive set of algorithms for diverse numerical computations. These computations include Fast Fourier Transform, linear system solving, and eigenvalue/eigenvector calculations.

## Table of Contents
- [Features](#features)
- [Compilation](#Compilation)
- [Getting Started](#getting-started)
- [License](#License)
- [Credits](#credits)

## Features

NAFPack incorporates the following numerical analysis algorithms:

- **Fast Fourier Transform (FFT):** Implementation of the Cooley-Tukey FFT algorithm and integration with the [**FFTW**](https://www.fftw.org/) library for signal processing and frequency domain analysis.

- **Linear System Solvers:** Comprehensive implementations of various methods for solving linear systems, including Gaussian elimination, LU decomposition, and iterative methods like Gauss Seidel.

- **Eigenvalue and Eigenvector Computations:** Algorithms designed to find eigenvalues and eigenvectors of matrices, such as the power iteration method and the QR algorithm.

## Compilation

Ensure you have [CMake](https://cmake.org/download/), [**FFTW**](https://www.fftw.org/), [fpm](https://fpm.fortran-lang.org/en/latest/installation/)  and [**Fortran Compiler**](https://fortran-lang.org/compilers/) (e.g. gfortran) installed.

1. Linux (Debian/Ubuntu)

    ```sh
    $ sudo apt-get update
    $ sudo apt-get install gfortran cmake pkg-config libfftw3-dev
    ```

    ```sh
    # Install fpm
    $ curl -Lo fpm https://github.com/fortran-lang/fpm/releases/download/fpm-linux-x86_64
    $ chmod +x fpm
    $ sudo mv fpm /usr/local/bin
    ```

2. Windows (MSYS2)

    Download and execute the MSYS2 installer, then update MSYS2 as explained on the site [MSYS2](https://www.msys2.org/)

    **Note:** All commands must be run from the MSYS2 MinGW terminal

    - Build tools (gfortran, python, make, pkgconf, ...): 
    ```sh
    $ pacman -S mingw-w64-ucrt-x86_64-toolchain base-devel
    ```
    - CMake:
    ```sh
    $ pacman -S mingw-w64-ucrt-x86_64-cmake
    ```
    - FFTW:
    ```sh
    $ pacman -S mingw-w64-ucrt-x86_64-fftw
    ```
    - fpm:
     ```sh
    $ pacman -S mingw-w64-ucrt-x86_64-fpm
    ```

    

#### Using CMake


1. Create a `build` directory:
    ```sh
    $ mkdir build && cd build
    ```

2. Run CMake to generate build files:
    - Linux
    ```sh
    $ cmake ..
    ```
    - Windows
    ```sh
    $ cmake -G "MSYS Makefiles" ..
    ```

3.  Compile the project with make :
    ```sh
    $ make
    ```

4. Run the tests:
    - Linux
    ```sh
    $ ./bin/main_test
    ```
    - Windows
    ```sh
    $ .\bin\main_test.exe
    ```

#### Using Fortran Package Manager (fpm)


1. Build the project:
    ```sh
    $ fpm build
    ```

2. Run the tests:
    ```sh
    $ fpm test
    ```

## Getting Started

-  To integrate NAFPack into your project, you can use the Fortran Package Manager [**fpm**](https://fpm.fortran-lang.org/).

    In your `fpm.toml` file, add the following dependency:
    ```toml
    [dependencies]
    NAFPack = { git = "https://github.com/Minard-Jules/NAFPack.git" }
    ```

# License

This project is licensed under the MIT License - see the [LICENSE](LICENSE.md) file for details.


# Credits

* [**Minard Jules**](https://github.com/Minard-Jules): Creator of the project.
