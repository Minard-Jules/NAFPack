# NAFPack - Fortran Numerical Analysis Package

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

Make sure you have a Fortran compiler (e.g., `gfortran`) installed on your system.

You can compile NAFPack using either `make` or the Fortran Package Manager (`fpm`):

### Using Make

1. Add your compiler to the [**Makefile**](Makefile)
    ```Makefile
    FC = your compiler (e.g., `gfortran`)
    ...
    ```

2. Clone the repository:
    ```sh
    git clone https://github.com/Minard-Jules/NAFPack.git
    cd NAFPack
    ```

3. Build the project:
    ```sh
    make 
    ```

4. Run the tests:
    ```sh
    make run_test
    ```

### Using Fortran Package Manager (fpm)

1. Ensure you have [fpm](https://fpm.fortran-lang.org/en/latest/installation/) installed.

2. Build the project:
    ```sh
    fpm build
    ```

3. Run the tests:
    ```sh
    fpm test
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
