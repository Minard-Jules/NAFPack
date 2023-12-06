# NAFPack - Fortran Numerical Analysis Package

NAFPack is a Fortran-based numerical analysis package, offering a comprehensive set of algorithms for diverse numerical computations. These computations include Fast Fourier Transform, linear system solving, and eigenvalue/eigenvector calculations.

## Table of Contents
- [Features](#features)
- [Getting Started](#getting-started)
- [Contributing](#contributing)
- [Credits](#credits)

## Features

NAFPack incorporates the following numerical analysis algorithms:

- **Fast Fourier Transform (FFT):** Implementation of the Cooley-Tukey FFT algorithm and integration with the [**FFTW**](https://www.fftw.org/) library for signal processing and frequency domain analysis.

- **Linear System Solvers:** Comprehensive implementations of various methods for solving linear systems, including Gaussian elimination, LU decomposition, and iterative methods like Gauss Seidel.

- **Eigenvalue and Eigenvector Computations:** Algorithms designed to find eigenvalues and eigenvectors of matrices, such as the power iteration method and the QR algorithm.

## Getting Started

To integrate NAFPack into your project, utilize the Fortran Package Manager [**fpm**](https://fpm.fortran-lang.org/).

1. In your `fpm.toml` file, add the following dependency:
    ```toml
    [dependencies]
    NAFPack = { git = "https://github.com/Minard-Jules/NAFPack.git" }
    ```

2. Enjoy the power of NAFPack for your numerical analysis tasks! If you encounter any questions, issues, or have suggestions, please feel free to open an issue on the [GitHub repository](https://github.com/Minard-Jules/NAFPack).

# Credits

* [**Minard Jules**](https://github.com/Minard-Jules): Creator of the project.
