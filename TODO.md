# ✅ NAFPack – TODO & Roadmap

## 📌 Objective
Develop a modular and high-performance Fortran package for numerical analysis, including classical and advanced methods, with parallelism, portability, and interoperability.

---

## ✅ Already Implemented

### 🔢 Linear Systems
- [x] Gaussian Elimination (no pivoting, partial pivoting, full pivoting)
- [x] LU Decomposition
- [x] LDU Decomposition
- [x] Cholesky Decomposition
- [x] QR Decomposition
  - [x] Householder
  - [x] Givens
  - [x] Gram-Schmidt
- [x] Jacobi Iterative Method
- [x] Gauss-Seidel Method
- [x] SOR (Successive Over-Relaxation) Method
- [x] TDMA (Tridiagonal Matrix Algorithm)

### 🧠 Eigenvalue Algorithms
- [x] Power Iteration
- [x] QR Algorithm (with Householder, Givens, and Gram-Schmidt)
- [x] Shifted QR Algorithm

### 🔁 Fourier Transforms
- [x] DFT (1D)
- [x] FFT (1D and 2D)
- [x] IFFT (1D and 2D)

### 📚 External Libraries Integrated
- [x] FFTW
- [x] LAPACK
- [x] BLAS

---

## 🛠️ Features to Add

### 🔧 Performance Improvements
- [ ] Use `DO CONCURRENT` for local parallelism
- [ ] Integrate OpenMP (shared-memory parallelism)
- [ ] Integrate MPI (distributed-memory computation)
- [ ] Modular error handling system (`error stop` or `naf_error_handling` module)

---

### 🧮 Linear Systems

#### 🔹 Direct Methods
- [ ] Crout LU decomposition
- [ ] Doolittle LU decomposition
- [ ] Bartels–Golub stable LU

#### 🔹 Iterative Methods
- [ ] Conjugate Gradient (CG)
- [ ] GMRES
- [ ] BiCGSTAB
- [ ] Preconditioned CG (PCG)

#### 🔹 Utilities
- [ ] Automatic detection of matrix properties: symmetry, positive-definiteness, sparsity

---

### 🧠 Eigenvalues and Eigenvectors

- [ ] Rayleigh quotient iteration
- [ ] Lanczos algorithm
- [ ] Arnoldi iteration
- [ ] Davidson algorithm
- [ ] Jacobi method for symmetric matrices

---

### 🔠 Additional Matrix Decompositions
- [ ] SVD (Singular Value Decomposition)
- [ ] Polar decomposition
- [ ] Schur decomposition

---

### 📐 Interpolation
- [ ] Lagrange interpolation
- [ ] Hermite interpolation
- [ ] Cubic splines
- [ ] Barycentric interpolation

---

### ∫ Numerical Integration

#### 🔹 Basic Methods
- [ ] Trapezoidal rule
- [ ] Simpson's rule
- [ ] Boole's rule

#### 🔹 Advanced Methods
- [ ] Adaptive quadrature
- [ ] Gauss-Legendre quadrature
- [ ] Monte Carlo integration (1D, multiD)

---

### 🧩 Nonlinear Equations
- [ ] Newton-Raphson (1D and multiD)
- [ ] Secant method
- [ ] Brent's method
- [ ] Relaxation method

---

### ⏱️ Differential Equations
- [ ] Runge-Kutta methods (RK2, RK4, RK45)
- [ ] Crank-Nicolson method
- [ ] Symplectic integrators
- [ ] Simple PDE solvers (e.g., 1D diffusion)

---

### 🧮 Advanced Numerical Algebra
- [ ] Square root (iterative methods)
- [ ] Log/exp (Taylor series, CORDIC)
- [ ] Matrix computations:
  - [ ] Condition number
  - [ ] Determinant
  - [ ] Rank
  - [ ] Norms (1, ∞, Frobenius)

---

### 📦 Interfaces and Extensions
- [ ] Python interface via `f2py`
- [ ] Matrix Market file support (.mtx)
- [ ] Matrix generator (Hilbert, diagonally dominant, etc.)
- [ ] CSV export for benchmarking results

---

### 🧪 Tests & Documentation
- [ ] Integrate `pFUnit` for unit testing
- [ ] Documentation generation (FORD or Doxygen)
- [ ] Add example use cases for each module

---

## 📅 Development Roadmap

### 🔹 Phase 1 – Core Optimizations (July)
- [ ] Apply `DO CONCURRENT` parallelism
- [ ] Refactor existing modules
- [ ] Start unit testing integration

### 🔹 Phase 2 – Advanced Methods (August - September)
- [ ] Add CG, GMRES, BiCGSTAB
- [ ] Add eigenvalue algorithms: Arnoldi, Lanczos
- [ ] Add SVD and matrix decompositions
- [ ] Add interpolation and numerical integration methods

### 🔹 Phase 3 – Interfaces & Extensions (October)
- [ ] Integrate Python interface via `f2py`
- [ ] Add Matrix Market format support
- [ ] Implement test matrix generators

### 🔹 Phase 4 – Symbolic Tools & ODEs (November)
- [ ] Add ODE solvers
- [ ] Implement tools for condition number, norm, stability analysis

### 🔹 Phase 5 – Release 1.0 (December)
- [ ] Finalize full documentation
- [ ] Stabilize unit test suite
- [ ] Publish to GitHub as stable version

---

## 📚 Useful Resources
- LAPACK/BLAS documentation
- FFTW Fortran bindings
- “Numerical Recipes” & “Matrix Computations” as references
