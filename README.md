# Matrix-Library (C++)

A comprehensive, **template-based C++ matrix algebra library** inspired by MIT OpenCourseWare **18.06 (Linear Algebra)**. It supports common matrix operations plus higher-level linear-algebra algorithms such as Gaussian elimination, LU/QR factorization, SVD, least-squares fitting, eigen computations, and FFT utilities. It also supports complex matrices.

## Features

### Matrix class (template)
- Template matrix type: `matrix<DataType>` (works with `int`, `float`, `double`, `long double`, and the included `complex`)
- **Basic operations**
  - Addition/subtraction: `operator+`, `operator-`
  - Matrix multiply/divide: `operator*`, `operator/` (division uses inverse)
  - Scalar multiply: `operator*(DataType)`
  - Power: `operator^(power)`
  - Transpose: `transpose()`
  - Dot product: `dot()`
- **Matrix properties**
  - `det()`, `trace()`, `rank()`
  - `is_square()`, `same_shape()`, `operator==`
  - Feature checks: `is_symmetric()`, `is_identity()`, `is_zero()`, `is_diagonal()`, `is_upper_tri()`, `is_lower_tri()`, `is_orthogonal()`, etc.
- **Solvers / factorizations / advanced algorithms**
  - Gaussian elimination: `gauss_down()`, `gauss_up()`, `rref()`
  - LU: `lu_fact(L, P, U)`
  - QR: `qr_fact(Q, R)`
  - SVD: `svd(U, S, VT)`
  - Inverse: `inverse()`
  - Linear system solve (augmented matrix): `solve()`
  - Least squares: `fit_least_squares(output)`
  - Projection matrix: `projection()`
  - Gram–Schmidt: `gram_shmidt()`
  - Eigenvalues / eigenvectors: `eigen_values(max_iteration, min_diff)`, `eigen_vectors(eigen_values)`
- **FFT**
  - `fft_col(dimension)`, `fft()`
  - Fourier helpers: `fourier_diagonal()`, `fourier_mat()`

### Matrix type detection & compression
The library can **detect special matrix structures** and store them more efficiently:
- Types include: `general`, `utri`, `ltri`, `diagonal`, `symmetric`, `anti_symmetric`, `constant`, `iden`, `orthonormal`
- Methods:
  - `fill_features(check_tol)`
  - `compress()`
  - `decompress()`

### Complex numbers
Includes a custom `complex` number type in `complex.h / complex.cpp` supporting:
- Arithmetic operators, scalar ops, power `^`, comparisons, polar angle `theta()`
- Helper functions like `abs()`, `conjugate()`, etc.
- Complex matrices are supported with their operations (hermitian inner product)
  
## Repository layout
- `matrix_algebra.h` / `matrix_algebra.cpp` — main matrix library
- `complex.h` / `complex.cpp` — complex number implementation
- `DOCUMENTATION.MD` — detailed API/feature documentation
- `DEMO.md` — demo link
- `README.md` — (this file)

## Getting started

### 1) Include the headers
```cpp
#include "matrix_algebra.h"
#include <iostream>
using namespace std;
```

### 2) Example usage
```cpp
#include "matrix_algebra.h"
#include <iostream>

using namespace std;

int main() {
    matrix<double> A(3, 3);
    matrix<double> B(3, 3);

    A.fill(2.0);
    B.set_identity();

    matrix<double> C = A + B;
    matrix<double> D = A * B;

    cout << "det(A) = " << A.det() << "\n";
    cout << "rank(A) = " << A.rank() << "\n";

    matrix<double> invA = A.inverse();

    // Solve Ax = b using appended matrix [A|b]
    matrix<double> b(3, 1);
    b.fill(1.0);
    matrix<double> Ab = A.append_cols(b);
    matrix<double> x = Ab.solve();

    A.show();
    return 0;
}
```

## Notes / behavior
- **Indexing is 0-based** (`at(row, col)`).
- Some operations return “error values” (e.g., `-1` or a `matrix(1,1,-1)`) when shapes are invalid or a square matrix is required.
- `compress()` can reduce storage for special matrix types; operations are intended to work with both compressed and uncompressed matrices.

## Demo
See `DEMO.md` for the YouTube demo link.

## Documentation
For the full list of functions and explanations, see `DOCUMENTATION.MD`.
