# Matrix Library Documentation

## Overview

This is a comprehensive C++ matrix algebra library implementing algorithms from MIT OpenCourseWare 18.06 (Linear Algebra). The library supports matrix arithmetic, least squares fitting, matrix compression, linear solvers, factorization, and advanced matrix operations including FFT and eigenvalue computations.

---

## Table of Contents

1. [Complex Number Class](#complex-number-class)
2. [Matrix Class](#matrix-class)
3. [Matrix Types & Compression](#matrix-types--compression)
4. [Core Operations](#core-operations)
5. [Advanced Algorithms](#advanced-algorithms)
6. [Utility Functions](#utility-functions)

---

## Complex Number Class

### File: `complex.h` / `complex.cpp`

The `complex` class provides complex number arithmetic with mathematical operations.

#### Class Definition

```cpp
class complex {
private:
    long double re, im;  // Real and imaginary parts
public:
    // Constructors
    complex(long double _re=0, long double _im=0);
    complex(int);
    complex(const complex&);
    
    // Getters
    long double get_re(void) const;
    long double get_im(void) const;
    
    // Arithmetic operators
    complex operator+(const complex&) const;
    complex operator-(const complex&) const;
    complex operator*(const complex&) const;
    complex operator/(const complex&) const;
    
    // Scalar operations
    complex operator*(const long double&) const;
    complex operator/(const long double&) const;
    
    // Compound assignment operators
    void operator+=(const complex&);
    void operator-=(const complex&);
    void operator*=(const complex&);
    void operator/=(const complex&);
    void operator^=(const long double& power);
    
    // Power operation
    complex operator^(const long double& power) const;
    
    // Comparison operators
    bool operator>(const complex&) const;
    bool operator<(const complex&) const;
    bool operator==(const complex&) const;
    bool operator!=(const complex&) const;
    
    // Polar form
    long double theta(void) const;  // Returns angle in radians
    
    // Stream output
    friend std::ostream& operator<<(std::ostream& os, const complex& c);
};
```

#### Key Functions

| Function | Purpose |
|----------|---------|
| `abs(const complex&)` | Magnitude of complex number |
| `conjugate(const complex&)` | Complex conjugate |
| `pow(const complex&, int)` | Integer power of complex number |
| `to_string(complex&)` | Convert to string representation |

#### Usage Example

```cpp
complex c1(3, 4);        // 3 + 4j
complex c2(1, 2);        // 1 + 2j
complex c3 = c1 + c2;    // Addition
complex c4 = c1 * c2;    // Multiplication
long double mag = abs(c1);  // Magnitude
```

---

## Matrix Class

### File: `matrix_algebra.h` / `matrix_algebra.cpp`

Template-based matrix class supporting multiple data types (int, float, double, complex, long double).

### Template Declaration

```cpp
template<typename DataType>
class matrix {
    // ... implementation
};
```

### Constructors

```cpp
matrix();                                           // Empty matrix
matrix(int rows, int cols, DataType value=no_init, bool compressed=false);
matrix(int rows, int cols, DataType* init_array, int array_size);
matrix(const matrix&);                             // Copy constructor
```

### Dimensions & Getters

```cpp
int get_rows() const;      // Number of rows
int get_cols() const;      // Number of columns
int get_size() const;      // Number of elements stored
int get_type() const;      // Matrix type (general, utri, ltri, etc.)
```

### Element Access

```cpp
DataType& at(int row, int col);           // Mutable access
const DataType& at(int row, int col) const; // Const access
bool is_valid_index(int row, int col) const;

// For all matrix types: general, utri (upper triangular), 
// ltri (lower triangular), diagonal, identity, symmetric, 
// anti-symmetric, constant, orthonormal
```

---

## Matrix Types & Compression

The library optimizes storage by detecting and compressing special matrix types.

### Matrix Types Enum

```cpp
enum {
    general = 0,        // Full matrix, no special structure
    utri,              // Upper triangular
    ltri,              // Lower triangular
    diagonal,          // Diagonal matrix
    symmetric,         // A = A^T
    anti_symmetric,    // A = -A^T (skew-symmetric)
    constant,          // All elements are the same
    iden,              // Identity matrix
    orthonormal        // Columns form orthonormal set
};
```

### Compression Methods

```cpp
void compress(void);                    // Auto-detect and compress
void decompress(void);                  // Restore to general form
void fill_features(double check_tol=1e-6);  // Analyze matrix properties

// Type detection functions
bool is_symmetric(void) const;
bool is_anti_symmetric(void) const;
bool is_diagonal(void) const;
bool is_identity(void) const;
bool is_upper_tri(void) const;
bool is_lower_tri(void) const;
bool is_zero(void) const;
bool is_scalar(void) const;
bool is_idempotent(void) const;      // A² = A
bool is_nilpotent(void) const;       // A^k = 0 for some k
bool is_involutory(void) const;      // A² = I
bool is_orthogonal(void) const;      // A·A^T = I
bool is_positive_definite(void) const;
bool is_independent(void) const;     // Column vectors linearly independent
```

### Storage Efficiency

- **General**: Full storage (m×n elements)
- **Triangular**: Only non-zero elements stored + mapping arrays
- **Symmetric/Anti-symmetric**: Upper triangle only + formula for lower
- **Diagonal**: Only diagonal elements
- **Identity**: Single element (constant 1)
- **Constant**: Single element (repeated value)

---

## Core Operations

### Basic Arithmetic

```cpp
matrix operator+(const matrix&) const;  // Addition
matrix operator-(const matrix&) const;  // Subtraction
matrix operator*(const matrix&) const;  // Matrix multiplication
matrix operator/(const matrix&) const;  // A/B = A·B^(-1)
matrix operator*(DataType scalar) const; // Scalar multiplication
matrix operator^(const DataType& power) const; // Matrix power A^k
```

### Matrix Properties

```cpp
bool is_square(void) const;
bool same_shape(const matrix&) const;
bool operator==(const matrix&) const;

DataType det(void) const;              // Determinant (square only)
DataType trace(void) const;            // Sum of diagonal (square only)
int rank(void) const;                  // Rank (via Gaussian elimination)
```

### Transpose & Adjoint

```cpp
matrix transpose(void) const;          // M^T (conjugate for complex)
```

### Vector Operations

```cpp
DataType dot(const matrix&) const;     // Inner product
DataType norm2(void) const;            // Euclidean norm (L2)
DataType length(void) const;           // Alias for norm2()
DataType theta(matrix&) const;         // Angle between vectors (in degrees)
bool is_orthogonal(const matrix&) const;  // Perpendicular check
bool is_parallel(const matrix&) const;    // Parallel check
DataType col_length(const int& col_i) const; // Length of column
```

### Display & String Conversion

```cpp
void show(void) const;                 // Print matrix to console
string mat_to_string(void) const;      // Convert to string
void fill(DataType value);             // Fill all elements with value
void set_identity(void);               // Convert to identity matrix
```

---

## Advanced Algorithms

### Gaussian Elimination

```cpp
// Downward elimination → Upper triangular form
matrix gauss_down(matrix<int>* pivots_indices=NULL,
                  int pivots_locations=new_locations) const;

// Upward elimination → Lower triangular form
matrix gauss_up(matrix<int>* pivots_indices=NULL) const;

// Row operations
void row_axpy(DataType alpha, int x_row, int y_row);  // row_y += α·row_x
bool switch_rows(int row1, int row2);                  // Swap rows
```

### LU Factorization

```cpp
void lu_fact(matrix& lower_fact, matrix& permutation, matrix& upper_fact) const;
// Decomposes A = P·L·U where:
//   P = permutation matrix (row exchanges)
//   L = lower triangular matrix
//   U = upper triangular matrix
```

### QR Factorization

```cpp
void qr_fact(matrix& q, matrix& r) const;
// Decomposes A = Q·R where:
//   Q = orthogonal matrix
//   R = upper triangular matrix
```

### SVD (Singular Value Decomposition)

```cpp
void svd(matrix& u, matrix& s, matrix& vt) const;
// Decomposes A = U·S·V^T where:
//   U = left singular vectors
//   S = diagonal singular values
//   V^T = right singular vectors (transposed)
```

### Gram-Schmidt Orthonormalization

```cpp
matrix gram_shmidt(void) const;
// Converts column vectors into orthonormal basis
// Returns orthonormal matrix
```

### Linear Solver

```cpp
matrix solve(void) const;
// Solves Ax = b for an appended matrix [A|b]
// Uses Gaussian elimination + back substitution
// Returns solution vector x
```

### Substitution Methods

```cpp
DataType back_sub(int selected_row, const matrix& solution_matrix) const;
// Back substitution on upper triangular system

DataType fwd_sub(int selected_row, const matrix& solution_matrix) const;
// Forward substitution on lower triangular system
```

### Matrix Inversion

```cpp
matrix inverse(void) const;
// Computes M^(-1) for square invertible matrices
// Returns identity matrix if orthogonal (efficient)
// Uses Gauss-Jordan elimination otherwise
```

### Eigenvalue & Eigenvector Computations

```cpp
matrix eigen_values(const int& max_iteration,
                   const double& min_diff=check_tolerance) const;
// Computes eigenvalues using QR factorization
// Returns column matrix of eigenvalues

matrix eigen_vectors(const matrix& eigen_values) const;
// Computes eigenvectors given eigenvalues
// Solves (A - λI)v = 0 for each eigenvalue λ

DataType cofactor(int row_i, int col_i) const;
matrix cofactors(void) const;
matrix SubLambdaI(DataType lambda) const;  // Returns A - λI
```

---

## Linear Algebra Transformations

### Least Squares Fitting

```cpp
matrix fit_least_squares(const matrix& output) const;
// Fits input→output linear transformation
// Solves: A^T·A·x = A^T·b
// Returns least squares solution x

matrix projection(void) const;
// Projection matrix P = A·(A^T·A)^(-1)·A^T
// Projects any vector b onto column space of A
```

### Basis & Dimension Operations

```cpp
int dim(void) const;                        // Dimension of column space
int dim_null_cols(void) const;              // Dimension of null space (columns)
int dim_null_rows(void) const;              // Dimension of null space (rows)
bool is_basis(int dimension) const;         // Check if vectors form basis

matrix basis_cols(void) const;              // Basis for column space
matrix basis_rows(void) const;              // Basis for row space
matrix null_cols(void) const;               // Null space basis (columns)
matrix null_rows(matrix* elementary=NULL) const;  // Null space basis (rows)

matrix extract_col(int index) const;        // Extract column as matrix
```

### Row Reduced Echelon Form

```cpp
matrix rref(matrix<int>* pivots_indices=NULL) const;
// Row Reduced Echelon Form
// Returns matrix in RREF with pivot locations
// Used for determining rank, null space, basis
```

### Pivot Analysis

```cpp
int is_pivot(int row_index, int col_index);
int is_pivot_up(int row_index, int col_index);
matrix get_pivots(matrix<int>* pivots_locations=NULL) const;
matrix<int> get_pindex(void);
void fix_pivots(void);
```

---

## Signal Processing & FFT

### Fast Fourier Transform

```cpp
matrix fft_col(int dimension) const;
// FFT on single column
// Requires power-of-2 dimension

matrix fft(void) const;
// FFT on entire matrix (column-wise)

matrix split(int half) const;
// Split matrix into halves (upper/lower)
// Used internally by FFT

// Helper functions
matrix<complex> fourier_diagonal(int dimension, int n);
matrix<complex> fourier_mat(int dimension);
```

---

## Matrix Reshaping & Manipulation

### Matrix Composition

```cpp
matrix append_cols(const matrix& src) const;
// Append source matrix columns to right
// Returns [A | B]

matrix append_rows(const matrix& src) const;
// Append source matrix rows below
// Returns [A]
//         [B]

matrix arrange(const matrix<int>& seq) const;
// Rearrange rows according to sequence
```

### Resizing & Region Operations

```cpp
matrix resize(int wanted_rows=get_rows(),
              int wanted_cols=get_cols(),
              DataType padding_value=0) const;
// Resize matrix with optional padding

void at_quarter(int region, const matrix& input);
// Place input matrix in specific quadrant
// region: upper_left, lower_left, lower_right, upper_right
```

### Filtering

```cpp
void filter(double filter_tolerance=check_tolerance);
// Remove elements below tolerance threshold
// Effectively zeros out small numerical errors
```

---

## Utility Functions

### Template Helpers

```cpp
template <typename DataType>
matrix<DataType> identity(int n);
// Creates n×n identity matrix

template <typename DataType>
matrix<DataType> rand(int rows, int cols, int max_val=INT_MAX);
// Creates matrix with random elements
```

### Global Functions

```cpp
template<typename DataType>
DataType* get_vec(int r, int c);
// Allocate memory for vector

template<typename DataType>
void fill_vec(DataType* vec, int size, DataType val=-1);
// Fill vector with value

template<typename DataType>
void copy_vec(DataType*& dest, const DataType*& src, int size);
// Copy vector
```

### Complex-Specific

```cpp
complex conjugate(const complex& val);
// Complex conjugate

// Overloads for all scalar types (return identity for non-complex)
float conjugate(const float& val);
double conjugate(const double& val);
long double conjugate(const long double& val);
// ... etc for all integer types
```

---

## Constants & Configuration

### Tolerance Values

```cpp
const long double tolerance = 1e-32;
// Used in: inverse, Gram-Schmidt, Gaussian elimination
// For numerical stability in precision calculations

const long double check_tolerance = 1e-6;
// Used in: is_identity, is_zero, feature detection
// More lenient for structural checks
```

### Conversion Constants

```cpp
const long double to_deg = 180/M_PI;   // Radians to degrees
const long double to_rad = M_PI/180;   // Degrees to radians
```

---

## Error Handling

The library uses return values to signal errors:

```cpp
string shape_error = "\nmatrices aren't the same shape default garbage value is -1\n";
string square_error = "\nmatrix must be square to perform this operation default garbage value is -1\n";
string uninit_error = "\nmatrix isn't initialized yet\n";
```

### Common Error Return Values

- **-1**: General error indicator in scalar return types
- **matrix(1,1,-1)**: Error matrix for matrix return types
- **Empty matrix**: Failed operations may return empty matrices

---

## Performance Characteristics

### Time Complexity

| Operation | Complexity |
|-----------|-----------|
| Addition/Subtraction | O(mn) |
| Scalar Multiplication | O(mn) |
| Matrix Multiplication | O(mnp) for m×n × n×p |
| Transpose | O(mn) |
| Determinant (Gaussian Elim.) | O(n³) |
| Inverse (Gauss-Jordan) | O(n³) |
| LU Factorization | O(n³) |
| QR Factorization | O(mn²) |
| SVD | O(mn²) |
| Rank | O(n³) |
| FFT | O(n log n) |

### Space Complexity

- **Uncompressed**: O(mn)
- **Compressed (special types)**: O(k) where k << mn
  - Diagonal: O(n)
  - Identity: O(1)
  - Constant: O(1)
  - Symmetric: O(n²/2 + n)

---

## Example Usage

```cpp
#include "matrix_algebra.h"
#include <iostream>

using namespace std;

int main() {
    // Create matrices
    matrix<double> A(3, 3);
    matrix<double> B(3, 3);
    
    // Fill with data
    A.fill(2.0);
    B.set_identity();
    
    // Basic operations
    matrix<double> C = A + B;      // Addition
    matrix<double> D = A * B;      // Multiplication
    double det_A = A.det();        // Determinant
    
    // Advanced operations
    matrix<double> inv_A = A.inverse();  // Inverse
    int r = A.rank();                     // Rank
    
    // Check properties
    if (A.is_symmetric()) {
        cout << "A is symmetric\n";
    }
    
    // Factorizations
    matrix<double> L, U, P;
    A.lu_fact(L, P, U);  // LU decomposition
    
    // Linear solver
    matrix<double> b(3, 1);
    b.fill(1.0);
    matrix<double> Av_b = A.append_cols(b);
    matrix<double> x = Av_b.solve();      // Solve Ax = b
    
    // Display results
    A.show();
    
    return 0;
}
```

---

## Notes

- All matrix indices are 0-based
- Matrices can be reused with the assignment operator
- Compressed matrices work transparently with uncompressed ones in operations
- Complex number support for all operations
- Numerical stability is prioritized with configurable tolerance values
