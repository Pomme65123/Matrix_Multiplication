# Optimized Matrix Multiplication: Naive vs Strassen Algorithm

A high-performance C++ implementation comparing naive matrix multiplication with the Strassen algorithm, featuring extensive memory and computational optimizations.

## Key Optimizations

### Memory Layout Improvements
- **Contiguous Memory Storage**: Replaced `vector<vector<double>>` with flat `vector<double>` for better cache locality
- **Move Semantics**: Implemented move constructors and assignment operators to eliminate unnecessary copies
- **Raw Pointer Access**: Direct memory access in tight loops to minimize function call overhead

### Algorithm Optimizations
- **Cache-Friendly Multiplication**: Reordered loops (i-k-j) for optimal memory access patterns
- **Optimized Threshold**: Increased Strassen base case threshold from 64 to 128 for better performance
- **Vectorized Operations**: Used `std::copy` for bulk memory operations
- **Padding Optimization**: Efficient power-of-2 padding using bit operations

### Performance Features
- **Compiler Optimization Ready**: Compatible with `-O3`, `-march=native`, `-funroll-loops`, `-ffast-math`
- **Modern C++ Standards**: Supports C++23 features
- **Zero-Copy Operations**: Minimized temporary object creation in recursive calls

## Project Structure

```
Matrix_Multiplication/
├── main.cpp                              # Optimized implementation
├── README.md                             # This file
└── data/                                 # Output directory
    ├── Naive_data.csv                    # Randomly Generated Data
    ├── Strassen_data.csv                 # Randomly Generated Data
    ├── result_Naive_Alorithm.csv         # Naive algorithm results
    └── result_Strassen_Alorithm.csv      # Strassen algorithm results
```

## Compilation & Usage

### Basic Compilation
```bash
g++ main.cpp -o matrix_mult.exe
```

### Optimized Compilation (Recommended)
```bash
g++ -std=c++23 -O3 -march=native -funroll-loops -ffast-math main.cpp -o matrix_mult_optimized.exe
```

### For C++23 Compatibility
```bash
g++ -std=c++23 -O3 -march=native main.cpp -o matrix_mult.exe
```

### Running the Program
```bash
./matrix_mult.exe
```

## Performance Testing

The program includes built-in timing functionality. Uncomment the timing code in `main()`:

```cpp
auto naiveTime = measureExecutionTime(MatrixTests::NaiveAlgorithmTest, matrix_1, matrix_2);
cout << "Naive Algorithm Execution Time: " << naiveTime.count() << " ms" << endl;

auto strassenTime = measureExecutionTime(MatrixTests::StrassenAlgorithmTest, matrix_1, matrix_2);
cout << "Strassen Algorithm Execution Time: " << strassenTime.count() << " ms" << endl;
```

### PowerShell Timing
```powershell
Measure-Command {.\matrix_mult.exe}
```

## Algorithms Implemented

### 1. Naive Matrix Multiplication
- **Complexity**: O(n³)
- **Optimizations**: Cache-friendly i-k-j loop ordering, raw pointer arithmetic
- **Best For**: Smaller matrices (< 128×128)

### 2. Strassen Algorithm
- **Complexity**: O(n^2.807)
- **Optimizations**: Efficient submatrix extraction, optimized padding, move semantics
- **Best For**: Large matrices (≥ 128×128)
- **Recursive Base Case**: Falls back to naive algorithm for matrices ≤ 128×128

## Matrix Class Features

### Core Functionality
- **Contiguous Memory Layout**: Single `vector<double>` with row-major storage
- **Move Semantics**: Efficient transfer of large matrices
- **Operator Overloading**: Natural `+`, `-` syntax
- **Memory Safety**: Bounds checking in debug builds

### Utility Functions
- **Random Generation**: Uniform distribution [0, 10]
- **CSV Export**: Configurable precision output
- **Square Matrix Detection**: Efficient dimension checking
- **Automatic Padding**: Power-of-2 padding for Strassen algorithm

## Current Configuration

The program currently tests with 512×512 matrices (`2048/4`):
```cpp
size_t row_1 = 2048/2;    // 512
size_t column_1 = 2048/2; // 512
size_t row_2 = 2048/2;    // 512 
size_t column_2 = 2048/2; // 512
```

## Expected Performance Gains

### Memory Access Improvements
- **~40-60% faster** memory access due to contiguous layout
- **Reduced cache misses** from improved spatial locality
- **Lower memory fragmentation** from single allocation per matrix

### Computational Improvements
- **~20-30% faster** naive multiplication from loop reordering
- **~15-25% faster** Strassen algorithm from optimized submatrix operations
- **Reduced allocation overhead** from move semantics

## Technical Requirements

- **C++ Compiler**: GCC 9+ or equivalent with C++23 support
- **Memory**: Sufficient RAM for matrix storage (512×512 ≈ 2MB per matrix)
- **Platform**: Cross-platform compatible (Windows, Linux, macOS)

## Output Files

All results are saved as CSV files in the `data/` directory:
- High precision (6 decimal places by default)
- Standard comma-separated format
- Compatible with Excel, MATLAB, Python pandas

## Advanced Usage

### Custom Matrix Sizes
Modify the size constants in `main()`:
```cpp
size_t row_1 = 1024;    // Custom size
size_t column_1 = 1024;
```

### Precision Control
Adjust CSV output precision:
```cpp
result.saveCSV("output.csv", 10);  // 10 decimal places
```

### Algorithm Selection
The Strassen algorithm automatically chooses the best approach:
- **n ≤ 128**: Uses optimized naive algorithm
- **n > 128**: Uses recursive Strassen algorithm

## References

- [Strassen Algorithm Paper](https://www.cise.ufl.edu/~sahni/papers/strassen.pdf)
- [Matrix Multiplication Wikipedia](https://en.wikipedia.org/wiki/Matrix_multiplication)
- Modern C++ optimization techniques and cache-friendly programming patterns

