Matrix Multiplication (Naive vs Strassen Algorithm)
This project implements and compares two matrix multiplication algorithms in C++: the Naive method and the Strassen algorithm. It benchmarks their performance and saves results to CSV files for further analysis.

Features
  -Matrix class supporting:
    -Dynamic allocation using STL vector
    -Random matrix generation with uniform distribution
    -Basic arithmetic operations (+, -)
    -Matrix dimension padding to the nearest power of two
    -CSV export functionality
  -Naive and Strassen algorithms for matrix multiplication
  -Automatic result saving to the data/ directory
  -Execution time measurement through a high-resolution clock
  -Configurable matrix dimensions

Algorithms Overview

Naive Matrix Multiplication
  -Classical approach with time complexity O(n<sup>3</sup>)
  -Iterates through all element combinations using three nested loops

Strassen Matrix Multiplication
  -Divide-and-conquer method with complexity O(n<sup>log⁡<sub>2</sub>7</sup>) ≈ O(n<sup>2.81</sup>)
  -Recursively splits matrices into quadrants
  -Reduces the number of multiplications from 8 to 7 per recursion step
  -Automatically switches to Naive algorithm for small matrices (size ≤ 64)

Clone:
```
git clone https://github.com/Pomme65123/Matrix_Multiplication.git
cd MatrixMultiplication
```
Compile:
```
g++ -std=c++17 main.cpp -o matrix_mult
```
