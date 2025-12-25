#include <chrono>
#include <iostream>
#include <random>
#include <vector>
#include <utility>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <functional>
#include <algorithm>
#include <algorithm>
using namespace std;

class Matrix {
public:
	vector<double> data;
	size_t rows, cols;
	Matrix() : data(), rows(0), cols(0) {}
	Matrix(size_t rows, size_t columns) : data(rows * columns, 0.0), rows(rows), cols(columns) {}
	
	Matrix(Matrix&& other) noexcept : data(std::move(other.data)), rows(other.rows), cols(other.cols) {
		other.rows = other.cols = 0;
	}
	
	Matrix& operator=(Matrix&& other) noexcept {
		if (this != &other) {
			data = std::move(other.data);
			rows = other.rows;
			cols = other.cols;
			other.rows = other.cols = 0;
		}
		return *this;
	}
	
	size_t numRows() const noexcept { return rows; }
	size_t numColumns() const noexcept { return cols; }
	bool empty() const noexcept { return data.empty(); }

	double& operator()(size_t i, size_t j) noexcept {
		return data[i * cols + j];
	}

	const double& operator()(size_t i, size_t j) const noexcept {
		return data[i * cols + j];
	}
	
	double* rawData() noexcept { return data.data(); }
	const double* rawData() const noexcept { return data.data(); }

    Matrix operator+(const Matrix &other) const {
    	if (numRows() != other.numRows() || numColumns() != other.numColumns()) {
			throw std::invalid_argument("Matrix dimensions must match for addition.");
		}

		Matrix result(numRows(), numColumns());

		for (auto col{0uz}; col < numColumns(); col++) {
			for (auto row{0uz}; row < numRows(); row++) {
				result(row, col) = (*this)(row, col) + other(row, col);
			}
		}

		return result;
    }

	Matrix operator-(const Matrix &other) const {
    	if (numRows() != other.numRows() || numColumns() != other.numColumns()) {
			throw std::invalid_argument("Matrix dimensions must match for subtraction.");
		}

		Matrix result(numRows(), numColumns());

		for (auto col{0uz}; col < numColumns(); col++) {
			for (auto row{0uz}; row < numRows(); row++) {
				result(row, col) = (*this)(row, col) - other(row, col);
			}
		}

		return result;
	}

	friend ostream &operator<<(ostream &out, const Matrix &other);

	//Fills the matrix with random values I got from cppreference
	void randomMatrix() {
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<double> dis(0, 10);

		double* ptr = data.data();
		const size_t size = data.size();
		for (auto i{0uz}; i < size; ++i) {
			ptr[i] = dis(gen);
		}
	}

	//Checks if it's a square matrix
	bool checkSquare() const noexcept {
		return rows == cols;
	}

	//Padding to square matrices
	static pair<Matrix,Matrix> padding(const Matrix &A, const Matrix &B) {

		if (A.cols != B.rows) {
		    throw std::invalid_argument("Invalid matrix multiplication");
		}

		size_t row_A = A.rows;
		size_t column_A = A.cols;
		size_t column_B = B.cols;

		size_t newSize = std::max({row_A, column_A, column_B});
		// Find next power of 2
		newSize = static_cast<size_t>(1) << static_cast<size_t>(ceil(log2(newSize)));

		Matrix padded_A(newSize, newSize);
		Matrix padded_B(newSize, newSize);

		const double* srcA = A.data.data();
		double* dstA = padded_A.data.data();
		for (auto row{0uz}; row < row_A; row++) {
			std::copy(srcA + row * column_A, srcA + (row + 1) * column_A, 
			         dstA + row * newSize);
		}

		const double* srcB = B.data.data();
		double* dstB = padded_B.data.data();
		for (auto col{0uz}; col < column_A; col++) {
			std::copy(srcB + col * column_B, srcB + (col + 1) * column_B, 
			         dstB + col * newSize);
		}

		return make_pair(std::move(padded_A), std::move(padded_B));
	}

	void saveCSV(const string &filename, int precision = 6) const {
		string directory = "data/";
		string fullPath = directory + filename;
		ofstream file(fullPath);
		if(!file.is_open()) {
			cerr << "File Bad No Open" << endl;
			return;
		}

		file << fixed << setprecision(precision);
		const double* ptr = data.data();

		for (auto col{0uz}; col < cols; col++) {
			for (auto row{0uz}; row < rows; row++) {
				file << ptr[row * cols + col];
				if (col != cols - 1) {
					file << ",";
				}
			}
			file << "\n";
		}
		file.close();
	}
};



ostream &operator<<(ostream &out, const Matrix &other) {
	for (auto col{0uz}; col < other.cols; col++) {
		for (auto row{0uz}; row < other.rows; row++) {
			out << "Position: (" << row << "," << col << "): " << other(row,col) << endl;
		}
	}
	return out;
}



class MatrixAlgorithms {
public:
	//Uses the Naive Algorithm for matrix multiplication
	static Matrix NaiveAlgorithm(const Matrix &A, const Matrix &B) {

		//https://en.wikipedia.org/wiki/Matrix_multiplication

		if (A.empty() || B.empty() || A.cols != B.rows) {
			cout << "Invalid matrix multiplication" << endl;
			return Matrix(0,0);
		}

		const size_t rowsA = A.rows;
		const size_t colsA = A.cols;
		const size_t colsB = B.cols;

		Matrix result(rowsA, colsB);
		

		const double* aPtr = A.data.data();
		const double* bPtr = B.data.data();
		double* resPtr = result.data.data();

		for (auto col_A{0uz}; col_A < colsA; col_A++) {
			for (auto row{0uz}; row < rowsA; row++) {
				const double aik = aPtr[row * colsA + col_A];
				const double* bRow = bPtr + col_A * colsB;
				double* resRow = resPtr + row * colsB;
				
				for (auto col_B{0uz}; col_B < colsB; col_B++) {
					resRow[col_B] += aik * bRow[col_B];
				}
			}
		}

		return result;
	}

	//Uses the Strassen Algorithm for matrix multiplication
	static Matrix StrassenAlgorithm(const Matrix &A, const Matrix &B) {

		//https://www.cise.ufl.edu/~sahni/papers/strassen.pdf
		//matrix size > 64

		if (A.empty() || B.empty() || A.cols != B.rows) {
			cout << "Invalid matrix multiplication" << endl;
			return Matrix(0,0);
		}

		auto [padded_A, padded_B] = Matrix::padding(A, B);

		size_t size_padded = padded_A.rows;

		if (size_padded <= 128) {  // Isize_paddedcreased threshold for better performance
			Matrix fullResult = NaiveAlgorithm(padded_A, padded_B);
			Matrix finalResult(A.rows, B.cols);
			
			// Optimized copying with raw pointers
			const double* srcPtr = fullResult.data.data();
			double* dstPtr = finalResult.data.data();
			for (auto row{0uz}; row < A.rows; row++) {
				std::copy(srcPtr + row * fullResult.cols, 
				         srcPtr + row * fullResult.cols + B.cols,
				         dstPtr + row * B.cols);
			}
			return finalResult;
		}
	
		// Optimized subMatrix with move semantics
		auto subMatrix = [](const Matrix &M, size_t row, size_t column, size_t size) {
			Matrix result(size, size);
			const double* srcPtr = M.data.data();
			double* dstPtr = result.data.data();
			
			for (auto i{0uz}; i < size; i++) {
				std::copy(srcPtr + (row + i) * M.cols + column,
				         srcPtr + (row + i) * M.cols + column + size,
				         dstPtr + i * size);
			}
			return result;
		};

		size_t k = size_padded/2;

		Matrix A11 = subMatrix(padded_A, 0, 0, k);
		Matrix A12 = subMatrix(padded_A, 0, k, k);
		Matrix A21 = subMatrix(padded_A, k, 0, k);
		Matrix A22 = subMatrix(padded_A, k, k, k);
		
		Matrix B11 = subMatrix(padded_B, 0, 0, k);
		Matrix B12 = subMatrix(padded_B, 0, k, k);
		Matrix B21 = subMatrix(padded_B, k, 0, k);
		Matrix B22 = subMatrix(padded_B, k, k, k);

		Matrix M1 = StrassenAlgorithm(A11 + A22, B11 + B22);
		Matrix M2 = StrassenAlgorithm(A21 + A22, B11);
		Matrix M3 = StrassenAlgorithm(A11, B12 - B22);
		Matrix M4 = StrassenAlgorithm(A22, B21 - B11);
		Matrix M5 = StrassenAlgorithm(A11 + A12, B22);
		Matrix M6 = StrassenAlgorithm(A21 - A11, B11 + B12);
		Matrix M7 = StrassenAlgorithm(A12 - A22, B21 + B22);

		Matrix C11 = M1 + M4 - M5 + M7;
		Matrix C12 = M3 + M5;
		Matrix C21 = M2 + M4;
		Matrix C22 = M1 - M2 + M3 + M6;

		Matrix result(padded_A.rows, padded_A.cols);
		
		double* resultPtr = result.data.data();
		const double* c11Ptr = C11.data.data();
		const double* c12Ptr = C12.data.data();
		const double* c21Ptr = C21.data.data();
		const double* c22Ptr = C22.data.data();
		
		for(auto i{0uz}; i < k; ++i) {
			std::copy(c11Ptr + i * k, c11Ptr + (i + 1) * k, 
			         resultPtr + i * result.cols);
			std::copy(c12Ptr + i * k, c12Ptr + (i + 1) * k, 
			         resultPtr + i * result.cols + k);
			         
			std::copy(c21Ptr + i * k, c21Ptr + (i + 1) * k, 
			         resultPtr + (i + k) * result.cols);
			std::copy(c22Ptr + i * k, c22Ptr + (i + 1) * k, 
			         resultPtr + (i + k) * result.cols + k);
		}

		Matrix finalResult(A.rows, B.cols);
		const double* srcPtr = result.data.data();
		double* dstPtr = finalResult.data.data();
		for (auto row{0uz}; row < A.rows; row++) {
			std::copy(srcPtr + row * result.cols, 
			         srcPtr + row * result.cols + B.cols,
			         dstPtr + row * B.cols);
		}

		return finalResult;
	}

};



class MatrixTests {
public:
	//Tests the Naive Algorithm
	static Matrix NaiveAlgorithmTest(const Matrix &A, const Matrix &B) {

		cout << "Begin: Naive" << endl;

		Matrix result = MatrixAlgorithms::NaiveAlgorithm(A, B);

		result.saveCSV("result_Naive_Alorithm.csv");

		return result;
	}

	//Tests the Strassen Algorithm
	static Matrix StrassenAlgorithmTest(const Matrix &A, const Matrix &B) {
	
		cout << "Begin: Strassen" << endl;

		Matrix result = MatrixAlgorithms::StrassenAlgorithm(A, B);

		result.saveCSV("result_Strassen_Alorithm.csv");

		return result;
	}

};


//Displays program execution time
template <typename Func, typename... Args>
auto measureExecutionTime(Func&& algorithmTest, Args &&... args) {

	auto start = std::chrono::steady_clock::now();
    std::invoke(std::forward<Func>(algorithmTest), std::forward<Args>(args)...);
    auto end = std::chrono::steady_clock::now();
    
    return std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
}

int main() {
	size_t row_1 = 2048/2;
	size_t column_1 = 2048/2;
	size_t row_2 = 2048/2;
	size_t column_2 = 2048/2;

	Matrix matrix_1(row_1, column_1);
	Matrix matrix_2(row_2, column_2);

	matrix_1.randomMatrix();
	matrix_2.randomMatrix();

	Matrix Naive_Matrix = MatrixTests::NaiveAlgorithmTest(matrix_1, matrix_2);
	Matrix Strassen_Matrix = MatrixTests::StrassenAlgorithmTest(matrix_1, matrix_2);

	Naive_Matrix.saveCSV("Naive_data.csv");
	Strassen_Matrix.saveCSV("Strassen_data.csv");

	// auto naiveTime = measureExecutionTime(MatrixTests::NaiveAlgorithmTest, matrix_1, matrix_2);
    // cout << "Naive Algorithm Execution Time: " << naiveTime.count() << " ms" << endl;

    // auto strassenTime = measureExecutionTime(MatrixTests::StrassenAlgorithmTest, matrix_1, matrix_2);
    // cout << "Strassen Algorithm Execution Time: " << strassenTime.count() << " ms" << endl;

	return 0;

}