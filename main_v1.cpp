//NOTE NEVER EVER EVER TRY AND DEBUG WITH A MATRIX GENERATOR INSIDE THE TEST CASE

#include <chrono>
#include <iostream>
#include <random>
#include <vector>
#include <utility>
#include <cassert>
#include <cmath>
#include <limits>
#include <fstream>
#include <iomanip>


using namespace std;

class Matrix {
private:
	vector<vector<double>> data;

public:
	Matrix() : data(0, vector<double>(0)) {}
	Matrix(size_t rows, size_t columns) : data(rows, vector<double>(columns, 0.0)) {}
	size_t numRows() const {return data.size();}
	size_t numColumns() const {return data[0].size();}

	vector<vector<double>> &getData() {return data;}
	const vector<vector<double>> &getData() const {return data;}

	double &operator()(size_t i, size_t j) {
        return data[i][j];
    }

    const double &operator()(size_t i, size_t j) const {
        return data[i][j];
    }

	//Fills the matrix with random values I got from cppreference
	void randomMatrix() {
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<double> dis(0, 10);

		for (size_t i = 0; i < numRows(); i++) {
			for (size_t j = 0; j < numColumns(); j++) {
			  data[i][j] = dis(gen);
			}
		}
	}

	//Prints a matrix
	void printMatrix() const {
		cout << "Matrix: " << endl;
		for (size_t i = 0; i < numRows(); i++) {
			for (size_t j = 0; j < numColumns(); j++) {
				cout << "Position: [" << i << "][" << j << "] = " << data[i][j] << endl;
			}
		}
	}

	//Checks if it's a square matrix
	bool checkSquare() const {
		if (numRows() != numColumns()) {
			return 0;
		}
		return 1;
	}

	//Padding to square matrices
	static pair<Matrix,Matrix> padding(const Matrix &A, const Matrix &B) {
		size_t row_A = A.numRows();
		size_t column_A = A.numColumns();
		size_t column_B = B.numColumns();

		size_t newSize = max(row_A, max(column_A, column_B));
		newSize = static_cast<size_t> (pow(2, ceil(log2(newSize))));
		//cout << "The Size of the Matrix is: " << newSize <<endl;

		Matrix padded_A(newSize, newSize);
		Matrix padded_B(newSize, newSize);

		for (size_t i = 0; i < row_A; i++) {
			for (size_t j = 0; j < column_A; j++) {
				padded_A(i,j) = A(i,j);
			}
		}

		for (size_t i = 0; i < column_A; i++) {
			for (size_t j = 0; j < column_B; j++) {
				padded_B(i,j) = B(i,j);
			}
		}

		return make_pair(padded_A, padded_B);
	}

	Matrix operator-(const Matrix &other) const {
    	if (numRows() != other.numRows() || numColumns() != other.numColumns()) {
			throw std::invalid_argument("Matrix dimensions must match for subtraction.");
		}
		Matrix result(numRows(), numColumns());
		for (size_t i = 0; i < numRows(); ++i) {
			for (size_t j = 0; j < numColumns(); ++j) {
				result(i, j) = (*this)(i, j) - other(i, j);
			}
		}
		return result;
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

		for (size_t i = 0; i < numRows(); i++) {
			for (size_t j = 0; j < numColumns(); j++) {
				file << (*this)(i, j);
				if (j != numColumns() - 1) {
					file << ",";
				}
			}
			file << "\n";
		}
		file.close();
	}
};



class MatrixAlgorithms {
public:
	//Uses the Naive Algorithm for matrix multiplication
	static Matrix NaiveAlgorithm(	const Matrix &A, 
									const Matrix &B) {

		//https://en.wikipedia.org/wiki/Matrix_multiplication

		if (A.getData().empty() || B.getData().empty() || A.numColumns() != B.numRows()) {
		cout << "Invalid matrix multiplication" << endl;
			return Matrix(0,0);
		}

		size_t rowsA = A.numRows();
		size_t colsA = A.numColumns();
		size_t colsB = B.numColumns();

		Matrix result(rowsA, colsB);

		for (size_t i = 0; i < rowsA; i++) {
			for (size_t k = 0; k < colsA; k++) {
				double temp = A(i,k);
			  for (size_t j = 0; j < colsB; j++)
			    result(i,j) += temp * B(k,j);
			}
		}

		return result;
	}

	//Uses the Strassen Algorithm for matrix multiplication
	static Matrix StrassenAlgorithm(const Matrix &A, const Matrix &B) {

		//https://www.cise.ufl.edu/~sahni/papers/strassen.pdf
		//matrix size > 64

		if (A.getData().empty() || B.getData().empty() || A.numColumns() != B.numRows()) {
			cout << "Invalid matrix multiplication" << endl;
			return Matrix(0,0);
		}

		auto [padded_A, padded_B] = Matrix::padding(A, B);

		// cout << "BEGIN THE PADDING" << endl;
		// padded_A.printMatrix();
		// padded_B.printMatrix();

		size_t n = padded_A.numRows();

		if (n <= 64) {
			Matrix fullResult = NaiveAlgorithm(padded_A, padded_B);
			Matrix finalResult(A.numRows(), B.numColumns());
			for (size_t i = 0; i < A.numRows(); i++) {
				for (size_t j = 0; j < B.numColumns(); j++) {
					finalResult(i, j) = fullResult(i, j);
				}
			}
			return finalResult;
		}
	
		//LAMBDAS
		auto subMatrix = [](const Matrix &M, size_t row, size_t column, size_t size) {
			assert(row + size <= M.numRows());
			assert(column + size <= M.numColumns());
			Matrix result(size, size);
			for (size_t i = 0; i < size; i++) {
				for (size_t j = 0; j < size; j++) {
					result(i,j) = M(row + i, column + j);
				}
			}
			return result;
		};

		auto add = [](const Matrix &M_1, const Matrix &M_2) {
			//numRows() = numColumns() but for clarity
			Matrix result(M_1.numRows(), M_1.numColumns());
			for (size_t i = 0; i < M_1.numRows(); i++) {
				for (size_t j = 0; j < M_1.numColumns(); j++) {
					result(i,j) = M_1(i,j) + M_2(i,j);
				}
			}
			return result;
		};

		auto subtract = [](const Matrix &M_1, const Matrix & M_2) {
			Matrix result(M_1.numRows(), M_1.numColumns());
			for (size_t i = 0; i < M_1.numRows(); i++) {
				for (size_t j = 0; j < M_1.numColumns(); j++) {
					result(i,j) = M_1(i,j) - M_2(i,j);
				}
			}
			return result;
		};

		//END LAMBDAS

		size_t k = n/2;

		Matrix A11 = subMatrix(padded_A, 0, 0, k);
		Matrix A12 = subMatrix(padded_A, 0, k, k);
		Matrix A21 = subMatrix(padded_A, k, 0, k);
		Matrix A22 = subMatrix(padded_A, k, k, k);
		
		Matrix B11 = subMatrix(padded_B, 0, 0, k);
		Matrix B12 = subMatrix(padded_B, 0, k, k);
		Matrix B21 = subMatrix(padded_B, k, 0, k);
		Matrix B22 = subMatrix(padded_B, k, k, k);

		//DEGBURRING
		// A.saveCSV("A_data.csv");
		// B.saveCSV("B_data.csv");
		
		// padded_A.saveCSV("padded_A_data.csv");
		// padded_B.saveCSV("padded_B_data.csv");

		// A11.saveCSV("A11_data.csv");
		// A12.saveCSV("A12_data.csv");
		// A21.saveCSV("A21_data.csv");
		// A22.saveCSV("A22_data.csv");
		// B11.saveCSV("B11_data.csv");
		// B12.saveCSV("B12_data.csv");
		// B21.saveCSV("B21_data.csv");
		// B22.saveCSV("B22_data.csv");
		//DEGBURRING

		Matrix M1 = StrassenAlgorithm(add(A11, A22), add(B11, B22));
		Matrix M2 = StrassenAlgorithm(add(A21, A22), B11);
		Matrix M3 = StrassenAlgorithm(A11, subtract(B12, B22));
		Matrix M4 = StrassenAlgorithm(A22, subtract(B21, B11));
		Matrix M5 = StrassenAlgorithm(add(A11, A12), B22);
		Matrix M6 = StrassenAlgorithm(subtract(A21, A11), add(B11, B12));
		Matrix M7 = StrassenAlgorithm(subtract(A12, A22), add(B21, B22));

		Matrix C11 = add(subtract(add(M1, M4), M5) , M7);
		Matrix C12 = add(M3, M5);
		Matrix C21 = add(M2, M4);
		Matrix C22 = add(add(subtract(M1, M2), M3), M6);

		//Again this must be a squre matrix
		Matrix result(padded_A.numRows(), padded_A.numColumns());

		for(size_t i = 0; i < k; i++) {
			for (size_t j = 0; j < k; j++) {
				result(i,j) = C11(i,j);
				result(i, j + k) = C12(i, j);
				result(i + k, j) = C21(i, j);
				result(i + k, j + k) = C22(i, j);
			}
		}

		Matrix finalResult(A.numRows(), B.numColumns());
		for (size_t i = 0; i < finalResult.numRows(); i++) {
			for (size_t j = 0; j < finalResult.numColumns(); j++) {
				finalResult(i, j) = result(i, j);
			}
		}

		return finalResult; // PLACEHOLDER
	}

};



class MatrixTests {
public:
	//Tests the Naive Algorithm
	static Matrix NaiveAlgorithmTest(const Matrix &A, const Matrix &B) {

		cout << "Begin: Naive" << endl;

		// A.printMatrix();
		// B.printMatrix();

		Matrix result = MatrixAlgorithms::NaiveAlgorithm(A, B);

		result.saveCSV("result_Naive_Alorithm.csv");

		return result;
	}

	//Tests the Strassen Algorithm
	static Matrix StrassenAlgorithmTest(const Matrix &A, const Matrix &B) {
	
		cout << "Begin: Strassen" << endl;

		// A.printMatrix();
		// B.printMatrix();

		Matrix result = MatrixAlgorithms::StrassenAlgorithm(A, B);

		result.saveCSV("result_Strassen_Alorithm.csv");

		return result;
	}

};




//Displays program execution time
template <typename Func>
void measureExecutionTime(Func algorithmTest, const Matrix &matrix_1, const Matrix &matrix_2) {

	auto start = std::chrono::high_resolution_clock::now();
	algorithmTest(matrix_1, matrix_2);
	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    cout << "Total Duration: " << duration.count() << "ms" << endl;
}


int main() {
	size_t row_1 = 512;
	size_t column_1 = 512;
	size_t row_2 = 512;
	size_t column_2 = 512;

	Matrix matrix_1(row_1, column_1);
	Matrix matrix_2(row_2, column_2);

	matrix_1.randomMatrix();
	matrix_2.randomMatrix();

	Matrix Naive_Matrix = MatrixTests::NaiveAlgorithmTest(matrix_1, matrix_2);
	Matrix Strassen_Matrix = MatrixTests::StrassenAlgorithmTest(matrix_1, matrix_2);

	Naive_Matrix.saveCSV("Naive_data.csv");
	Strassen_Matrix.saveCSV("Strassen_data.csv");

	Matrix difference = Naive_Matrix - Strassen_Matrix;
	difference.printMatrix();


	//Naive_Matrix.printMatrix();
	//Strassen_Matrix.printMatrix();

	// measureExecutionTime(MatrixTests::NaiveAlgorithmTest, matrix_1, matrix_2);
	// measureExecutionTime(MatrixTests::StrassenAlgorithmTest, matrix_1, matrix_2);

	//cout << (Naive_Matrix.isEqual(Strassen_Matrix)) << endl;
	//Matrix difference = Naive_Matrix - Strassen_Matrix;
	//difference.printMatrix();

	return 0;

}