#include <chrono>
#include <iostream>
#include <random>
#include <vector>

using namespace std;

class Matrix {
private:
	vector<vector<double>> data;

public:
	Matrix(size_t rows, size_t columns) : data(rows, vector<double>(columns, 0.0)) {}
	size_t numRows() const {return data.size();}
	size_t numColumns() const {return data[0].size();}

	vector<vector<double>> &getData() {return data;}
	const vector<vector<double>> &getData() const {return data;}

	double& operator()(size_t i, size_t j) {
        return data[i][j];
    }

    const double& operator()(size_t i, size_t j) const {
        return data[i][j];
    }

	//Fills the matrix with random values I got from cppreference
	void randomMatrix() {
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<double> dis(0, 100);

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

		Matrix result(A.numRows(), B.numColumns());

		for (size_t i = 0; i < A.numRows(); i++) {
			for (size_t j = 0; j < B.numColumns(); j++) {
			  for (size_t k = 0; k < A.numColumns(); k++)
			    result(i,j) += A(i,k) * B(k,j);
			}
		}

		return result;
	}

};



class MatrixTests {
public:
	//Tests the Naive Algorithm
	static void NaiveAlgorithmTest(const size_t row_1, const size_t row_2,
							const size_t column_1, const size_t column_2) {

		Matrix matrix_1(row_1, column_1);
		Matrix matrix_2(row_2, column_2);

		matrix_1.randomMatrix();
		matrix_2.randomMatrix();

		matrix_1.printMatrix();
		matrix_2.printMatrix();

		Matrix result = MatrixAlgorithms::NaiveAlgorithm(matrix_1, matrix_2);

		result.printMatrix();
	}

};





int main() {

	MatrixTests::NaiveAlgorithmTest(50, 40, 40, 30);

	return 0;
}