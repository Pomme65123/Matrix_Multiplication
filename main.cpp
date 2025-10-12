#include <chrono>
#include <iostream>
#include <random>
#include <vector>

using namespace std;

//Makes matrix of size [rows][columns]
vector<vector<double>> createMatrix(size_t rows, size_t cols) {
	return vector<vector<double>>(rows, vector<double>(cols, 0.0));
}

//Fills the matrix with random values I got from cppreference
void randomMatrix(vector<vector<double>> &matrix) {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis(0, 100);

	for (size_t i = 0; i < matrix.size(); i++) {
		for (size_t j = 0; j < matrix[i].size(); j++) {
		  matrix[i][j] = dis(gen);
		}
	}
}

//Prints a matrix
void printMatrix(const vector<vector<double>> &matrix) {
	cout << "Matrix: " << endl;
	for (size_t i = 0; i < matrix.size(); i++) {
		for (size_t j = 0; j < matrix[0].size(); j++) {
			cout << "Position: [" << i << "][" << j << "] = " << matrix[i][j] << endl;
		}
	}
}

//Uses the Naive Algorithm for matrix multiplication
vector<vector<double>> NaiveAlgorithm(	const vector<vector<double>> &matrix_1, 
										const vector<vector<double>> &matrix_2) {

	//https://en.wikipedia.org/wiki/Matrix_multiplication

	if (matrix_1.empty() || matrix_2.empty() || matrix_1[0].size() != matrix_2.size()) {
	cout << "Invalid matrix multiplication" << endl;
		return {};
	}

	vector<vector<double>> result = createMatrix(matrix_1.size(), matrix_2[0].size());

	for (size_t i = 0; i < matrix_1.size(); i++) {
		for (size_t j = 0; j < matrix_2[0].size(); j++) {
		  for (size_t k = 0; k < matrix_1[0].size(); k++)
		    result[i][j] += matrix_1[i][k] * matrix_2[k][j];
		}
	}

	return result;
}

//Tests the Naive Algorithm
void NaiveAlgorithmTest(const size_t row_1, const size_t row_2,
						const size_t column_1, const size_t column_2) {

	vector<vector<double>> matrix_1 = createMatrix(row_1, column_1);
	vector<vector<double>> matrix_2 = createMatrix(row_2, column_2);

	randomMatrix(matrix_1);
	randomMatrix(matrix_2);

	printMatrix(matrix_1);
	printMatrix(matrix_2);

	vector<vector<double>> result = NaiveAlgorithm(matrix_1, matrix_2);

	printMatrix(result);
}




int main() {

	NaiveAlgorithmTest(5, 4, 4, 3);

	return 0;
}