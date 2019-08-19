// qr_decomposition.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>

template<typename _In = std::vector<std::vector<std::double_t>>,
		 typename _Out = std::vector<std::vector<std::double_t>>>
_Out qr_transpose(_In matrix)
{
	std::vector<std::vector<std::double_t>> result(matrix.size());
	for (std::size_t row = 0; row < matrix.size(); row++)
	{
		result[row].resize(matrix[row].size());
		for (std::size_t col = 0; col < matrix[row].size(); col++)
			result[row][col] = matrix[col][row];
	}

	return result;
}

template<typename _InX = std::vector<std::vector<std::double_t>>,
	   	 typename _InY = std::vector<std::vector<std::double_t>>,
		 typename _InZ = std::vector<std::vector<std::double_t>>>
void qr_matmul(_InX matrix1, _InY matrix2, _InZ & result)
{
	result = std::vector<std::vector<std::double_t>>();
	result.resize(matrix1.size());

	for (std::size_t row = 0; row < matrix1.size(); row++)
		result[row].resize(matrix1[row].size());

	for (std::size_t row = 0; row < matrix1.size(); row++)
		for (std::size_t col = 0; col < matrix1[row].size(); col++)
			for (std::size_t k = 0; k < matrix1[row].size(); k++)
				result[row][col] += matrix1[row][k] * matrix2[k][col];
}

template<typename _InX = std::vector<std::vector<std::double_t>>,
		 typename _InY = std::vector<std::vector<std::double_t>>>
void qr_q_initialize(_InX & matrix, _InY & result)
{
	std::vector<std::vector<std::double_t>> \
		matrix_t = qr_transpose(matrix); 
	
	result.resize(matrix_t.size());
	for (std::size_t row = 0; row < matrix_t.size(); row++)
		result[row].resize(matrix_t[row].size());

	for (std::size_t col = 0; col < result[0].size(); col++)
		result[0][col] = matrix_t[0][col];
}

template<typename _InX = std::vector<std::double_t>,
	     typename _InY = std::vector<std::double_t>>
std::double_t qr_get_proj_uv(_InX & vector_u, _InY & vector_v)
{
	std::double_t scalar = .0f;
	for (std::size_t row = 0; row < vector_v.size(); row++)
		scalar += vector_u[row] * vector_v[row];

	return scalar;
}

template<typename _In = std::vector<std::double_t>>
std::double_t qr_get_norm_v(_In & vector)
{
	std::double_t scalar = .0f;
	for (std::size_t row = 0; row < vector.size(); row++)
		scalar += std::powl(vector[row], 2.0f);

	return std::sqrtl(scalar);
}

template<typename _InX = std::vector<std::vector<std::double_t>>,
	typename _InY = std::vector<std::vector<std::double_t>>>
void qr_get_q_matrix(_InX matrix, _InY & result)
{
	qr_q_initialize(matrix, result);

	for (std::size_t row = 1; row < result.size(); row++)
	{
		std::vector<std::double_t> proj(result[row].size());
		for (std::size_t row_p = 0; row_p < row; row_p++)
		{
			std::double_t proj_uu = \
				qr_get_proj_uv(result[row_p], result[row_p]);
			std::double_t proj_uv = \
				qr_get_proj_uv(matrix[row], result[row_p]);

			std::double_t c = proj_uv / proj_uu;
			for (std::size_t col = 0; col < result[row_p].size(); col++)
				proj[col] += c * result[row_p][col];
		}

		for (std::size_t col = 0; col < result[0].size(); col++)
			result[row][col] = matrix[row][col] - proj[col];
	}

	for (std::size_t row = 0; row < result.size(); row++)
	{
		std::double_t norm = qr_get_norm_v(result[row]);
		for (std::size_t col = 0; col < result[row].size(); col++)
			result[row][col] /= norm;
	}
}

template<typename _In = \
	std::vector<std::vector<std::double_t>>>
bool qr_is_right_triangular(_In matrix)
{
	for (std::size_t row = 0; row < matrix.size(); row++)
		for (std::size_t col = 0; col < matrix[row].size(); col++)
			if ((row > col) && (std::fabsl(matrix[row][col]) > 10e-6)) 
				return false;

	return true;
}

template<typename _In  = std::vector<std::vector<std::double_t>>,
		 typename _Out = std::vector<std::vector<std::double_t>>>
_Out qr_decompose_matrix(_In matrix)
{
	std::vector<std::vector<std::double_t>> Q, R;
	do {
		qr_get_q_matrix(matrix, Q); 
		qr_matmul(Q, matrix, R); 
		qr_matmul(R, qr_transpose(Q), matrix);
	} while (!qr_is_right_triangular(matrix));

	return matrix;
}

std::vector<std::vector<std::double_t>> \
	qr_read_matrix_from_csv_file(const std::string filename)
{
	// Declare a vector of transactions
	std::vector<std::vector<std::double_t>> matrix;
	// Create an input file stream
	std::ifstream ifs(filename, std::ifstream::in);
	// Allocate buffer for each line of file
	char* buf = (char*)malloc(sizeof(char) * 4096);
	// Perform a check if this is not the end of file
	while (ifs.eof() == false && ifs.peek() >= 0)
	{
		// Declare a tokens buffer
		std::vector<std::double_t> row;
		// Read the current line from file
		const char* delim = ","; ifs.getline(buf, 4096);
		// Extract the first token from the current string
		char* token = std::strtok(buf, delim);
		// Extract all tokens separated by ',' character
		while (token != nullptr)
		{
			// Add the current token (i.e. item) to the vector of tokens
			row.push_back(std::atof(token));
			// Extract the next token
			token = std::strtok(nullptr, delim);
		}

		// Add all items in the current vector to the vector of transactions
		matrix.push_back(row);
	}

	// Deallocate buffer and close the input file stream
	free(buf); ifs.close();

	// Return the vector of transactions
	return matrix;
}

template<typename _In = std::vector<std::vector<std::double_t>>>
void qr_print_matrix(_In matrix)
{
	for (std::size_t row = 0; row < matrix.size(); row++)
	{
		for (std::size_t col = 0; col < matrix[row].size(); col++)
			std::cout << matrix[row][col] << " ";

		std::cout << "\n";
	}
}

int main()
{
	std::cout << "QR-Decomposition v.0.11a by Arthur V. Ratz @ CodeProject.Com\n";
	std::cout << "============================================================\n";

	std::string filename = "\0";
	std::cout << "\nInput (*.csv) file: "; std::cin >> filename; std::cout << "\n";

	std::vector<std::vector<std::double_t>> A \
		= qr_read_matrix_from_csv_file(filename);

	qr_print_matrix(A); std::cout << "\n";

	std::vector<std::vector<std::double_t>> E \
		= qr_decompose_matrix(A);

	std::cout << "Eigenvalues:\n\n";

	std::vector<std::double_t> eigvals;
	for (std::size_t row = 0; row < E.size(); row++)
		eigvals.push_back(E[row][row]);

	std::sort(eigvals.begin(), eigvals.end(), 
		[&](const std::double_t& v1, std::double_t& v2) {
			return v1 > v2;
		});

	for (std::size_t row = 0; row < E.size(); row++)
		std::cout << "E(" << row + 1 << ") = " << eigvals[row] << "\n";

	std::cout << "\n";

	std::cin.get();
	std::cin.get();
}