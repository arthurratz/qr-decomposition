// qr_gen.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <fstream>

#pragma warning(disable:26451)

template<typename _In = std::vector<std::vector<std::double_t>>>
void qr_gen_symm_matrix(_In & matrix, const std::size_t size, \
	const std::int32_t min_val, const std::int32_t max_val)
{
	matrix.resize(size);
	for (std::size_t row = 0; row < size; row++)
		matrix[row].resize(size);

	std::srand((unsigned)std::time(nullptr));

	for (std::size_t row = 0; row < size; row++)
	{
		for (std::size_t col = row; col < size; col++)
			matrix[row][col] = std::rand() % max_val + min_val;

		for (std::size_t row2 = row; row2 < size; row2++)
			matrix[row2][row] = matrix[row][row2];
	}
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

template<typename _In = std::vector<std::vector<std::double_t>>>
void qr_write_matrix_to_csv_file(_In matrix, const std::string filename)
{
	std::ofstream ofs(filename, std::ofstream::out);
	for (std::size_t row = 0; row < matrix.size(); row++)
	{
		for (std::size_t col = 0; col < matrix[row].size(); col++)
		{
			ofs << matrix[row][col];
			if (col != matrix[row].size() - 1)
				ofs << ",";
		}

		if (row != matrix.size() - 1)
			ofs << "\n";
	}

	ofs.close();
}

int main()
{
	std::cout << "QR-Symmetric Matrix Generator v.0.11a by Arthur V. Ratz @ CodeProject.Com\n";
	std::cout << "=========================================================================\n";

	std::int32_t v_min = 0, v_max = 0;
	std::string filename = "\0"; std::size_t size = 0L;
	std::cout << "\nMatrix Size (NxN) N = "; std::cin >> size;
	std::cout << "\nMinimum value v_min = "; std::cin >> v_min;
	std::cout << "\nMaximum value v_max = "; std::cin >> v_max;
	std::cout << "\nOutput (*.csv) file: "; std::cin >> filename; std::cout << "\n";

	std::vector<std::vector<std::double_t>> matrix;
	qr_gen_symm_matrix(matrix, size, v_min, v_max);

	qr_write_matrix_to_csv_file(matrix, filename);
}