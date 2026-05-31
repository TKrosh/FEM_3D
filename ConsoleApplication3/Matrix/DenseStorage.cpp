#include "DenseStorage.h"

void DenseStorage::multiplicationByVector(std::vector<double>* v)
{
	std::vector<double> tmp = *v;
	for (int i = 0; i < size; i++)
	{
		double sum = 0.0;
		for (int j = 0; j < size; j++)
		{
			sum += matrix[i][j] * tmp[j];
		}
		(*v)[i] = sum;
	}
}

void DenseStorage::matrixClear()
{
	for (int i = 0; i < size; i++)
	{
		std::fill(matrix[i].begin(), matrix[i].end(), 0);
	}
	std::fill(rightPart.begin(), rightPart.end(), 0);
}

void DenseStorage::test()
{
	std::cout << "Global Matrix:" << std::endl;
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
			std::cout << " " << matrix[i][j] << " ";
		std::cout << std::endl;
	}
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << "Global right Part:" << std::endl;
	for (int i = 0; i < size; i++)
	{
		std::cout << " " << rightPart[i] << " ";
	}

}

void DenseStorage::AddLocalMatrix(TreeLinearLagrange elem, double** LocalMatrix)
{
}

void DenseStorage::AddToRightPart(TreeLinearLagrange elem, std::vector<double> v)
{
}

//void DenseStorage::AddDirihletBoundary(int index, double ug)
//{
//	for (int i = 0; i < size; i++)
//	{
//		matrix[index][i] = 0.0;
//	}
//	matrix[index][index] = 1.0;
//	rightPart[index] = ug;
//}
//
//void DenseStorage::AddNeumannBoundary(int index, double theta)
//{
//	rightPart[index] += theta;
//}
//
//void DenseStorage::AddRobinBoundary(int index, double u_Betta, double betta)
//{
//	matrix[index][index] += betta;
//	rightPart[index] += betta * u_Betta;
//}

