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
		for (int j = 0; j <= i; j++)
		{
			if (matrix[i][j] == 0)
				std::cout  << matrix[i][j]<< " ";
			else
				std::cout << matrix[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << "Global right Part:" << std::endl;
	for (int i = 0; i < size; i++)
	{
		std::cout << " " << rightPart[i] << " ";
	}
	std::cout << std::endl;
}

void DenseStorage::AddLocalMatrix(TreeLinearLagrange& elem, double** LocalMatrix)
{
	for (int i = 0; i < elem.n; ++i)
	{
		int global_i = elem.CoordsIndexes[i];
		int free_i = FreeIndex[global_i];

		for (int j = 0; j < elem.n; ++j)
		{
			int global_j = elem.CoordsIndexes[j];
			int free_j = FreeIndex[global_j];
			double val = LocalMatrix[i][j];

			if (free_i != -1 && free_j != -1)
			{
				matrix[free_i][free_j] += val;
			}
			else if (free_i != -1 && free_j == -1)
			{
				double u_j = DirichletValues[global_j];
				rightPart[free_i] -= val * u_j;
			}
		}
	}
}

void DenseStorage::AddToRightPart(TreeLinearLagrange& elem, std::vector<double> v)
{
	for (int i = 0; i < elem.n; ++i)
	{
		int global_i = elem.CoordsIndexes[i];
		int free_i = FreeIndex[global_i];
		if (free_i != -1)
		{
			rightPart[free_i] += v[i];
		}
	}
}


void DenseStorage::AddDirihletBoundary(BiLinearLagrange& elem, std::vector<double> v)
{
	double weight = 1;
	for (int i = 0; i < elem.n; i++)
	{
		int elemIndex = elem.CoordsIndexes[i];
		for (int j = 0; j < size; j++)
		{
			matrix[elemIndex][j] = 0.0;
		}
		matrix[elemIndex][elemIndex] = weight;
		rightPart[elemIndex] = v[i] * weight;
	}
}

void DenseStorage::GetDiag(std::vector<double>* d)
{
	for (int i = 0; i < size; i++)
		(*d)[i] = matrix[i][i];
}

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

