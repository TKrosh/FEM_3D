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
				std::cout << "*" << " ";
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

}

void DenseStorage::AddLocalMatrix(TreeLinearLagrange& elem, double** LocalMatrix)
{
	for (int i = 0; i < elem.n; i++)
	{
		for (int j = 0; j < elem.n; j++)
		{
			matrix[elem.CoordsIndexes[i]][elem.CoordsIndexes[j]] += LocalMatrix[i][j];
		}
	}
}

void DenseStorage::AddToRightPart(TreeLinearLagrange& elem, std::vector<double> v)
{
	for (int i = 0; i < elem.n; i++)
	{
		rightPart[elem.CoordsIndexes[i]] += v[i];
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

