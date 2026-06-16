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
void DenseStorage::multiplicationByVector_T(std::vector<double>* v)
{
	std::vector<double> tmp = *v;
	for (int i = 0; i < size; i++)
	{
		double sum = 0.0;
		for (int j = 0; j < size; j++)
		{
			sum += matrix[j][i] * tmp[j];
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
		{
			if (matrix[i][j] == 0)
				std::cout << matrix[i][j] << " ";
			else
				std::cout << "* ";// matrix[i][j] << " ";
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

void DenseStorage::AddLocalMatrix_Garmonic(TreeLinearLagrange& elem, double** LocalMatrix)
{
	for (int i = 0; i < elem.n; ++i)
	{
		int global_i = elem.CoordsIndexes[i];
		int s_ind_i = 2 * global_i;
		int c_ind_i = s_ind_i + 1;
		int s_free_i = FreeIndex[s_ind_i];
		int c_free_i = FreeIndex[c_ind_i];

		for (int j = 0; j < elem.n; ++j)
		{
			int global_j = elem.CoordsIndexes[j];
			int s_ind_j = 2 * global_j;
			int c_ind_j = s_ind_j + 1;
			int s_free_j = FreeIndex[s_ind_j];
			int c_free_j = FreeIndex[c_ind_j];

			// блок ss (p_ij)
			double p_ij = LocalMatrix[2 * i][2 * j];
			if (s_free_i != -1 && s_free_j != -1)
				matrix[s_free_i][s_free_j] += p_ij;
			else if (s_free_i != -1 && s_free_j == -1)
				rightPart[s_free_i] -= p_ij * DirichletValues[s_ind_j];
			//else if (s_free_i == -1 && s_free_j != -1)
			//	rightPart[s_free_j] -= p_ij * DirichletValues[s_ind_i];

			// блок sc (-c_ij)
			double mc_ij = LocalMatrix[2 * i][2 * j + 1];
			if (s_free_i != -1 && c_free_j != -1)
				matrix[s_free_i][c_free_j] += mc_ij;
			else if (s_free_i != -1 && c_free_j == -1)
				rightPart[s_free_i] -= mc_ij * DirichletValues[c_ind_j];
			//else if (s_free_i == -1 && c_free_j != -1)
			//	rightPart[c_free_j] -= mc_ij * DirichletValues[s_ind_i];

			// блок cs (+c_ij)
			double pc_ij = LocalMatrix[2 * i + 1][2 * j];
			if (c_free_i != -1 && s_free_j != -1)
				matrix[c_free_i][s_free_j] += pc_ij;
			else if (c_free_i != -1 && s_free_j == -1)
				rightPart[c_free_i] -= pc_ij * DirichletValues[s_ind_j];
			//else if (c_free_i == -1 && s_free_j != -1)
			//	rightPart[s_free_j] -= pc_ij * DirichletValues[c_ind_i];

			// блок cc (p_ij)
			double p_ij2 = LocalMatrix[2 * i + 1][2 * j + 1];
			if (c_free_i != -1 && c_free_j != -1)
				matrix[c_free_i][c_free_j] += p_ij2;
			else if (c_free_i != -1 && c_free_j == -1)
				rightPart[c_free_i] -= p_ij2 * DirichletValues[c_ind_j];
			//else if (c_free_i == -1 && c_free_j != -1)
			//	rightPart[c_free_j] -= p_ij2 * DirichletValues[c_ind_i];
		}
	}
}

void DenseStorage::AddToRightPart_Garmonic(TreeLinearLagrange& elem, std::vector<double> v)
{
	for (int i = 0; i < elem.n; ++i)
	{
		int global_i = elem.CoordsIndexes[i];
		int s_free_i = FreeIndex[2 * global_i];
		if (s_free_i != -1)
		{
			rightPart[s_free_i] += v[2 * i];
		}


		int c_free_i = FreeIndex[2 * global_i + 1];
		if (c_free_i != -1)
		{
			rightPart[c_free_i] += v[2 * i + 1];
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

void DenseStorage::LU_decompose()
{
	for (int i = 0; i < size; i++)
	{
		for (int j = i + 1; j < size; j++)
		{
			double sumL = 0;
			double sumU = 0;

			for (int k = 0; k < i; k++)
			{
				sumU += matrix[i][k] * matrix[k][j];
				sumL += matrix[j][k] * matrix[k][i];
			}
			matrix[i][j] = matrix[i][j] - sumU;
			matrix[j][i] = (matrix[j][i] - sumL) / matrix[i][i];
			matrix[j][j] -= matrix[i][j] * matrix[j][i];
		}
	}
}

void DenseStorage::Solve_L(std::vector<double>& y)
{
	for (int i = 0; i < size; i++)
	{
		double sum = 0.0;
		for (int j = 0; j < i; j++)
		{
			sum += matrix[i][j] * y[j];
		}
		y[i] = rightPart[i] - sum;
	}
}

void DenseStorage::Solve_U(std::vector<double>& x)
{
	for (int i = size - 1; i >= 0; i--)
	{
		double sum = 0.0;
		for (int j = i + 1; j < size; j++)
		{
			sum += matrix[i][j] * x[j];
		}
		x[i] = (x[i] - sum) / matrix[i][i];
	}
}
