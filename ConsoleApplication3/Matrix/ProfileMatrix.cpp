#include "ProfileMatrix.h"

ProfileMatrix::ProfileMatrix(int n,
	std::vector<int> indexStart,
	std::map<int, double> dirichletValues,
	std::vector<int> freeIndex)
{

	DirichletValues = dirichletValues;
	FreeIndex = freeIndex;

	size = n;
	di.resize(n, 0.0);
	rightPart.resize(n, 0.0);
	ia.resize(n + 1, 0.0);

	for (int i = 1; i <= n; i++)
	{
		ia[i] = ia[i - 1] + (i - 1 - indexStart[i - 1]);
	}
	al.resize(ia[n], 0.0);
}

void ProfileMatrix::multiplicationByVector(std::vector<double>* v)
{
	std::vector<double> tmp = *v;
	for (int i = 0; i < size; i++)
	{
		double sum = di[i] * tmp[i];

		int start = ia[i];
		int end = ia[i + 1];
		int lineStart = i - (end - start);
		for (int pos = start; pos < end; pos++)
		{
			int j = lineStart + (pos - start);
			sum += al[pos] * tmp[j];
		}
		
		for (int j = i + 1; j < size; j++)
		{
			int js = ia[j];
			int je = ia[j + 1];
			int jLineStart = j - (je - js);
			if (i >= jLineStart)
			{
				sum += al[js + (i - jLineStart)] * tmp[j];
			}
		}
		(*v)[i] = sum;
	}
}

void ProfileMatrix::matrixClear()
{

	for (int i = 0; i < size; i++)
	{
		di[i] = 0;
		rightPart[i] = 0;
	}
	for (int i = 0; i < ia[size]; i++) al[i] = 0;
}

void ProfileMatrix::test()
{
	std::cout << "Global Matrix:" << std::endl;

	//first diag elem
	std::cout << di[0] << " " << std::endl;
	//for (int j = 1; j < size; j++) std::cout << 0 << " ";

	for (int i = 1; i < size; i++)
	{
		int CountOfElems_InLine = ia[i + 1] - ia[i];
		int startIndex = i - CountOfElems_InLine;
		int start = ia[i];
		for (int j = 0; j < i; j++)
		{
			if (j >= startIndex)
			{
				if (al[start + j - startIndex] == 0)
					std::cout << 0 << " ";
				else
					std::cout << al[start + j - startIndex] << " ";
			}
			else
			{
				std::cout << 0 << " ";
			}
		}
		std::cout << di[i] << " " << std::endl;
	}
	std::cout << std::endl;
	for (double i : al)
	{
		std::cout << i << " ";
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

void ProfileMatrix::AddLocalMatrix(TreeLinearLagrange& elem, double** LocalMatrix)
{
	for (int i = 0; i < elem.n; ++i)
	{
		int global_i = elem.CoordsIndexes[i];
		int free_i = FreeIndex[global_i];

		// ƒиагональный вклад
		if (free_i != -1) {
			di[free_i] += LocalMatrix[i][i];
		}

		// Ќижний треугольник (j < i)
		for (int j = 0; j < i; ++j)
		{
			int global_j = elem.CoordsIndexes[j];
			int free_j = FreeIndex[global_j];
			double val = LocalMatrix[i][j];   // симметрична€ матрица, берЄм один раз

			if (free_i != -1 && free_j != -1)   // оба свободны
			{
				if (free_i > free_j)
				{
					// добавл€ем в строку free_i, столбец free_j
					int CountOfElems_InLine = ia[free_i + 1] - ia[free_i];
					int lineStartInd = free_i - CountOfElems_InLine;
					if (free_j >= lineStartInd) {
						int offset = free_j - lineStartInd;
						al[ia[free_i] + offset] += val;
					}
				}
				else  // free_i < free_j
				{
					// добавл€ем симметрично в строку free_j, столбец free_i
					int CountOfElems_InLine = ia[free_j + 1] - ia[free_j];
					int lineStartInd = free_j - CountOfElems_InLine;
					if (free_i >= lineStartInd) {
						int offset = free_i - lineStartInd;
						al[ia[free_j] + offset] += val;
					}
				}
			}
			else if (free_i != -1 && free_j == -1)  // i свободен, j Ц Dirichlet
			{
				double u_j = DirichletValues[global_j];
				rightPart[free_i] -= val * u_j;
			}
			else if (free_i == -1 && free_j != -1)  // j свободен, i Ц Dirichlet
			{
				double u_i = DirichletValues[global_i];
				rightPart[free_j] -= val * u_i;
			}
			// если оба -1 Ц ничего не делаем
		}
	}
}

void ProfileMatrix::AddToRightPart(TreeLinearLagrange& elem, std::vector<double> v)
{
	for (int i = 0; i < elem.n; ++i)
	{
		int global_i = elem.CoordsIndexes[i];
		int free_i = FreeIndex[global_i];
		if (free_i != -1) {
			rightPart[free_i] += v[i];
		}
		// вклад в закреплЄнные узлы не нужен
	}
}

void ProfileMatrix::AddDirihletBoundary(BiLinearLagrange& elem, std::vector<double> v)
{

}
