#include "ProfileMatrix.h"

ProfileMatrix::ProfileMatrix(int n, std::vector<int> indexStart)
{
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
					std::cout << "*" << " ";
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
}

void ProfileMatrix::AddLocalMatrix(TreeLinearLagrange& elem, double** LocalMatrix)
{
	for (int i = 0; i < elem.n; i++)
	{
		int globalLineStartind = ia[elem.CoordsIndexes[i]];
		int CountOfElems_InLine = (ia[elem.CoordsIndexes[i] + 1] - ia[elem.CoordsIndexes[i]]);
		int lineStartInd = elem.CoordsIndexes[i] - CountOfElems_InLine;

		di[elem.CoordsIndexes[i]] += LocalMatrix[i][i];
		for (int j = 0; j < i; j++)
		{
			if (globalLineStartind == 0)
			{
				int fdsaf = elem.CoordsIndexes[j] - lineStartInd;
			}
			if (globalLineStartind >= ia[elem.CoordsIndexes[j]])
			{
				int col_ind = elem.CoordsIndexes[j] - lineStartInd;
				al[globalLineStartind + col_ind] += LocalMatrix[i][j];
				//std::cout << globalLineStartind + col_ind << "[" << i << " " << j << "]"
				//	<< "[" << elem.CoordsIndexes[i] << " " << elem.CoordsIndexes[j] << "]" << std::endl;
			}
			else
			{
				int globalLineStartind = ia[elem.CoordsIndexes[j]];
				int CountOfElems_InLine = (ia[elem.CoordsIndexes[j] + 1] - ia[elem.CoordsIndexes[j]]);
				int lineStartInd = elem.CoordsIndexes[j] - CountOfElems_InLine;

				int col_ind = elem.CoordsIndexes[i] - lineStartInd;
				al[globalLineStartind + col_ind] += LocalMatrix[j][i];
				//std::cout << globalLineStartind + col_ind << "[" << j << " " << i << "]"
				//	<< "[" << elem.CoordsIndexes[j] << " " << elem.CoordsIndexes[i] << "]" << std::endl;
			}
		}
	}
	//std::cout << std::endl;
	//std::cout << std::endl;
}

void ProfileMatrix::AddToRightPart(TreeLinearLagrange& elem, std::vector<double> v)
{
	for (int i = 0; i < elem.n; i++)
	{
		rightPart[elem.CoordsIndexes[i]] += v[i];
	}
}

void ProfileMatrix::AddDirihletBoundary(BiLinearLagrange& elem, std::vector<double> v)
{
	double weight = 1;
	for (int i = 0; i < elem.n; i++)
	{
		int elemIndex = elem.CoordsIndexes[i];
		//for (int j = 0; j < size; j++)
		//{
		//	matrix[elemIndex][j] = 0.0;
		//}
		di[elemIndex] = weight;
		rightPart[elemIndex] = v[i] * weight;
	}
	
	
}
