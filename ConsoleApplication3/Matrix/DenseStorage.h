#pragma once
#include <vector>
#include "../Elements/TreeLinearLagrange.h"
#include "../Elements/BiLinearLagrange.h"
#include "../Elements/LinearElemTemplate.h"
#include <map>
class DenseStorage
{
public:
	int size, type = 0;
	std::vector<std::vector<double>> matrix;
	std::vector<double> rightPart;

	std::map<int, double> DirichletValues;
	std::vector<int> FreeIndex;

	DenseStorage() : size(0), matrix(size, std::vector<double>(size, 0.0)), rightPart(size, 0.0) {};
	DenseStorage(int n,
		std::vector<int> indexStart,
		std::map<int, double> dirichletValues,
		std::vector<int> freeIndex) : 
		size(n), matrix(n, std::vector<double>(n, 0.0)), 
		rightPart(size, 0.0),
		DirichletValues(dirichletValues),
		FreeIndex(freeIndex) {};

	void ReservedMemory_SpareStorage(std::vector<std::vector<int>> rowCols);

	void multiplicationByVector(std::vector<double>* v);
	void multiplicationByVector_T(std::vector<double>* v);
	void matrixClear();
	void test();

	void AddLocalMatrix(TreeLinearLagrange& elem, double** LocalMatrix);
	void AddToRightPart(TreeLinearLagrange& elem, std::vector<double> v);

	void AddLocalMatrix_Garmonic(TreeLinearLagrange& elem, double** LocalMatrix);
	void AddToRightPart_Garmonic(TreeLinearLagrange& elem, std::vector<double> v);

	void AddDirihletBoundary(BiLinearLagrange& elem, std::vector<double> v);
	void AddNeumannBoundary(int index, double theta); // in process
	void AddRobinBoundary(int index, double u_Betta, double betta); // in process

	void GetDiag(std::vector<double>* d);

	void LU_decompose();
	void Solve_L(std::vector<double>& y);
	void Solve_U(std::vector<double>& x);
};

