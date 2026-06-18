#pragma once
#include <vector>
#include "../Elements/TreeLinearLagrange.h"
#include "../Elements/BiLinearLagrange.h"
#include "../Elements/LinearElemTemplate.h"
#include <map>
#include <algorithm>

class RowColStorage
{
public:
	int size, type = 2;
	std::vector<double> rightPart;
	// diag
	std::vector<double> di;
	// indexes
	std::vector<int> ia;
	std::vector<int> ja;
	// low triangle
	std::vector<double> al;
	// upper triangle
	std::vector<double> au;
	std::map<int, double> DirichletValues;
	std::vector<int> FreeIndex;

	// for preconditioning
	std::vector<double> L;
	std::vector<double> D;
	std::vector<double> U;

	RowColStorage(int n,
		std::vector<int> indexStart,
		std::map<int, double> dirichletValues,
		std::vector<int> freeIndex);
	RowColStorage() : size(0), di(size, 0.0), ia(size, 0.0), al(size, 0.0), au(size, 0.0), rightPart(size, 0.0) {};

	void ReservedMemory_SpareStorage(std::vector<std::vector<int>> rowCols);

	void multiplicationByVector(std::vector<double>* v);
	void multiplicationByVector_T(std::vector<double>* v);
	void matrixClear();
	void test();

	void AddLocalMatrix(TreeLinearLagrange& elem, double** LocalMatrix);
	void AddToRightPart(TreeLinearLagrange& elem, std::vector<double> v);

	void AddLocalMatrix_Garmonic(TreeLinearLagrange& elem, double** LocalMatrix);
	void AddToRightPart_Garmonic(TreeLinearLagrange& elem, std::vector<double> v);

	void GetDiag(std::vector<double>* d);

	void LU_preconditioning();
	void Solve_L(std::vector<double>& y);
	void Solve_U(std::vector<double>& x);
	void Solve_L_Transope(std::vector<double>& y);
	void Solve_U_Transope(std::vector<double>& x);

	// we can't make that, because LU decompisition breaks profile
	// so it doesn't make sense.
	void LU_decompose() {};

};

