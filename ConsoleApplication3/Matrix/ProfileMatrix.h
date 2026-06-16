#pragma once
#include <vector>
#include "../Elements/TreeLinearLagrange.h"
#include "../Elements/BiLinearLagrange.h"
#include "../Elements/LinearElemTemplate.h"
#include <map>

class ProfileMatrix
{
public:
	int size;
	std::vector<double> rightPart;
	// diag
	std::vector<double> di;
	// indexes
	std::vector<int> ia;
	// low triangle
	std::vector<double> al;
	// upper triangle
	std::vector<double> au;
	std::map<int, double> DirichletValues;
	std::vector<int> FreeIndex;



	ProfileMatrix(int n,
		std::vector<int> indexStart,
		std::map<int, double> dirichletValues,
		std::vector<int> freeIndex);
	ProfileMatrix() : size(0), di(size, 0.0), ia(size, 0.0), al(size, 0.0), au(size, 0.0), rightPart(size, 0.0) {};

	void multiplicationByVector(std::vector<double>* v);
	void multiplicationByVector_T(std::vector<double>* v);
	void matrixClear();
	void test();

	void AddLocalMatrix(TreeLinearLagrange& elem, double** LocalMatrix);
	void AddToRightPart(TreeLinearLagrange& elem, std::vector<double> v);

	void AddLocalMatrix_Garmonic(TreeLinearLagrange& elem, double** LocalMatrix);
	void AddToRightPart_Garmonic(TreeLinearLagrange& elem, std::vector<double> v);

	void GetDiag(std::vector<double>* d);
	//void AddNeumannBoundary(int index, double theta); // in process
	//void AddRobinBoundary(int index, double u_Betta, double betta); // in process

	void LU_decompose();
	void Solve_L(std::vector<double>& y);
	void Solve_U(std::vector<double>& x);
};

