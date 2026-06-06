#pragma once
#include <vector>
#include "../Elements/TreeLinearLagrange.h"
#include "../Elements/BiLinearLagrange.h"
#include "../Elements/LinearElemTemplate.h"

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

	// upper triangle, we would need this if the matrix were not symmetric.
	//std::vector<double> au;

	ProfileMatrix(int n, std::vector<int> indexStart);
	ProfileMatrix() : size(0), di(size, 0.0), ia(size, 0.0), al(size, 0.0), rightPart(size, 0.0) {};

	void multiplicationByVector(std::vector<double>* v);
	void matrixClear();
	void test();

	void AddLocalMatrix(TreeLinearLagrange& elem, double** LocalMatrix);
	void AddToRightPart(TreeLinearLagrange& elem, std::vector<double> v);

	void AddDirihletBoundary(BiLinearLagrange& elem, std::vector<double> v);
	//void AddNeumannBoundary(int index, double theta); // in process
	//void AddRobinBoundary(int index, double u_Betta, double betta); // in process
};

