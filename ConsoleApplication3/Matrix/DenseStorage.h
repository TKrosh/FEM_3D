#pragma once
#include <vector>
#include "../Elements/TreeLinearLagrange.h"
class DenseStorage
{
public:
	int size;
	std::vector<std::vector<double>> matrix;
	std::vector<double> rightPart;

	DenseStorage() : size(0), matrix(size, std::vector<double>(size, 0.0)), rightPart(size, 0.0) {};
	DenseStorage(int n) : size(n), matrix(n, std::vector<double>(n, 0.0)), rightPart(size, 0.0) {};

	void multiplicationByVector(std::vector<double>* v);
	void matrixClear();
	void test();

	void AddLocalMatrix(TreeLinearLagrange& elem, double** LocalMatrix);
	void AddToRightPart(TreeLinearLagrange& elem, std::vector<double> v);

	//void AddDirihletBoundary(int index, double ug);
	//void AddNeumannBoundary(int index, double theta);
	//void AddRobinBoundary(int index, double u_Betta, double betta);
};

