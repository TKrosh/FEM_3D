#pragma once
#include <vector>
#include "../Elements/TreeLinearLagrange.h"
#include "../Elements/BiLinearLagrange.h"
#include "../Elements/LinearElemTemplate.h"
#include <map>
class DenseStorage
{
public:
	int size;
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

	void multiplicationByVector(std::vector<double>* v);
	void matrixClear();
	void test();

	void AddLocalMatrix(TreeLinearLagrange& elem, double** LocalMatrix);
	void AddToRightPart(TreeLinearLagrange& elem, std::vector<double> v);

	void AddDirihletBoundary(BiLinearLagrange& elem, std::vector<double> v);
	void AddNeumannBoundary(int index, double theta); // in process
	void AddRobinBoundary(int index, double u_Betta, double betta); // in process
};

