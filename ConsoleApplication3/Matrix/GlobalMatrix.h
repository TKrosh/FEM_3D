#pragma once
#include "../Elements/TreeLinearLagrange.h"
#include "../Elements/BiLinearLagrange.h"
#include <map>

template<typename StoragePolicy>
class GlobalMatrix
{
public:
	StoragePolicy storage;

	GlobalMatrix() : storage() {};
	GlobalMatrix(int n, 
		std::vector<int> indexStart, 
		std::map<int, double> dirichletValues,
		std::vector<int> freeIndex) :
		storage(n, indexStart, dirichletValues, freeIndex) {};

	void multiplicationByVector(std::vector<double>* v)
	{
		storage.multiplicationByVector(v);
	};

	void matrixClear()
	{
		storage.matrixClear();
	};

	void test()
	{
		storage.test();
	};

	void AddLocalMatrix(TreeLinearLagrange& elem, double** LocalMatrix)
	{
		storage.AddLocalMatrix(elem, LocalMatrix);
	};

	void AddToRightPart(TreeLinearLagrange& elem, std::vector<double> v)
	{
		storage.AddToRightPart(elem, v);
	};

	void AddDirihletBoundary(BiLinearLagrange& elem, std::vector<double> v)
	{
		storage.AddDirihletBoundary(elem, v);
	};

	void AddNeumannBoundary(int index, double theta)
	{
		storage.AddNeumannBoundary(index, theta);
	};

	void AddRobinBoundary(int index, double u_Betta, double betta)
	{
		storage.AddRobinBoundary(index, u_Betta, betta);
	};

	void GetDiag(std::vector<double>* d)
	{
		storage.GetDiag(d);
	};
};

