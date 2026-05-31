#pragma once
#include "../Elements/TreeLinearLagrange.h"

template<typename StoragePolicy>
class GlobalMatrix
{
public:
	StoragePolicy storage;

	GlobalMatrix() : storage() {};
	GlobalMatrix(int n) : storage(n) {};

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

	void AddLocalMatrix(TreeLinearLagrange elem, double** LocalMatrix)
	{
		storage.AddLocalMatrix(elem, LocalMatrix);
	};

	void AddToRightPart(TreeLinearLagrange elem, std::vector<double> v)
	{
		storage.AddToRightPart(elem, v);
	};

	//void AddDirihletBoundary(int index, double ug)
	//{
	//	storage.AddDirihletBoundary(index, ug);
	//};

	//void AddNeumannBoundary(int index, double theta)
	//{
	//	storage.AddNeumannBoundary(index, theta);
	//};

	//void AddRobinBoundary(int index, double u_Betta, double betta)
	//{
	//	storage.AddRobinBoundary(index, u_Betta, betta);
	//};
};

