#pragma once
#include <vector>
#include "Vector3D.h"
#include "TreeLinearLagrange.h"
#include "BiLinearLagrange.h"
#include <iostream>

class Mesh
{
private:
	int CountOfVertexes = 0;
	std::vector<Vector3D> Vertexes;

	int CountOfElements = 0;
	std::vector<TreeLinearLagrange> Elements;

	int CountOfBoundaryElements = 0;
	std::vector<BiLinearLagrange> BoundaryElements;
public:
	Mesh() {};
	void ShowInfo();

	void SetNumberOfVertexes(int n);
	void AddVertex(double x, double y, double z);

	void SetNumberOfElements(int n);
	void AddElement(std::vector<int> indexes, int material);

	void SetNumberOfBoundaryElements(int n);
	void AddBoundaryElement(std::vector<int> indexes, int material, int BoundaryIndex);
};

