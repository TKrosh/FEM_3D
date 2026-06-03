#pragma once
#include <vector>
#include "Vector3D.h"
#include "../Elements/TreeLinearLagrange.h"
#include "../Elements/BiLinearLagrange.h"
#include "../Elements/LinearElemTemplate.h"
#include <iostream>


class Mesh
{
public:

	int CountOfVertexes = 0;
	std::vector<Vector3D> Vertexes;

	int CountOfElements = 0;
	std::vector<TreeLinearLagrange> Elements;

	int CountOfBoundaryElements = 0;
	std::vector<BiLinearLagrange> BoundaryElements;

	Mesh() {};
	void ShowInfo();

	void SetNumberOfVertexes(int n);
	void AddVertex(double x, double y, double z);

	void SetNumberOfElements(int n);
	void AddElement(std::vector<int> indexes, int material, double hx, double hy, double hz);

	void SetNumberOfBoundaryElements(int n);
	void AddBoundaryElement(std::vector<int> indexes, int material, int BoundaryIndex,
		Vector3D& _e1, Vector3D& _e2);

	Vector3D PointToGlobal(Vector3D p, LinearElemTemplate& elem);
};

