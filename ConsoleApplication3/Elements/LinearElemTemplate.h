#pragma once
#include <vector>

class LinearElemTemplate
{
public:
	int n, m, minCoord;
	std::vector<int> CoordsIndexes;
	int MaterialIndex;

	virtual double GetBasis(double x, double y, double z, int index) = 0;
};

