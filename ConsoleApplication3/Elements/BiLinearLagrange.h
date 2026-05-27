#pragma once
#include <vector>
#include <iostream>
class BiLinearLagrange
{
private:
	std::vector<int> CoordsIndexes;
	int MaterialIndex;
	//BoundaryIndex - index in massive of boundary condition
	int BoundaryIndex;

public:
	BiLinearLagrange(std::vector<int> coordsIndexes, int material, int boundaryIndex);

	friend std::ostream& operator<<(std::ostream& os, BiLinearLagrange& elem)
	{
		os << "Boundary elem Coord Indexes: ";
		for (int i = 0; i < elem.CoordsIndexes.size(); i++) os << elem.CoordsIndexes[i] << " ";
		os << "(" << "material :" << elem.MaterialIndex
			<< "; BoundaryIndex: " << elem.BoundaryIndex
			<< ")";
		return os;
	}
};

