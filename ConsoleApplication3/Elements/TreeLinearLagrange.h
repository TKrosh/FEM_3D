#pragma once
#include <iostream>
#include <vector>
class TreeLinearLagrange
{
private:
	std::vector<int> CoordsIndexes;
	int MaterialIndex;

public:
	TreeLinearLagrange(std::vector<int> coordsIndexes, int material);

	friend std::ostream& operator<<(std::ostream& os, TreeLinearLagrange& elem)
	{
		os << "Elem Coord Indexes: ";
		for (int i = 0; i < elem.CoordsIndexes.size(); i++) os << elem.CoordsIndexes[i] << " ";
		os << "(Material: " << elem.MaterialIndex << ")";
		return os;
	}
};

