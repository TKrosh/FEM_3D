#include "TreeLinearLagrange.h"

TreeLinearLagrange::TreeLinearLagrange(std::vector<int> coordsIndexes, int material)
{
	int CountOfPoints = coordsIndexes.size();
	CoordsIndexes.reserve(CountOfPoints);

	for (int i = 0; i < CountOfPoints; i++)
	{
		CoordsIndexes.emplace_back(coordsIndexes[i]);
	}

	MaterialIndex = material;
}