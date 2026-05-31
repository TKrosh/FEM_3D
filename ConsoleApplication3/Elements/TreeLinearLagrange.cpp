#include "TreeLinearLagrange.h"

TreeLinearLagrange::TreeLinearLagrange(std::vector<int> coordsIndexes, int material, double hx, double hy, double hz)
{
	int CountOfPoints = coordsIndexes.size();
	CoordsIndexes.reserve(CountOfPoints);

	for (int i = 0; i < CountOfPoints; i++)
	{
		CoordsIndexes.emplace_back(coordsIndexes[i]);
	}

	Hx = hx; 
	Hy = hy; 
	Hz = hz; 
	jacobian = Hx * Hy * Hz;

	MaterialIndex = material;
}

double TreeLinearLagrange::GetBasis(double x, double y, double z, int index)
{
	double val;
	switch (index)
	{
	case 0:
		val = BasisF1_X(x) * BasisF1_Y(y) * BasisF1_Z(z);
		break;
	case 1:
		val = BasisF2_X(x) * BasisF1_Y(y) * BasisF1_Z(z);
		break;
	case 2:
		val = BasisF2_X(x) * BasisF2_Y(y) * BasisF1_Z(z);
		break;
	case 3:
		val = BasisF1_X(x) * BasisF2_Y(y) * BasisF1_Z(z);
		break;
	case 4:
		val = BasisF1_X(x) * BasisF1_Y(y) * BasisF2_Z(z);
		break;
	case 5:
		val = BasisF2_X(x) * BasisF1_Y(y) * BasisF2_Z(z);
		break;
	case 6:
		val = BasisF2_X(x) * BasisF2_Y(y) * BasisF2_Z(z);
		break;
	case 7:
		val = BasisF1_X(x) * BasisF2_Y(y) * BasisF2_Z(z);
		break;
	}
	return val;
}

double TreeLinearLagrange::GetBasisDerivativeValue(double x, double y, double z, int index)
{
	double val;
	switch (index)
	{
	case 0:
		val = d_BasisF1_X(x) * d_BasisF1_Y(y) * d_BasisF1_Z(z);
		break;
	case 1:
		val = d_BasisF2_X(x) * d_BasisF1_Y(y) * d_BasisF1_Z(z);
		break;
	case 2:
		val = d_BasisF2_X(x) * d_BasisF2_Y(y) * d_BasisF1_Z(z);
		break;
	case 3:
		val = d_BasisF1_X(x) * d_BasisF2_Y(y) * d_BasisF1_Z(z);
		break;
	case 4:
		val = d_BasisF1_X(x) * d_BasisF1_Y(y) * d_BasisF2_Z(z);
		break;
	case 5:
		val = d_BasisF2_X(x) * d_BasisF1_Y(y) * d_BasisF2_Z(z);
		break;
	case 6:
		val = d_BasisF2_X(x) * d_BasisF2_Y(y) * d_BasisF2_Z(z);
		break;
	case 7:
		val = d_BasisF1_X(x) * d_BasisF2_Y(y) * d_BasisF2_Z(z);
		break;
	}
	return val;
}