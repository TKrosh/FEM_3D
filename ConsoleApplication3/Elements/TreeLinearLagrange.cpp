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

	n = 8;
	m = w.size(); // amount of weight;
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

double TreeLinearLagrange::GetBasisDerivativeValue(double x, double y, double z, int i, int j)
{
	Vector3D fi_i = Vector3D();
	Vector3D fi_j = Vector3D();

	switch (i)
	{
	case 0:
		fi_i.X = d_BasisF1_X(x) * BasisF1_Y(y) * BasisF1_Z(z) / Hx;
		fi_i.Y = BasisF1_X(x) * d_BasisF1_Y(y) * BasisF1_Z(z) / Hy;
		fi_i.Z = BasisF1_X(x) * BasisF1_Y(y) * d_BasisF1_Z(z) / Hz;
		break;
	case 1:
		fi_i.X = d_BasisF2_X(x) * BasisF1_Y(y) * BasisF1_Z(z) / Hx;
		fi_i.Y = BasisF2_X(x) * d_BasisF1_Y(y) * BasisF1_Z(z) / Hy;
		fi_i.Z = BasisF2_X(x) * BasisF1_Y(y) * d_BasisF1_Z(z) / Hz;
		break;
	case 2:
		fi_i.X = d_BasisF2_X(x) * BasisF2_Y(y) * BasisF1_Z(z) / Hx;
		fi_i.Y = BasisF2_X(x) * d_BasisF2_Y(y) * BasisF1_Z(z) / Hy;
		fi_i.Z = BasisF2_X(x) * BasisF2_Y(y) * d_BasisF1_Z(z) / Hz;
		break;
	case 3:
		fi_i.X = d_BasisF1_X(x) * BasisF2_Y(y) * BasisF1_Z(z) / Hx;
		fi_i.Y = BasisF1_X(x) * d_BasisF2_Y(y) * BasisF1_Z(z) / Hy;
		fi_i.Z = BasisF1_X(x) * BasisF2_Y(y) * d_BasisF1_Z(z) / Hz;
		break;
	case 4:
		fi_i.X = d_BasisF1_X(x) * BasisF1_Y(y) * BasisF2_Z(z) / Hx;
		fi_i.Y = BasisF1_X(x) * d_BasisF1_Y(y) * BasisF2_Z(z) / Hy;
		fi_i.Z = BasisF1_X(x) * BasisF1_Y(y) * d_BasisF2_Z(z) / Hz;
		break;
	case 5:
		fi_i.X = d_BasisF2_X(x) * BasisF1_Y(y) * BasisF2_Z(z) / Hx;
		fi_i.Y = BasisF2_X(x) * d_BasisF1_Y(y) * BasisF2_Z(z) / Hy;
		fi_i.Z = BasisF2_X(x) * BasisF1_Y(y) * d_BasisF2_Z(z) / Hz;
		break;
	case 6:
		fi_i.X = d_BasisF2_X(x) * BasisF2_Y(y) * BasisF2_Z(z) / Hx;
		fi_i.Y = BasisF2_X(x) * d_BasisF2_Y(y) * BasisF2_Z(z) / Hy;
		fi_i.Z = BasisF2_X(x) * BasisF2_Y(y) * d_BasisF2_Z(z) / Hz;
		break;
	case 7:
		fi_i.X = d_BasisF1_X(x) * BasisF2_Y(y) * BasisF2_Z(z) / Hx;
		fi_i.Y = BasisF1_X(x) * d_BasisF2_Y(y) * BasisF2_Z(z) / Hy;
		fi_i.Z = BasisF1_X(x) * BasisF2_Y(y) * d_BasisF2_Z(z) / Hz;
		break;
	}


	switch (j)
	{
	case 0:
		fi_j.X = d_BasisF1_X(x) * BasisF1_Y(y) * BasisF1_Z(z) / Hx;
		fi_j.Y = BasisF1_X(x) * d_BasisF1_Y(y) * BasisF1_Z(z) / Hy;
		fi_j.Z = BasisF1_X(x) * BasisF1_Y(y) * d_BasisF1_Z(z) / Hz;
		break;
	case 1:
		fi_j.X = d_BasisF2_X(x) * BasisF1_Y(y) * BasisF1_Z(z) / Hx;
		fi_j.Y = BasisF2_X(x) * d_BasisF1_Y(y) * BasisF1_Z(z) / Hy;
		fi_j.Z = BasisF2_X(x) * BasisF1_Y(y) * d_BasisF1_Z(z) / Hz;
		break;
	case 2:
		fi_j.X = d_BasisF2_X(x) * BasisF2_Y(y) * BasisF1_Z(z) / Hx;
		fi_j.Y = BasisF2_X(x) * d_BasisF2_Y(y) * BasisF1_Z(z) / Hy;
		fi_j.Z = BasisF2_X(x) * BasisF2_Y(y) * d_BasisF1_Z(z) / Hz;
		break;
	case 3:
		fi_j.X = d_BasisF1_X(x) * BasisF2_Y(y) * BasisF1_Z(z) / Hx;
		fi_j.Y = BasisF1_X(x) * d_BasisF2_Y(y) * BasisF1_Z(z) / Hy;
		fi_j.Z = BasisF1_X(x) * BasisF2_Y(y) * d_BasisF1_Z(z) / Hz;
		break;
	case 4:
		fi_j.X = d_BasisF1_X(x) * BasisF1_Y(y) * BasisF2_Z(z) / Hx;
		fi_j.Y = BasisF1_X(x) * d_BasisF1_Y(y) * BasisF2_Z(z) / Hy;
		fi_j.Z = BasisF1_X(x) * BasisF1_Y(y) * d_BasisF2_Z(z) / Hz;
		break;
	case 5:
		fi_j.X = d_BasisF2_X(x) * BasisF1_Y(y) * BasisF2_Z(z) / Hx;
		fi_j.Y = BasisF2_X(x) * d_BasisF1_Y(y) * BasisF2_Z(z) / Hy;
		fi_j.Z = BasisF2_X(x) * BasisF1_Y(y) * d_BasisF2_Z(z) / Hz;
		break;
	case 6:
		fi_j.X = d_BasisF2_X(x) * BasisF2_Y(y) * BasisF2_Z(z) / Hx;
		fi_j.Y = BasisF2_X(x) * d_BasisF2_Y(y) * BasisF2_Z(z) / Hy;
		fi_j.Z = BasisF2_X(x) * BasisF2_Y(y) * d_BasisF2_Z(z) / Hz;
		break;
	case 7:
		fi_j.X = d_BasisF1_X(x) * BasisF2_Y(y) * BasisF2_Z(z) / Hx;
		fi_j.Y = BasisF1_X(x) * d_BasisF2_Y(y) * BasisF2_Z(z) / Hy;
		fi_j.Z = BasisF1_X(x) * BasisF2_Y(y) * d_BasisF2_Z(z) / Hz;
		break;
	}

	return fi_i * fi_j;
}