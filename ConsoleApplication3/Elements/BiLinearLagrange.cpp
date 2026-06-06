#include "BiLinearLagrange.h"

BiLinearLagrange::BiLinearLagrange(std::vector<int> coordsIndexes, int material, int boundaryIndex,
	Vector3D& _e1, Vector3D& _e2)
{
	int CountOfPoints = coordsIndexes.size();
	CoordsIndexes.reserve(CountOfPoints);
	minCoord = coordsIndexes[0];

	for (int i = 0; i < CountOfPoints; i++)
	{
		CoordsIndexes.emplace_back(coordsIndexes[i]);
		if (minCoord > coordsIndexes[i]) minCoord = coordsIndexes[i];
	}

	e1 = _e1;
	e2 = _e2;

	jacobian = sqrt(e1 * e1) * sqrt(e2 * e2);
	m = w.size(); // amount of weight;
	n = 4;

	MaterialIndex = material;
	BoundaryIndex = boundaryIndex;
}

double BiLinearLagrange::GetBasis(double x, double y, double z, int index)
{
	double val;
	switch (index)
	{
	case 0:
		val = BasisF1_X(x) * BasisF1_Y(y);
		break;
	case 1:
		val = BasisF2_X(x) * BasisF1_Y(y);
		break;
	case 2:
		val = BasisF2_X(x) * BasisF2_Y(y);
		break;
	case 3:
		val = BasisF1_X(x) * BasisF2_Y(y);
		break;
	}
	return val;
}

double BiLinearLagrange::GetBasisDerivativeValue(double x, double y, double z, int i, int j)
{
	Vector3D fi_i = Vector3D();
	Vector3D fi_j = Vector3D();

	double h_ksi = sqrt(e1 * e1);
	double h_eta = sqrt(e1 * e1);

	switch (i)
	{
	case 0:
		fi_i.X = d_BasisF1_X(x) * BasisF1_Y(y) / h_ksi;
		fi_i.Y = BasisF1_X(x) * d_BasisF1_Y(y) / h_eta;
		break;
	case 1:
		fi_i.X = d_BasisF2_X(x) * BasisF1_Y(y) / h_ksi;
		fi_i.Y = BasisF2_X(x) * d_BasisF1_Y(y) / h_eta;
		break;
	case 2:
		fi_i.X = d_BasisF2_X(x) * BasisF2_Y(y) / h_ksi;
		fi_i.Y = BasisF2_X(x) * d_BasisF2_Y(y) / h_eta;
		break;
	case 3:
		fi_i.X = d_BasisF1_X(x) * BasisF2_Y(y) / h_ksi;
		fi_i.Y = BasisF1_X(x) * d_BasisF2_Y(y) / h_eta;
		break;
	}


	switch (j)
	{
	case 0:
		fi_j.X = d_BasisF1_X(x) * BasisF1_Y(y) / h_ksi;
		fi_j.Y = BasisF1_X(x) * d_BasisF1_Y(y) / h_eta;
		break;
	case 1:
		fi_j.X = d_BasisF2_X(x) * BasisF1_Y(y) / h_ksi;
		fi_j.Y = BasisF2_X(x) * d_BasisF1_Y(y) / h_eta;
		break;
	case 2:
		fi_j.X = d_BasisF2_X(x) * BasisF2_Y(y) / h_ksi;
		fi_j.Y = BasisF2_X(x) * d_BasisF2_Y(y) / h_eta;
		break;
	case 3:
		fi_j.X = d_BasisF1_X(x) * BasisF2_Y(y) / h_ksi;
		fi_j.Y = BasisF1_X(x) * d_BasisF2_Y(y) / h_eta;
		break;
	}


	fi_j.Z = 0.0;
	fi_i.Z = 0.0;

	return fi_i * fi_j;
}