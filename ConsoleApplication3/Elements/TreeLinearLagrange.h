#pragma once
#include "../Core/Vector3D.h"
#include "../Elements/LinearElemTemplate.h"
#include <iostream>
#include <vector>
class TreeLinearLagrange : public LinearElemTemplate
{
public:
	std::vector<double> w = { 0.125, 0.125, 0.125, 0.125,
							  0.125, 0.125, 0.125, 0.125 };

	double a = 1.0 / sqrt(3);
	std::vector<Vector3D> quadraturePoints = {
		Vector3D((1 - a) / 2, (1 - a) / 2, (1 - a) / 2),
		Vector3D((1 - a) / 2, (1 - a) / 2, (1 + a) / 2),
		Vector3D((1 - a) / 2, (1 + a) / 2, (1 - a) / 2),
		Vector3D((1 - a) / 2, (1 + a) / 2, (1 + a) / 2),
		Vector3D((1 + a) / 2, (1 - a) / 2, (1 - a) / 2),
		Vector3D((1 + a) / 2, (1 - a) / 2, (1 + a) / 2),
		Vector3D((1 + a) / 2, (1 + a) / 2, (1 - a) / 2),
		Vector3D((1 + a) / 2, (1 + a) / 2, (1 + a) / 2)
	};

	//int m = w.size(); // amount of weight;
	//int n = 8;
	//std::vector<int> CoordsIndexes;
	//int MaterialIndex;

	double Hx, Hy, Hz, jacobian;

	TreeLinearLagrange(std::vector<int> coordsIndexes, int material, double hx, double hy, double hz);

	friend std::ostream& operator<<(std::ostream& os, TreeLinearLagrange& elem)
	{
		os << "Elem Coord Indexes: ";
		for (int i = 0; i < elem.CoordsIndexes.size(); i++) os << elem.CoordsIndexes[i] << " ";
		os << std::endl << "(Material: " << elem.MaterialIndex << ") ";
		os << "hx = " << elem.Hx << " hy = " << elem.Hy << " hz = " << elem.Hz << std::endl;
		os << "min index  = " << elem.minCoord << std::endl;
		os << "--------------------";

		return os;
	}

	double GetBasis(double x, double y, double z, int index);

	// grad fi_i * grad fi_j
	double GetBasisDerivativeValue(double x, double y, double z, int i, int j);

	//basis
	double BasisF1_X(double x)
	{
		return 1 - x;
	};
	double BasisF2_X(double x)
	{
		return x;
	};

	double BasisF1_Y(double y)
	{
		return 1 - y;
	};
	double BasisF2_Y(double y)
	{
		return y;
	};

	double BasisF1_Z(double z)
	{
		return 1 - z;
	};
	double BasisF2_Z(double z)
	{
		return z;
	};

	//Basis derivative
	double d_BasisF1_X(double x)
	{
		return -1;
	};
	double d_BasisF2_X(double x)
	{
		return 1;
	};

	double d_BasisF1_Y(double y)
	{
		return -1;
	};
	double d_BasisF2_Y(double y)
	{
		return 1;
	};

	double d_BasisF1_Z(double z)
	{
		return -1;
	};
	double d_BasisF2_Z(double z)
	{
		return 1;
	};
};

