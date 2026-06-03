#pragma once
#include "../Core/Vector3D.h"
#include "../Elements/LinearElemTemplate.h"
#include <vector>
#include <iostream>
class BiLinearLagrange : public LinearElemTemplate
{
public:
	std::vector<double> w = { 0.25, 0.25, 0.25, 0.25 };

	double a = 1.0 / sqrt(3);
	std::vector<Vector3D> quadraturePoints = {
		Vector3D((1 - a) / 2, (1 - a) / 2, 0),
		Vector3D((1 - a) / 2, (1 + a) / 2, 0),
		Vector3D((1 + a) / 2, (1 - a) / 2, 0),
		Vector3D((1 + a) / 2, (1 + a) / 2, 0)
	};

	//int m = w.size(); // amount of weight;
	//int n = 4;
	//std::vector<int> CoordsIndexes;
	//int MaterialIndex;

	Vector3D e1, e2;
	double jacobian;


	//BoundaryIndex - index in massive of boundary condition
	int BoundaryIndex;

	BiLinearLagrange(std::vector<int> coordsIndexes, int material, int boundaryIndex, 
		Vector3D& _e1, Vector3D& _e2);

	friend std::ostream& operator<<(std::ostream& os, BiLinearLagrange& elem)
	{
		os << "Boundary elem Coord Indexes: ";
		for (int i = 0; i < elem.CoordsIndexes.size(); i++) os << elem.CoordsIndexes[i] << " ";
		os << "(" << "material :" << elem.MaterialIndex
			<< "; BoundaryIndex: " << elem.BoundaryIndex
			<< ")";
		return os;
	}

	double GetBasis(double x, double y, double z, int index);
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
};

