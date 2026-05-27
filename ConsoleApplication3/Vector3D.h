#pragma once
#include <iostream>
#include <vector>
class Vector3D
{
public:
	double X;
	double Y;
	double Z;

	Vector3D(double x, double y, double z);
	Vector3D();

	
	friend std::ostream& operator<<(std::ostream& os, Vector3D& v)
	{
		os << '{' << v.X << ", " << v.Y << ", " << v.Z << '}';
		return os;
	}
};