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

	friend double operator* (const Vector3D& v1, const Vector3D& v2)
	{
		return v1.X * v2.X + v1.Y * v2.Y + v1.Z * v2.Z;
	}

	friend Vector3D operator- (const Vector3D& v1, const Vector3D& v2)
	{
		Vector3D res = Vector3D();
		res.X = v1.X - v2.X;
		res.Y = v1.Y - v2.Y;
		res.Z = v1.Z - v2.Z;
		return res;
	}

	friend std::ostream& operator<<(std::ostream& os, Vector3D& v)
	{
		os << '{' << v.X << ", " << v.Y << ", " << v.Z << '}';
		return os;
	}
};