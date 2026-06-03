#pragma once
#include "Vector3D.h"

class BoundryCondition
{
public:
	int type;
	double Betta;
	double (*Func)(Vector3D& x, double t);
	virtual ~BoundryCondition() = default;
	virtual double GetValue(Vector3D& x, double);
	virtual int GetType();
};

class Dirichlet : public BoundryCondition
{
public:
	Dirichlet(double (*func)(Vector3D& x, double t));
};

class Neumann : public BoundryCondition
{
public:
	Neumann(double (*func)(Vector3D& x, double t));
};

class Robin : public BoundryCondition
{
public:
	Robin(double (*func)(Vector3D& x, double t), double betta);
};

class TimeStart : public BoundryCondition
{
public:
	TimeStart(double (*func)(Vector3D& x, double t));
};
