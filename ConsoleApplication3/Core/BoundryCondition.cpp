#include "BoundryCondition.h"


Dirichlet::Dirichlet(double(*func)(Vector3D& x, double t))
{
	Func = func;
	type = 1;
}

double BoundryCondition::GetValue(Vector3D& x, double t)
{
	return Func(x, t);
}

int BoundryCondition::GetType()
{
	return type;
}

Neumann::Neumann(double(*func)(Vector3D& x, double t))
{
	Func = func;
	type = 2;
}

Robin::Robin(double(*func)(Vector3D& x, double t), double betta)
{
	Betta = betta;
	Func = func;
	type = 3;
}

TimeStart::TimeStart(double(*func)(Vector3D& x, double t))
{
	Func = func;
	type = 4;
}