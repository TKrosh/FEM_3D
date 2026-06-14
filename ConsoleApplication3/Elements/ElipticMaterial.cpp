#include "ElipticMaterial.h"

ElipticMaterial::ElipticMaterial(double (*lambda)(Vector3D&),
	double (*gamma)(Vector3D&),
	double (*f)(Vector3D&, double))
{
	Lambda = lambda;
	Gamma = gamma;
	F = f;
}