#include "GiperbolicMaterialGarmonic.h"

GiperbolicMaterial::GiperbolicMaterial(double (*lambda)(Vector3D&),
	double (*sigma)(Vector3D&),
	double (*xi)(Vector3D&),
	double (*f_s)(Vector3D&, double t),
	double (*f_c)(Vector3D&, double t), 
	double omega)
{
	Lambda = lambda;
	Sigma = sigma;
	Xi = xi;
	F_s = f_s;
	F_c = f_c;
	Omega = omega;
}