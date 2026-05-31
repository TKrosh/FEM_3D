#include "Assembler.h"

void Assembler::CalculateStiffness(double (*f)(Vector3D),
	TreeLinearLagrange elem,
	double**& localStiffness)
{
	for (int i = 0; i < elem.n; i++)
		for (int j = 0; j < elem.n; j++)
			localStiffness[i][j] = 0; 

	for (int i = 0; i < elem.n; i++)
	{
		for (int j = 0; j < elem.n; j++)
		{
			double integral = 0.0;
			for (int q = 0; q < elem.m; q++)
			{
				Vector3D QuadP = elem.quadraturePoints[q];
				double weight = elem.w[q];
				Vector3D glopbalP = TaskMesh.PointToGlobal(QuadP, elem);
				double lambda_val = f(glopbalP);
				double dphi_i = elem.GetBasisDerivativeValue(QuadP.X, QuadP.Y, QuadP.Z, i);
				double dphi_j = elem.GetBasisDerivativeValue(QuadP.X, QuadP.Y, QuadP.Z, j);
				integral += lambda_val * dphi_i * dphi_j * weight / elem.jacobian;
			}
			localStiffness[i][j] = integral;
		}
	}
}

std::vector<double> Assembler::CalculateLoad(double (*f)(Vector3D), TreeLinearLagrange elem)
{
	std::vector<double> resultVector;
	resultVector.resize(elem.n, 0);

	for (int i = 0; i < elem.n; i++)
	{
		double value = 0.0;
		for (int q = 0; q < elem.m; q++)
		{
			Vector3D QuadP = elem.quadraturePoints[q];
			double psi = elem.GetBasis(QuadP.X, QuadP.Y, QuadP.Z, i);
			Vector3D global_point = TaskMesh.PointToGlobal(QuadP, elem);
			value += f(global_point) * psi * elem.w[q] * elem.jacobian;
		}
		resultVector[i] = value;
	}
	return resultVector;
}