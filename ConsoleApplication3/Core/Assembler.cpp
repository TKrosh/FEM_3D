#include "Assembler.h"

void Assembler::CalculateStiffness(double (*f)(Vector3D&),
	TreeLinearLagrange& elem,
	double** localStiffness)
{
	for (int i = 0; i < elem.n; i++)
	{
		for (int j = 0; j < elem.n; j++)
		{
			double integral = 0.0;
			for (int q = 0; q < elem.m; q++)
			{
				Vector3D QuadP = elem.quadraturePoints[q];
				Vector3D glopbalP = TaskMesh.PointToGlobal(QuadP, elem);
				double lambda_val = f(glopbalP);
				double dphi_i_dphi_j = elem.GetBasisDerivativeValue(QuadP.X, QuadP.Y, QuadP.Z, i, j);
				integral += lambda_val * dphi_i_dphi_j * elem.w[q] * elem.jacobian;
			}
			localStiffness[i][j] = integral;
		}
	}
}

void Assembler::CalculateMass(double (*f)(Vector3D&),
	TreeLinearLagrange& elem,
	double** localMass)
{
	double Basis[8];

	for (int i = 0; i < elem.n; i++)
		for (int j = 0; j < elem.n; j++)
			localMass[i][j] = 0.0;

	for (int q = 0; q < elem.m; q++)
	{
		Vector3D QuadP = elem.quadraturePoints[q];
		Vector3D GlobalP = TaskMesh.PointToGlobal(QuadP, elem);
		double FVal = f(GlobalP);
		double w_J = elem.w[q] * elem.jacobian;

		for (int i = 0; i < elem.n; i++)
			Basis[i] = elem.GetBasis(QuadP.X, QuadP.Y, QuadP.Z, i);
		for (int i = 0; i < elem.n; i++)
		{
			double fi_i = Basis[i];
			for (int j = 0; j < elem.n; j++)
			{
				localMass[i][j] += FVal * fi_i * Basis[j] * w_J;
			}
		}
	}
}

std::vector<double> Assembler::CalculateLoad(double (*f)(Vector3D&), 
	TreeLinearLagrange& elem)
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

std::vector<double> Assembler::CallculateDirihletVector(double(*f)(Vector3D&, double t), 
	BiLinearLagrange& elem)
{
	std::vector<double> res;


	for (int i : elem.CoordsIndexes)
	{
		res.push_back(f(TaskMesh.Vertexes[i], 0));
	}

	return res;
}
