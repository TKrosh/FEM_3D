#include "Mesh.h"
#pragma once
//#include "../Elements/TreeLinearLagrange.h"
class Assembler
{
public:
	Mesh TaskMesh;

	Assembler(Mesh& taskMesh) { TaskMesh = taskMesh; };
	Assembler() { };

	void CalculateStiffness(double (*f)(Vector3D&),
		TreeLinearLagrange& elem,
		double ** localStiffness);

	void CalculateMass(double (*f)(Vector3D&),
		TreeLinearLagrange& elem,
		double ** localMass);

	std::vector<double> CalculateLoad(double (*f)(Vector3D&, double),
		TreeLinearLagrange& elem);

	std::vector<double> CallculateDirihletVector(double (*f)(Vector3D&, double t),
		BiLinearLagrange& elem);

	void CalculateStiffness_Garmonic(double (*lambda)(Vector3D&),
		double (*xi)(Vector3D&),
		TreeLinearLagrange& elem,
		double** localStiffness,
		double omega);

	void CalculateMass_Garmonic(double (*f)(Vector3D&),
		TreeLinearLagrange& elem,
		double** localMass,
		double omega);
};

