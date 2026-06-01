#pragma once

#include "../Matrix/GlobalMatrix.h"
#include "../Elements/ElipticMaterial.h"
#include <vector>
#include "Mesh.h"
#include "Assembler.h"
#include "../Elements/TreeLinearLagrange.h"

template <typename Storage>
class ElipticProblem
{
public:
	std::vector<ElipticMaterial> Materials;
	Mesh TaskMesh;
	Assembler assembler;
	GlobalMatrix<Storage> GlobalSLAU;
	std::vector<double> q;

	ElipticProblem(std::vector<ElipticMaterial> materials, Mesh& taskMesh)
	{
		Materials = materials;
		TaskMesh = taskMesh;
		assembler = Assembler(taskMesh);
		GlobalSLAU = GlobalMatrix<Storage>(TaskMesh.CountOfVertexes);
		q.resize(GlobalSLAU.storage.size, 0.0);
	};

	void Solve()
	{
		int n = 8; //size of local matrix
		double** G = new double* [n];
		for (int i = 0; i < n; i++)
		{
			G[i] = new double[n];
		}

		double** M = new double* [n];
		for (int i = 0; i < n; i++)
		{
			M[i] = new double[n];
		}

		for (int i = 0; i < TaskMesh.CountOfElements; i++)
		{
			TreeLinearLagrange elem = TaskMesh.Elements[i];
			assembler.CalculateStiffness(Materials[elem.MaterialIndex].Lambda, elem, G);
			assembler.CalculateMass(Materials[elem.MaterialIndex].Gamma, elem, M);
			std::vector<double> RightPart = assembler.CalculateLoad(Materials[elem.MaterialIndex].F, elem);

			GlobalSLAU.AddLocalMatrix(elem, G);
			GlobalSLAU.AddLocalMatrix(elem, M);
			GlobalSLAU.AddToRightPart(elem, RightPart);

			GlobalSLAU.test();
		}

		for (int i = 0; i < TaskMesh.CountOfBoundaryElements; i++)
		{

		}
		
		for (int r = 0; r < 2; r++) delete[] G[r];
		delete[] G;

		for (int r = 0; r < 2; r++) delete[] M[r];
		delete[] M;
	}
};
