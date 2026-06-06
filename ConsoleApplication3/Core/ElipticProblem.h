#pragma once
#include "../Matrix/GlobalMatrix.h"
#include "../Elements/ElipticMaterial.h"
#include <vector>
#include "Mesh.h"
#include "Assembler.h"
#include "../Elements/TreeLinearLagrange.h"
#include "../Elements/BiLinearLagrange.h"
#include "BoundryCondition.h"
#include "../Matrix/MatrixSolver.h"

template <typename Storage>
class ElipticProblem
{
public:
	std::vector<BoundryCondition> BoundryConditions;
	std::vector<ElipticMaterial> Materials;
	Mesh TaskMesh;
	Assembler assembler;
	GlobalMatrix<Storage> GlobalSLAU;
	std::vector<double> q;

	ElipticProblem(std::vector<ElipticMaterial> materials, Mesh& taskMesh, std::vector<BoundryCondition> boundryConditions)
	{
		Materials = materials;
		TaskMesh = taskMesh;
		assembler = Assembler(taskMesh);
		q.resize(TaskMesh.CountOfVertexes, 0.0);
		BoundryConditions = boundryConditions;
		BuildMatrixPortret();
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
		std::vector<double> RightPart;

		for (int i = 0; i < TaskMesh.CountOfElements; i++)
		{
			TreeLinearLagrange elem = TaskMesh.Elements[i];
			assembler.CalculateStiffness(Materials[elem.MaterialIndex].Lambda, elem, G);
			assembler.CalculateMass(Materials[elem.MaterialIndex].Gamma, elem, M);
			RightPart = assembler.CalculateLoad(Materials[elem.MaterialIndex].F, elem);

			GlobalSLAU.AddLocalMatrix(elem, G);
			GlobalSLAU.AddLocalMatrix(elem, M);
			GlobalSLAU.AddToRightPart(elem, RightPart);

			//GlobalSLAU.test(); // show matrix changes
		}
		for (int r = 0; r < n; r++) delete[] G[r];
		delete[] G;

		for (int r = 0; r < n; r++) delete[] M[r];
		delete[] M;
		
		GlobalSLAU.test();

		std::vector<double> BoundaryRightPart;
		for (int i = 0; i < TaskMesh.CountOfBoundaryElements; i++)
		{
			BiLinearLagrange BoundaryElem = TaskMesh.BoundaryElements[i];
			int type = BoundryConditions[BoundaryElem.BoundaryIndex].type;
			switch (type)
			{
			case 1:
				BoundaryRightPart = assembler.CallculateDirihletVector(
					BoundryConditions[BoundaryElem.BoundaryIndex].Func, BoundaryElem);
				GlobalSLAU.AddDirihletBoundary(BoundaryElem, BoundaryRightPart);
				break;
			}
		}

		MatrixSolver<Storage> solver = MatrixSolver<Storage>(10000, 1e-15, GlobalSLAU);
		q = solver.LOC();

		for (int i = 0; i < TaskMesh.CountOfVertexes; i++)
		{
			Vector3D p = TaskMesh.Vertexes[i];
			bool specialPoint = true; //p.X == 0.5 && p.Y == 0.5 && p.Z == 0.5;
			if (specialPoint) {
				double res_FEM = q[i];
				double res_real = BoundryConditions[0].Func(p, 0);
				std::cout << p << ") " << res_FEM << " " << res_real << " " << res_FEM - res_real << std::endl;
			}
		}
	}

	void BuildMatrixPortret()
	{
		std::vector<int> line_start;
		line_start.resize(TaskMesh.CountOfVertexes, -1);
		for (int i = 0; i < TaskMesh.CountOfElements; i++)
		{
			TreeLinearLagrange elem = TaskMesh.Elements[i];

			for (int i : elem.CoordsIndexes)
			{
				if (line_start[i] < 0 || elem.minCoord < line_start[i]) line_start[i] = elem.minCoord;
			}
		}

		//for (int i : line_start)
		//{
		//	std::cout << i << std::endl;
		//}

		GlobalSLAU = GlobalMatrix<Storage>(TaskMesh.CountOfVertexes, line_start);
	}
};
