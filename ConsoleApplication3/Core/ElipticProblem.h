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
#include <map>

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
	std::map<int, double> dirichletValues;
	std::vector<int> freeIndex;

	ElipticProblem(std::vector<ElipticMaterial> materials, 
		Mesh& taskMesh, 
		std::vector<BoundryCondition> boundryConditions)
	{
		Materials = materials;
		TaskMesh = taskMesh;
		assembler = Assembler(taskMesh);
		BoundryConditions = boundryConditions;

		q.resize(TaskMesh.CountOfVertexes, 0.0);
		freeIndex.resize(TaskMesh.CountOfVertexes, -2);
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
		
		//GlobalSLAU.test();
		//std::vector<double> BoundaryRightPart;
		//for (int i = 0; i < TaskMesh.CountOfBoundaryElements; i++)
		//{
		//	BiLinearLagrange BoundaryElem = TaskMesh.BoundaryElements[i];
		//	int type = BoundryConditions[BoundaryElem.BoundaryIndex].type;
		//	switch (type)
		//	{
		//	case 1:
		//		BoundaryRightPart = assembler.CallculateDirihletVector(
		//			BoundryConditions[BoundaryElem.BoundaryIndex].Func, BoundaryElem);
		//		GlobalSLAU.AddDirihletBoundary(BoundaryElem, BoundaryRightPart);
		//		break;
		//	}
		//}

		MatrixSolver<Storage> solver = MatrixSolver<Storage>(10000, 1e-15, GlobalSLAU);
		std::vector<double> Free_q = solver.CGM_D();

		for (int i = 0; i < TaskMesh.CountOfVertexes; i++) {
			if (freeIndex[i] != -1) {
				q[i] = Free_q[freeIndex[i]];
			}
			else {
				q[i] = dirichletValues[i];
			}
		}

		double L2 = 0.0;

		for (int i = 0; i < TaskMesh.CountOfVertexes; i++)
		{
			Vector3D p = TaskMesh.Vertexes[i];
			bool specialPoint = p.X == 0.5 && p.Y == 0.5 && p.Z == 0.5;
			double res_FEM = q[i];
			double res_real = BoundryConditions[0].Func(p, 0);
			if (specialPoint) {
				std::cout << p << ") " << res_FEM << " " << res_real << " " << res_FEM - res_real << std::endl;
			}
			
			L2 += (res_FEM - res_real) * (res_FEM - res_real) * TaskMesh.Elements[0].jacobian;
		}

		std::cout << "L2 = " << sqrt(L2);
	}

	void BuildMatrixPortret()
	{
		int nFree, nTotal = TaskMesh.CountOfVertexes;
		AccountingDirihletBoundaris(nFree, nTotal);
		std::vector<int> line_start_free(nFree, -1);
		for (int i = 0; i < TaskMesh.CountOfElements; i++) {
			TreeLinearLagrange& elem = TaskMesh.Elements[i];
			int minFree = nTotal; // большое число
			// Ищем минимальный свободный индекс среди узлов элемента
			for (int idx : elem.CoordsIndexes) {
				int freeIdx = freeIndex[idx];
				if (freeIdx != -1 && freeIdx < minFree) minFree = freeIdx;
			}
			// Записываем начало строки для всех свободных узлов элемента
			if (minFree < nFree) {
				for (int idx : elem.CoordsIndexes) {
					int freeIdx = freeIndex[idx];
					if (freeIdx != -1) {
						if (line_start_free[freeIdx] < 0 || minFree < line_start_free[freeIdx])
							line_start_free[freeIdx] = minFree;
					}
				}
			}
		}
		GlobalSLAU = GlobalMatrix<Storage>(nFree, line_start_free, dirichletValues, freeIndex);
	}

	void AccountingDirihletBoundaris(int& nFree, int& nTotal)
	{
		int countDirihletDOF = 0;

		for (int i = 0; i < TaskMesh.CountOfBoundaryElements; i++) {
			BiLinearLagrange& be = TaskMesh.BoundaryElements[i];
			int type = BoundryConditions[be.BoundaryIndex].type;
			if (type != 1) continue;

			auto& bcFunc = BoundryConditions[be.BoundaryIndex].Func;

			for (int localIdx = 0; localIdx < be.n; localIdx++) {
				int globalIdx = be.CoordsIndexes[localIdx];
				if (freeIndex[globalIdx] == -1)
					continue;

				Vector3D point = TaskMesh.Vertexes[globalIdx];
				double val = bcFunc(point, 0.0);
				freeIndex[globalIdx] = -1;
				countDirihletDOF++;
				dirichletValues[globalIdx] = val;
			}
		}

		nFree = nTotal - countDirihletDOF;

		int cnt = 0;
		for (int i = 0; i < nTotal; i++) {
			if (freeIndex[i] == -2) {
				freeIndex[i] = cnt;
				cnt++;
			}
		}
	}
};
