#pragma once
#include "../Matrix/GlobalMatrix.h"
#include "../Elements/GiperbolicMaterialGarmonic.h"
#include <vector>
#include "Mesh.h"
#include "Assembler.h"
#include "../Elements/TreeLinearLagrange.h"
#include "../Elements/BiLinearLagrange.h"
#include "BoundryCondition.h"
#include "../Matrix/MatrixSolver.h"
#include <map>
#include <chrono>
#include <algorithm>

template <typename Storage>
class GiperbolicProblem
{
public:
	std::vector<BoundryCondition> BoundryConditions_s;
	std::vector<BoundryCondition> BoundryConditions_c;
	std::vector<GiperbolicMaterial> Materials;
	Mesh TaskMesh;
	Assembler assembler;
	GlobalMatrix<Storage> GlobalSLAU;
	std::vector<double> q;
	std::map<int, double> dirichletValues;
	std::vector<int> freeIndex;
	MatrixSolver<Storage> solver;

	GiperbolicProblem(std::vector<GiperbolicMaterial> materials,
		Mesh& taskMesh,
		std::vector<BoundryCondition> boundryConditions_s, 
		std::vector<BoundryCondition> boundryConditions_c)
	{
		Materials = materials;
		TaskMesh = taskMesh;
		assembler = Assembler(taskMesh);
		BoundryConditions_s = boundryConditions_s;
		BoundryConditions_c = boundryConditions_c;

		q.resize(TaskMesh.CountOfVertexes * 2, 0.0);
		freeIndex.resize(TaskMesh.CountOfVertexes * 2, -2);// i'm not really sure here

		int MatrixType = GlobalSLAU.storage.type;
		switch (MatrixType)
		{
		case 0:
		case 1:
			BuildMatrixPortret();
			break;
		case 2:
			BuildMatrixPortret_SpareMatrix();
		}

	}

	void Solve()
	{
		int n = 8; //size of local matrix
		double** G_s = new double* [n];
		double** M_s = new double* [n];
		double** M_c = new double* [n];
		double** G_c = new double* [n];

		// matrix for combining matrix
		double** A = new double* [n * 2];

		for (int i = 0; i < n; i++)
		{
			G_s[i] = new double[n];
			M_s [i] = new double[n];
			G_c[i] = new double[n];
			M_c [i] = new double[n];
			A[2 * i] = new double[n * 2];
			A[2 * i + 1] = new double[n * 2];
		}
		std::vector<double> RightPart_s;
		std::vector<double> RightPart_c;
		std::vector<double> RightPart(n * 2);

		for (int k = 0; k < TaskMesh.CountOfElements; k++)
		{
			TreeLinearLagrange elem = TaskMesh.Elements[k];

			//sin
			assembler.CalculateStiffness_Garmonic(Materials[elem.MaterialIndex].Lambda, 
				Materials[elem.MaterialIndex].Xi,
				elem, G_s, Materials[elem.MaterialIndex].Omega);

			assembler.CalculateMass_Garmonic(Materials[elem.MaterialIndex].Sigma, elem, M_s, Materials[elem.MaterialIndex].Omega);

			RightPart_s = assembler.CalculateLoad(Materials[elem.MaterialIndex].F_s, elem);
			//----------------------------------------------------------------------------------
			//cos
			assembler.CalculateStiffness_Garmonic(Materials[elem.MaterialIndex].Lambda, 
				Materials[elem.MaterialIndex].Xi,
				elem, G_c, Materials[elem.MaterialIndex].Omega);

			assembler.CalculateMass_Garmonic(Materials[elem.MaterialIndex].Sigma, elem, M_c, Materials[elem.MaterialIndex].Omega);
			
			RightPart_c = assembler.CalculateLoad(Materials[elem.MaterialIndex].F_c, elem);

			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < n; j++)
				{
					A[2 * i][2 * j] = G_s[i][j];
					A[2 * i][2 * j + 1] = -M_s[i][j];
					A[2 * i + 1][2 * j] = M_c[i][j];
					A[2 * i + 1][2 * j + 1] = G_c[i][j];
				}
			}

			for (int i = 0; i < n; i++)
			{
				RightPart[i * 2] = RightPart_s[i];
				RightPart[i * 2 + 1] = RightPart_c[i];
			}

			GlobalSLAU.AddLocalMatrix_Garmonic(elem, A);
			GlobalSLAU.AddToRightPart_Garmonic(elem, RightPart);
		}

		//GlobalSLAU.test(); // show matrix changes
		std::cout << std::endl << "Start Solving SLAU..." << std::endl;

		for (int r = 0; r < n; r++)
		{
			delete[] G_s[r];
			delete[] M_s[r];
			delete[] G_c[r];
			delete[] M_c[r];
			delete[] A[2 * r];
			delete[] A[2 * r + 1];
		}
		delete[] G_s;
		delete[] M_s;
		delete[] G_c;
		delete[] M_c;
		delete[] A;

		solver = MatrixSolver<Storage>(10000, 1e-15, GlobalSLAU);
		SolveSLAU(1);
	};

	void SolveSLAU(int SolverType)
	{
		std::vector<double> Free_q;
		auto start = std::chrono::high_resolution_clock::now();
		switch (SolverType)
		{
		case 0:
			Free_q = solver.LU();
			break;
		case 1:
			Free_q = solver.LOC_D();
			break;
		case 2:
			Free_q = solver.CGM_D();
			break;
		}

		auto stop = std::chrono::high_resolution_clock::now();
		auto duration = duration_cast<std::chrono::microseconds>(stop - start);
		std::cout << "duration: " << duration.count() << " microseconds" << std::endl;
		// сборка финального вектора
		for (int i = 0; i < TaskMesh.CountOfVertexes * 2; i++) {
			if (freeIndex[i] != -1) {
				q[i] = Free_q[freeIndex[i]];
			}
			else {
				q[i] = dirichletValues[i];
			}
		}

		// Вывод результатов
		double l2_s = 0.0, l2_c = 0.0;
		for (int i = 0; i < TaskMesh.CountOfVertexes; ++i) {
			Vector3D p = TaskMesh.Vertexes[i];
			// FEM решение (s и c компоненты)
			double fem_s = q[2 * i];
			double fem_c = q[2 * i + 1];
			// Аналитическое решение (из краевых условий)
			double real_s = BoundryConditions_s[0].Func(p, 0.0);
			double real_c = BoundryConditions_c[0].Func(p, 0.0);
			double diff_s = fem_s - real_s;
			double diff_c = fem_c - real_c;

			// Норма L2 по узлам (для справки)
			double jac = TaskMesh.Elements[0].jacobian; // якобиан одного элемента (равномерная сетка)
			l2_s += diff_s * diff_s * jac;
			l2_c += diff_c * diff_c * jac;

			// Форматный вывод
			//std::cout << "(" << p.X << ", " << p.Y << ", " << p.Z << ") | "
			//	<< fem_s << ", " << real_s << ", " << diff_s << " | "
			//	<< fem_c << ", " << real_c << ", " << diff_c << "\n";
		}
		//std::cout << "-----------------------------------------------------------------------------\n";
		std::cout << "L2 error (s): " << sqrt(l2_s) << "\n";
		std::cout << "L2 error (c): " << sqrt(l2_c) << "\n";
		std::cout << "error mid: " << (sqrt(l2_c) + sqrt(l2_s)) / 2.0 << "\n";
	}

	void BuildMatrixPortret()
	{
		int nTotalDOF, nFree;
		AccountingDirihletBoundaris(nFree, nTotalDOF);

		std::vector<int> line_start_free(nFree, -1);
		for (int iel = 0; iel < TaskMesh.CountOfElements; iel++) {
			TreeLinearLagrange& elem = TaskMesh.Elements[iel];

			// собираем все свободные DOF элемента
			std::vector<int> dofs;
			int minFree = nTotalDOF;
			for (int idx : elem.CoordsIndexes) {
				int s_free = freeIndex[2 * idx];
				int c_free = freeIndex[2 * idx + 1];
				if (s_free != -1) { dofs.push_back(s_free); if (s_free < minFree) minFree = s_free; }
				if (c_free != -1) { dofs.push_back(c_free); if (c_free < minFree) minFree = c_free; }
			}

			if (minFree < nFree) {
				for (int d : dofs) {
					if (line_start_free[d] < 0 || minFree < line_start_free[d])
						line_start_free[d] = minFree;
				}
			}
		}
		GlobalSLAU = GlobalMatrix<Storage>(nFree, line_start_free, dirichletValues, freeIndex);
	}

	void BuildMatrixPortret_SpareMatrix()
	{
		int nTotalDOF, nFree;
		AccountingDirihletBoundaris(nFree, nTotalDOF);

		// Вектор множеств для накопления столбцов каждой строки
		std::vector<std::vector<int>> rowCols(nFree);

		// Временный массив для отметки добавленных столбцов в текущей строке (для исключения дублей)
		std::vector<int> marker(nFree, -1); // инициализируем -1, будем записывать туда индекс текущей строки

		for (int iel = 0; iel < TaskMesh.CountOfElements; ++iel) {
			TreeLinearLagrange& elem = TaskMesh.Elements[iel];

			std::vector<int> dofs;
			for (int node : elem.CoordsIndexes) {
				int s_free = freeIndex[2 * node];
				int c_free = freeIndex[2 * node + 1];
				if (s_free != -1) dofs.push_back(s_free);
				if (c_free != -1) dofs.push_back(c_free);
			}

			for (int row : dofs) {
				for (int col : dofs) {
					if (marker[col] != row) {
						marker[col] = row;
						rowCols[row].push_back(col);
					}
				}
			}
		}

		for (auto& cols : rowCols)
			std::sort(cols.begin(), cols.end());

		std::vector<int> empty = {};
		GlobalSLAU = GlobalMatrix<Storage>(nFree, empty, dirichletValues, freeIndex);
		GlobalSLAU.ReservedMemory_SpareStorage(rowCols);
	}

	void AccountingDirihletBoundaris(int& nFree, int& nTotalDOF)
	{
		int countDirihletDOF = 0;
		nTotalDOF = TaskMesh.CountOfVertexes * 2;
		freeIndex.assign(nTotalDOF, -2);
		dirichletValues.clear();

		for (int ib = 0; ib < TaskMesh.CountOfBoundaryElements; ib++) {
			BiLinearLagrange& be = TaskMesh.BoundaryElements[ib];
			int bndIdx = be.BoundaryIndex;
			int type_s = BoundryConditions_s[bndIdx].type;
			int type_c = BoundryConditions_c[bndIdx].type;
			if (type_s != 1 && type_c != 1) continue;

			auto& func_s = BoundryConditions_s[bndIdx].Func;
			auto& func_c = BoundryConditions_c[bndIdx].Func;

			for (int loc = 0; loc < be.n; loc++) {
				int node = be.CoordsIndexes[loc];
				Vector3D p = TaskMesh.Vertexes[node];

				int s_node = 2 * node;
				int c_node = s_node + 1;

				if (type_s == 1 && freeIndex[s_node] != -1) {
					freeIndex[s_node] = -1;
					dirichletValues[s_node] = func_s(p, 0.0);
					countDirihletDOF++;
				}
				if (type_c == 1 && freeIndex[c_node] != -1) {
					freeIndex[c_node] = -1;
					dirichletValues[c_node] = func_c(p, 0.0);
					countDirihletDOF++;
				}
			}
		}

		nFree = nTotalDOF - countDirihletDOF;
		int cnt = 0;
		for (int i = 0; i < nTotalDOF; i++) {

			if (freeIndex[i] == -2) {
				freeIndex[i] = cnt;
				cnt++;
			}
		}
	};
	
};

