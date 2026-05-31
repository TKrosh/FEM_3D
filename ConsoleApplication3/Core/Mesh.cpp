#include "Mesh.h"

//Points, really just points in area
void Mesh::SetNumberOfVertexes(int n)
{
	CountOfVertexes = n;
	Vertexes.reserve(n);
}
void Mesh::AddVertex(double x, double y, double z)
{
	Vertexes.emplace_back(Vector3D(x, y, z));
}

//Elements that build by points that we inicialisez hieght
void Mesh::SetNumberOfElements(int n)
{
	CountOfElements = n;
	Elements.reserve(n);
}
void Mesh::AddElement(std::vector<int> indexes, int material, double hx, double hy, double hz)
{
	Elements.emplace_back(TreeLinearLagrange(indexes, material, hx, hy, hz));
}

//Boundary Elements
void Mesh::SetNumberOfBoundaryElements(int n)
{
	CountOfBoundaryElements = n;
	BoundaryElements.reserve(n);
}

void Mesh::AddBoundaryElement(std::vector<int> indexes, int material, int BoundaryIndex)
{
	BoundaryElements.emplace_back(BiLinearLagrange(indexes, material, BoundaryIndex));
}

void Mesh::ShowInfo()
{
	for (int i = 0; i < CountOfVertexes; i++)
	{
		std::cout << Vertexes[i] << std::endl;
	}

	for (int i = 0; i < CountOfElements; i++)
	{
		std::cout << Elements[i] << std::endl;
	}
	
	for (int i = 0; i < CountOfBoundaryElements; i++)
	{
		std::cout << BoundaryElements[i] << std::endl;
	}
}

Vector3D Mesh::PointToGlobal(Vector3D p, TreeLinearLagrange& elem)
{
	Vector3D result(0.0, 0.0, 0.0);
	for (int i = 0; i < elem.n; ++i)
	{
		double N = elem.GetBasis(p.X, p.Y, p.Z, i); 
		int globalNodeIndex = elem.CoordsIndexes[i];
		const Vector3D& node = Vertexes[globalNodeIndex];
		result.X += N * node.X;
		result.Y += N * node.Y;
		result.Z += N * node.Z;
	}
	return result;
}
