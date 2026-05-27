#include "Reader.h"

void Reader::ReadMeshFromFile(const char* FileName, Mesh& TaskMesh)
{
	std::ifstream InfoFile;
	InfoFile.open(FileName);

	int n, material, Bindex;
	// time to start read Points
	InfoFile >> n;
	TaskMesh.SetNumberOfVertexes(n);
	for (int i = 0; i < n; i++)
	{
		double x, y, z;
		InfoFile >> x >> y >> z;
		TaskMesh.AddVertex(x, y, z);
	}

	// time to start read elements
	std::vector<int> indexes{0, 0, 0, 0, 0, 0, 0, 0};
	InfoFile >> n;
	TaskMesh.SetNumberOfElements(n);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < 8; j++) InfoFile >> indexes[j];
		InfoFile >> material;
		TaskMesh.AddElement(indexes, material);
	}

	// time to start read boundary elements
	std::vector<int> Boundaryindexes{ 0, 0, 0, 0};
	InfoFile >> n;
	TaskMesh.SetNumberOfBoundaryElements(n);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < 4; j++) InfoFile >> Boundaryindexes[j];
		InfoFile >> material >> Bindex;
		TaskMesh.AddBoundaryElement(Boundaryindexes, material, Bindex);
	}
}
