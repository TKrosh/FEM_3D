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

		/* we need some way to find : hx, hy, hz;
		 let's do that little bit si,ple
		if it doesn't work will see
		 
			7------6
           /      /|
		  /      / |
		 4------5  |
		 |  3---|--2
	   hx| / hy | /        z  y
		 |/ hx  |/         |/
		 0------1          *-> x
		 */
		Vector3D P0 = TaskMesh.Vertexes[indexes[0]];
		double hx = TaskMesh.Vertexes[indexes[1]].X - P0.X;
		double hy = TaskMesh.Vertexes[indexes[3]].Y - P0.Y;
		double hz = TaskMesh.Vertexes[indexes[4]].Z - P0.Z;

		TaskMesh.AddElement(indexes, material, hx, hy, hz);
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
