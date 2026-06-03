#pragma once
#include <vector>
#include "Mesh.h"
#include "fstream"
class Reader
{
public:

	Reader() {};

	void ReadMeshFromFile(const char* FileName, Mesh& TaskMesh);
	void BuildcubeMesh(const char* FileName, int n);
};

