#pragma once
#include "Mesh.h"
#include "fstream"
class Reader
{
public:

	Reader() {};

	void ReadMeshFromFile(const char* FileName, Mesh& TaskMesh);
};

