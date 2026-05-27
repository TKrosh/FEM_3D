#include <iostream>
#include "Core/Reader.h"
#include "Core/Mesh.h"

int main()
{
    const char * FileName = "SimpleMesh.txt";
    Reader reader = Reader();

    Mesh TaskMesh;
    reader.ReadMeshFromFile(FileName, TaskMesh);

    TaskMesh.ShowInfo();
}