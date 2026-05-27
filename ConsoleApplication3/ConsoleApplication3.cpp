#include <iostream>
#include "Reader.h"
#include "Mesh.h"

int main()
{
    const char * FileName = "SimpleMesh.txt";
    Reader reader = Reader();

    Mesh TaskMesh;
    reader.ReadMeshFromFile(FileName, TaskMesh);

    TaskMesh.ShowInfo();
}