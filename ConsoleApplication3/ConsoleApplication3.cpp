#include "Core/Reader.h"
#include "Core/ElipticProblem.h"
#include "Matrix/DenseStorage.h"

double ug(Vector3D& p)
{
    return p.X + p.Y + p.Z;
}

double div_grad(Vector3D& p)
{
    return 0.0;
}

double lambda(Vector3D& p)
{
    return 1.0;
}

double gamma(Vector3D& p)
{
    return 0.0;
}

double f(Vector3D& p)
{
    return (-lambda(p) * div_grad(p)) + gamma(p) * ug(p);
}

int main()
{
    //const char * FileName = "SimpleMesh.txt";
    const char * FileName = "twoCube.txt";
    Reader reader = Reader();
    Mesh TaskMesh;

    reader.ReadMeshFromFile(FileName, TaskMesh);
    //TaskMesh.ShowInfo();

    ElipticProblem<DenseStorage> Solver = ElipticProblem<DenseStorage>(
        {
            ElipticMaterial(lambda, gamma, f)
        },
        TaskMesh
    );

    Solver.Solve();
}