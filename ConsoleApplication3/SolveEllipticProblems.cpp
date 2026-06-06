#include "Core/Reader.h"
#include "Core/ElipticProblem.h"
#include "Matrix/DenseStorage.h"
#include "Matrix/ProfileMatrix.h"
#include "Core/BoundryCondition.h"

double ug(Vector3D& p, double t)
{
    return cos(p.X);
    //return p.X + p.Y + p.Z;
    //return p.X * p.X + p.Y * p.Y + p.Z * p.Z;
}

double div_grad(Vector3D& p)
{
    //return 0.0;
    //return 6.0;
    return -cos(p.X);
}

double lambda(Vector3D& p)
{
    return 1.0;
}

double gamma(Vector3D& p)
{
    return 1.0;
}

double f(Vector3D& p)
{
    return (-lambda(p) * div_grad(p)) + gamma(p) * ug(p, 0);
}

int SolveEllipticProblems()
{
    //const char * FileName = "SimpleMesh.txt";
    //const char * FileName = "twoCube.txt";

    Reader reader = Reader();
    Mesh TaskMesh;

    if (false)
    {
        const char* splitFileName = "tmpcube4.txt";
        reader.BuildcubeMesh(splitFileName, 16);
        return 0;
    }

    //const char* FileName = "SimpleMesh.txt";
    const char* FileName = "cube1.txt";
    reader.ReadMeshFromFile(FileName, TaskMesh);
    //TaskMesh.ShowInfo();

    ElipticProblem<DenseStorage> Solver = ElipticProblem<DenseStorage>(
        {
            ElipticMaterial(lambda, gamma, f)
        },
        TaskMesh,
        {
            Dirichlet(ug)
        }
    );
    Solver.Solve();
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "Enother Matrix";
    ElipticProblem<ProfileMatrix> Solver2 = ElipticProblem<ProfileMatrix>(
        {
            ElipticMaterial(lambda, gamma, f)
        },
        TaskMesh,
        {
            Dirichlet(ug)
        }
    );
    Solver2.Solve();
    
    return 0.0;
}