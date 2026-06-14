#include "Core/Reader.h"
#include "Core/GiperbolicProblemGarmonic.h"
#include "Matrix/DenseStorage.h"
#include "Matrix/ProfileMatrix.h"
#include "Core/BoundryCondition.h"
#include "Elements/GiperbolicMaterialGarmonic.h"

//koeficents
double omega = 1;
double lambda(Vector3D& p)
{
    return 1.0;
}
double sigma(Vector3D& p)
{
    return 1.0;
}
double xi(Vector3D& p)
{
    return 1.0;
}

//u
double ug_c(Vector3D& p, double t)
{
    //double val = p.X + p.Y + p.Z;
    double val = cos(p.Z);
    return val;
}
double ug_s(Vector3D& p, double t)
{
    //double val = p.X + p.Y + p.Z;
    double val = cos(p.X);
    return val;
}

//div(grad)
double div_grad_s(Vector3D& p, double t)
{
    //double div_grad_val = 0.0;
    double div_grad_val = -cos(p.X);
    return div_grad_val;
}
double div_grad_c(Vector3D& p, double t)
{
    //double div_grad_val = 0.0;
    double div_grad_val = -cos(p.Z);
    return div_grad_val;
}

//right parts
double f_s(Vector3D& p, double t)
{
    return (-lambda(p) * div_grad_s(p, t)) 
        - omega * sigma(p) * ug_c(p, t)
        - omega * omega * xi(p) * ug_s(p, t);
}

double f_c(Vector3D& p, double t)
{
    return (-lambda(p) * div_grad_c(p, t))
        + omega * sigma(p) * ug_s(p, t)
        - omega * omega * xi(p) * ug_c(p, t);
}



int SolveGiperbolicGarmonicProblems()
{
    Reader reader = Reader();
    Mesh TaskMesh;

    if (false)
    {
        const char* splitFileName = "tmpcube4.txt";
        reader.BuildcubeMesh(splitFileName, 16);
        return 0;
    }

    //const char* FileName = "SimpleMesh.txt";
    const char* FileName = "cube4.txt";
    reader.ReadMeshFromFile(FileName, TaskMesh);

    std::vector<double> time = { 0, 0.5, 1, 1.5, 2 };
    GiperbolicProblem<DenseStorage> Solver = GiperbolicProblem<DenseStorage>(
        {
            GiperbolicMaterial(lambda, sigma, xi, f_s, f_c, omega)
        },
        TaskMesh,
        {
            Dirichlet(ug_s)
        },
        {
            Dirichlet(ug_c)
        },
        time
    );
    Solver.Solve();

    return 0.0;
}