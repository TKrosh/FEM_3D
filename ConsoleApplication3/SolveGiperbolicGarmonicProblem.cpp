#include "Core/Reader.h"
#include "Core/GiperbolicProblemGarmonic.h"
#include "Matrix/DenseStorage.h"
#include "Matrix/ProfileMatrix.h"
#include "Core/BoundryCondition.h"
#include "Elements/GiperbolicMaterialGarmonic.h"

//koeficents
double omega = 0.5;
double lambda(Vector3D& p)
{
    return 2.0;
}
double sigma(Vector3D& p)
{
    return 3.0;
}
double xi(Vector3D& p)
{
    return 1.0;
}

//u
double ug_s(Vector3D& p, double t)
{
    //double val = p.X * p.X + p.Y * p.Y + p.Z * p.Z + 3 * p.X - p.Y + 12 * p.Z - 11;
    double val = cos(p.X) + cos(p.Y) + cos(p.Z);
    return val;
}
double ug_c(Vector3D& p, double t)
{
    //double val = p.X * p.X + p.Y * p.Y + p.Z * p.Z + p.X + 2 * p.Y + p.Z + 5;
    double val = sin(p.X) + sin(p.Y) + sin(p.Z);
    return val;
}

//div(grad)
double div_grad_s(Vector3D& p, double t)
{
    //double div_grad_val = 6.0;
    double div_grad_val = -(cos(p.X) + cos(p.Y) + cos(p.Z));
    return div_grad_val;
}
double div_grad_c(Vector3D& p, double t)
{
    //double div_grad_val = 6.0;
    double div_grad_val = -(sin(p.X) + sin(p.Y) + sin(p.Z));
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
        const char* splitFileName = "cubeTest3.txt";
        reader.BuildcubeMesh(splitFileName, 3);
        // cube5 содержит 27 элемента на грани
        return 0;
    }

    //const char* FileName = "SimpleMesh.txt";
    const char* FileName = "cube32.txt";
    reader.ReadMeshFromFile(FileName, TaskMesh);

    //GiperbolicProblem<DenseStorage> Solver = GiperbolicProblem<DenseStorage>(
    //GiperbolicProblem<ProfileMatrix> Solver = GiperbolicProblem<ProfileMatrix>(
    GiperbolicProblem<ProfileMatrix> Solver = GiperbolicProblem<ProfileMatrix>(
        {
            GiperbolicMaterial(lambda, sigma, xi, f_s, f_c, omega)
        },
        TaskMesh,
        {
            Dirichlet(ug_s)
        },
        {
            Dirichlet(ug_c)
        }
    );
    Solver.Solve();

    return 0.0;
}