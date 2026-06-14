#pragma once
#include "../Core/Vector3D.h"

class GiperbolicMaterial
{
public:

    GiperbolicMaterial(double (*lambda)(Vector3D&),
        double (*sigma)(Vector3D&),
        double (*xi)(Vector3D&),
        double (*f_s)(Vector3D&, double),
        double (*f_c)(Vector3D&, double), 
        double omega);

    double (*Lambda)(Vector3D& p),
        (*Sigma)(Vector3D& p),
        (*Xi)(Vector3D& p),
        (*F_s)(Vector3D& p, double t),
        (*F_c)(Vector3D& p, double t);

    double Omega;
};

