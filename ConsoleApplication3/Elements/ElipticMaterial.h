#pragma once
#include "../Core/Vector3D.h"
class ElipticMaterial
{
public:
    ElipticMaterial(double (*lambda)(Vector3D&),
        double (*gamma)(Vector3D&),
        double (*f)(Vector3D&));

    double (*Lambda)(Vector3D& p),
        (*Gamma)(Vector3D& p),
        (*F)(Vector3D& p);
};

