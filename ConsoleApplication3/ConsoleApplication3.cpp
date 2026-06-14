#include "Core/Reader.h"
#include "Core/ElipticProblem.h"
#include "Matrix/DenseStorage.h"
#include "Core/BoundryCondition.h"
#include "SolveEllipticProblems.h"
#include "SolveGiperbolicGarmonicProblem.h"


int main()
{
    SolveGiperbolicGarmonicProblems();
    //SolveEllipticProblems();
    return 0;
}