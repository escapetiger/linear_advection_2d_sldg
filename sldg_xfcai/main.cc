#include "globals.h"
#include "solver.h"

int main(int argc, char const *argv[])
{
    Solver solver;
    solver.Setup();
    solver.Allocate();
    solver.Init();
    solver.Debug();
    solver.Output("profile.dat");
    solver.Destroy();
    return 0;
}

