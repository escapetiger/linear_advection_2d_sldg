#include "solver.h"

int main(int argc, char const *argv[])
{
    const char *profile_2D_file = "data/profile_2D.dat";
    const char *residual_file = "data/residual.dat";
    const char *node_e_file = "data/node_e.dat";
    const char *node_u_file = "data/node_u.dat";
    const char *edge_e_file = "data/edge_e.dat";
    const char *edge_u_file = "data/edge_u.dat";
    const char *elem_e_file = "data/elem_e.dat";
    const char *elem_u_file = "data/elem_u.dat";

    double e[3];
    Solver solver;

    std::cout << "===== Options setup: started ====" << std::endl;
    solver.ReadOptions(argc, argv);
    std::cout << "===== Options setup: finished ====" << std::endl;
    std::cout << '\n';

    std::cout << "===== Mesh setup: started ====" << std::endl;
    solver.ReadMesh();
    std::cout << "===== Mesh setup: finished ====" << std::endl;
    std::cout << '\n';

    std::cout << "===== DG setup: started ====" << std::endl;
    solver.ReadDG();
    std::cout << "===== DG setup: finished ====" << std::endl;
    std::cout << '\n';

    std::cout << "===== Simulate: started ====" << std::endl;
    solver.Simulate();
    std::cout << "===== Simulate: finished ====" << std::endl;
    std::cout << '\n';

    std::cout << "===== Postprocess: started ====" << std::endl;
    // residual
    solver.CalcResidual();
    solver.GetResidual(3, e);
    solver.OutputResidual(residual_file);
    std::cout << "residual:"
              << ' ' << std::scientific << std::setprecision(4) << e[0]
              << ' ' << std::scientific << std::setprecision(4) << e[1]
              << ' ' << std::scientific << std::setprecision(4) << e[2]
              << '\n';
    // 2D profile
    solver.OutputSolution2D(profile_2D_file);
    // elem infos
    solver.PrintNodeE(node_e_file);
    solver.PrintNodeU(node_u_file);
    solver.PrintEdgeE(edge_e_file);
    solver.PrintEdgeU(edge_u_file);
    solver.PrintElemE(elem_e_file);
    solver.PrintElemU(elem_u_file);
    std::cout << "===== Postprocess: finished ====" << std::endl;
    std::cout << '\n';

    std::cout << "==== Clean: started ====" << '\n';
    solver.Destroy();
    std::cout << "==== Clean: finished ====" << '\n';
    std::cout << '\n';
}