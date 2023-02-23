#include "solver.h"
#include "io.h"

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
    const char *poi_buf_e_file = "data/poi_buf_e.dat";
    const char *seg_buf_e_file = "data/seg_buf_e.dat";
    const char *sub_buf_e_file = "data/sub_buf_e.dat";
    const char *poi_buf_u_file = "data/poi_buf_u.dat";
    const char *seg_buf_u_file = "data/seg_buf_u.dat";
    const char *sub_buf_u_file = "data/sub_buf_u.dat";
    const char *quad_file = "data/quad.dat";

    double e[3];
    Solver solver;
    IO io(solver);

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
    io.PrintNodeE(node_e_file);
    io.PrintNodeU(node_u_file);
    io.PrintEdgeE(edge_e_file);
    io.PrintEdgeU(edge_u_file);
    io.PrintElemE(elem_e_file);
    io.PrintElemU(elem_u_file);
    io.PrintPOIBufE(poi_buf_e_file);
    io.PrintSegmentBufE(seg_buf_e_file);
    io.PrintSubElemBufE(sub_buf_e_file);
    io.PrintPOIBufU(poi_buf_u_file);
    io.PrintSegmentBufU(seg_buf_u_file);
    io.PrintSubElemBufU(sub_buf_u_file);
    io.PrintQuadRule(quad_file);
    std::cout << "===== Postprocess: finished ====" << std::endl;
    std::cout << '\n';

    std::cout << "==== Clean: started ====" << '\n';
    solver.Destroy();
    std::cout << "==== Clean: finished ====" << '\n';
    std::cout << '\n';
}