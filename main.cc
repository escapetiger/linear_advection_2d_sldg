#include "solver.h"
#include "io.h"

void Debug(int argc, char const *argv[]);
void ErrorOrder();
void RunSolver(SolverOptions *opt, double *e);

int main(int argc, char const *argv[])
{
    ErrorOrder();
    // Debug(argc, argv);
}

void Debug(int argc, char const *argv[])
{
    const char *profile_2D_file = "data/profile_2D.dat";
    const char *residual_file = "data/residual.dat";
    const char *node_e_file = "data/node_e.dat";
    const char *node_u_file = "data/node_u.dat";
    const char *edge_e_file = "data/edge_e.dat";
    const char *edge_u_file = "data/edge_u.dat";
    const char *elem_e_file = "data/elem_e.dat";
    const char *elem_u_file = "data/elem_u.dat";
    const char *sub_buf_e_file = "data/sub_buf_e.dat";
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
    io.PrintSubElemBufE(sub_buf_e_file);
    io.PrintSubElemBufU(sub_buf_u_file);
    io.PrintQuadRule(quad_file);
    std::cout << "===== Postprocess: finished ====" << std::endl;
    std::cout << '\n';

    std::cout << "==== Clean: started ====" << '\n';
    solver.Destroy();
    std::cout << "==== Clean: finished ====" << '\n';
    std::cout << '\n';
}

void ErrorOrder()
{
    const char *error_order_file = "data/error_order.dat";
    const int N = 5;
    SolverOptions opt[N];
    double Nx[N], Ny[N], e[N][3], o[N][3];
    for (int i = 0; i < N; i++)
    {
        Nx[i] = 10 * pow(2, i);
        Ny[i] = Nx[i];

        opt[i].xmin = 0.;           ///< minimum value of x
        opt[i].xmax = 1.;           ///< maximum value of x
        opt[i].xmin = 0.;           ///< minimum value of y
        opt[i].ymax = 1.;           ///< maximum value of y
        opt[i].tmax = 1.0;          ///< final simulation time
        opt[i].CFL = 3.0;            ///< CFL number
        opt[i].ax = 1.;             ///< wave speed in x
        opt[i].ay = 1.;             ///< wave speed in y
        opt[i].Nx = 10 * pow(2, i); ///< size of x-partition
        opt[i].Ny = opt[i].Nx;      ///< size of y-partition
        opt[i].order = 1;           ///< approximation order of numerical scheme
    }

    for (int i = 0; i < N; i++)
    {
        RunSolver(&opt[i], e[i]);
    }
    for (int i = 0; i < N; i++)
    {
        if (i == 0)
        {
            o[i][0] = o[i][1] = o[i][2] = 0.;
        }
        else
        {
            o[i][0] = log2(e[i - 1][0] / e[i][0]);
            o[i][1] = log2(e[i - 1][1] / e[i][1]);
            o[i][2] = log2(e[i - 1][2] / e[i][2]);
        }
    }

    std::ofstream out;
    out.open(error_order_file);
    for (int i = 0; i < N; i++)
    {
        out << std::fixed << std::setprecision(0) << Nx[i] << '\t'
            << std::fixed << std::setprecision(0) << Ny[i] << '\t'
            << std::scientific << std::setprecision(4) << e[i][0] << '\t'
            << std::fixed << std::setprecision(2) << o[i][0] << '\t'
            << std::scientific << std::setprecision(4) << e[i][1] << '\t'
            << std::fixed << std::setprecision(2) << o[i][1] << '\t'
            << std::scientific << std::setprecision(4) << e[i][2] << '\t'
            << std::fixed << std::setprecision(2) << o[i][2] << '\n';
    }
    out.close();
}

void RunSolver(SolverOptions *opt, double *e)
{
    Solver solver;
    IO io(solver);

    std::cout << "===== Options setup: started ====" << std::endl;
    solver.ReadOptions(opt);
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
    std::cout << "residual:"
              << ' ' << std::scientific << std::setprecision(4) << e[0]
              << ' ' << std::scientific << std::setprecision(4) << e[1]
              << ' ' << std::scientific << std::setprecision(4) << e[2]
              << '\n';
    std::cout << "===== Postprocess: finished ====" << std::endl;
    std::cout << '\n';

    std::cout << "==== Clean: started ====" << '\n';
    solver.Destroy();
    std::cout << "==== Clean: finished ====" << '\n';
    std::cout << '\n';
}