#include "solver.h"

using mcm::DeleteMatrix;
using mcm::DeleteVector;
using mcm::Index2D;
using mcm::NewMatrix;
using mcm::NewVector;
using mcm::kernels::FuzzyZero;

void Solver::ReadOptions(int argc, char const *argv[])
{
    CommandLineReader reader(argc, argv);
    reader.AddOption(&xmin, "-xa", "xmin", "Minimum value of x", false);
    reader.AddOption(&xmax, "-xb", "xmax", "Maximum value of x", false);
    reader.AddOption(&ymin, "-ya", "ymin", "Minimum value of y", false);
    reader.AddOption(&ymax, "-yb", "ymax", "Maximum value of y", false);
    reader.AddOption(&tmax, "-T", "tmax", "Simulation time", true);
    reader.AddOption(&ax, "-ax", "ax", "wave speed in x", false);
    reader.AddOption(&ay, "-ay", "ay", "wave speed in y", false);
    reader.AddOption(&Nx, "-Nx", "Nx", "Cardinal Number of x-partitioin", true);
    reader.AddOption(&Ny, "-Ny", "Ny", "Cardinal Number of y-partitioin", true);
    reader.AddOption(&CFL, "-CFL", "CFL", "CFL number", true);
    reader.AddOption(&order, "-O", "order", "Approximation order", true);
    reader.Parse();
    reader.PrintOptions();
}

void Solver::ReadMesh()
{
    int ix, iy, io, id[4];
    const int Gx = Nx + 1;
    const int Gy = Ny + 1;
    N_node = Gx * Gy;
    N_edge = Nx * Gy + Gx * Ny;
    N_elem = Nx * Ny;

    // allocation
    x_grid = NewVector<double>(Nx + 1);
    y_grid = NewVector<double>(Ny + 1);
    node_e = NewVector<NodeE>(N_node);
    node_u = NewVector<NodeU>(N_node);
    edge_e = NewVector<EdgeE>(N_edge);
    edge_u = NewVector<EdgeU>(N_edge);
    elem_e = NewVector<ElemE>(N_elem);
    elem_u = NewVector<ElemU>(N_elem);

    // make grid lines
    hx = (xmax - xmin) / (double)Nx;
    hy = (xmax - xmin) / (double)Ny;
    for (ix = 0; ix < Gx; ix++)
    {
        x_grid[ix] = xmin + (double)(ix)*hx;
    }
    for (iy = 0; iy < Gy; iy++)
    {
        y_grid[iy] = ymin + (double)(iy)*hy;
    }

    // make NodeE buffer
    for (iy = 0; iy < Gy; iy++)
    {
        for (ix = 0; ix < Gx; ix++)
        {
            io = Index2D(ix, iy, Gx, Gy);
            node_e[io].pos = {x_grid[ix], y_grid[iy]};
        }
    }

    // make EdgeE buffer
    for (iy = 0; iy < Gy; iy++)
    {
        for (ix = 0; ix < Nx; ix++)
        {
            io = Index2D(ix, iy, Nx, Gy);
            id[0] = Index2D(ix, iy, Gx, Gy);
            id[1] = Index2D(ix + 1, iy, Gx, Gy);
            if (iy % 2 == 0)
            {
                edge_e[io].beg = id[0];
                edge_e[io].end = id[1];
            }
            else
            {
                edge_e[io].beg = id[1];
                edge_e[io].end = id[0];
            }
        }
    }
    for (iy = 0; iy < Ny; iy++)
    {
        for (ix = 0; ix < Gx; ix++)
        {
            io = Nx * Gy + Index2D(ix, iy, Gx, Ny);
            id[0] = Index2D(ix, iy, Gx, Gy);
            id[1] = Index2D(ix, iy + 1, Gx, Gy);
            if (ix % 2 == 0)
            {
                edge_e[io].beg = id[1];
                edge_e[io].end = id[0];
            }
            else
            {
                edge_e[io].beg = id[0];
                edge_e[io].end = id[1];
            }
        }
    }

    // make ElemE buffer
    for (iy = 0; iy < Ny; iy++)
    {
        for (ix = 0; ix < Nx; ix++)
        {
            io = Index2D(ix, iy, Nx, Ny);
            elem_e[io].nvert = 4;
            elem_e[io].nedge = 4;
            elem_e[io].ndof = N_dof_loc;
            elem_e[io].vert_buf[0] = Index2D(ix, iy, Gx, Gy);
            elem_e[io].vert_buf[1] = Index2D(ix + 1, iy, Gx, Gy);
            elem_e[io].vert_buf[2] = Index2D(ix + 1, iy + 1, Gx, Gy);
            elem_e[io].vert_buf[3] = Index2D(ix, iy + 1, Gx, Gy);
            elem_e[io].edge_buf[0] = Index2D(ix, iy, Nx, Gy);
            elem_e[io].edge_buf[1] = Nx * Gy + Index2D(ix + 1, iy, Gx, Ny);
            elem_e[io].edge_buf[2] = Index2D(ix, iy + 1, Nx, Gy);
            elem_e[io].edge_buf[3] = Nx * Gy + Index2D(ix, iy, Gx, Ny);
        }
    }
}

void Solver::ReadDG()
{
    int ix, iy, io, j, k;
    double x[2];

    P = order - 1;
    // N_dof_loc = (P + 1) * (P + 1); // Qk polynomial
    N_dof_loc = (P + 2) * (P + 1) / 2; // Pk polynomial
    N_GL_1D = P + 1;
    x_GL_1D = NewVector<double>(N_GL_1D);
    w_GL_1D = NewVector<double>(N_GL_1D);
    GetGaussLegendre(N_GL_1D, x_GL_1D, w_GL_1D);

    N_GL_2D = (P + 1) * (P + 1);
    x_GL_2D = NewVector<double>(N_GL_2D);
    y_GL_2D = NewVector<double>(N_GL_2D);
    w_GL_2D = NewVector<double>(N_GL_2D);
    v_GL_2D = NewMatrix<double>(N_dof_loc, N_GL_2D);
    v_XL_GL_2D = NewMatrix<double>(N_dof_loc, N_GL_1D);
    v_XR_GL_2D = NewMatrix<double>(N_dof_loc, N_GL_1D);
    v_GL_YL_2D = NewMatrix<double>(N_dof_loc, N_GL_1D);
    v_GL_YR_2D = NewMatrix<double>(N_dof_loc, N_GL_1D);
    v_XC_YC_2D = NewVector<double>(N_dof_loc);
    v_XL_YL_2D = NewVector<double>(N_dof_loc);
    v_XR_YL_2D = NewVector<double>(N_dof_loc);
    v_XL_YR_2D = NewVector<double>(N_dof_loc);
    v_XR_YR_2D = NewVector<double>(N_dof_loc);
    GetGaussLegendre(2, P + 1, x_GL_2D, y_GL_2D, w_GL_2D);

    for (j = 0; j < N_dof_loc; j++)
    {
        for (k = 0; k < N_GL_2D; k++)
        {
            x[0] = x_GL_2D[k];
            x[1] = y_GL_2D[k];
            v_GL_2D[j][k] = CalcPkBasis(2, P, j, x);
        }
        for (k = 0; k < N_GL_1D; k++)
        {
            x[0] = 0.0;
            x[1] = x_GL_1D[k];
            v_XL_GL_2D[j][k] = CalcPkBasis(2, P, j, x);
            x[0] = 1.0;
            v_XR_GL_2D[j][k] = CalcPkBasis(2, P, j, x);
            x[0] = x_GL_1D[k];
            x[1] = 0.0;
            v_GL_YL_2D[j][k] = CalcPkBasis(2, P, j, x);
            x[1] = 1.0;
            v_GL_YR_2D[j][k] = CalcPkBasis(2, P, j, x);
        }
        x[0] = 0.5;
        x[1] = 0.5;
        v_XC_YC_2D[j] = CalcPkBasis(2, P, j, x);
        x[0] = 0.0;
        x[1] = 0.0;
        v_XL_YL_2D[j] = CalcPkBasis(2, P, j, x);
        x[0] = 1.0;
        v_XR_YL_2D[j] = CalcPkBasis(2, P, j, x);
        x[1] = 1.0;
        v_XR_YR_2D[j] = CalcPkBasis(2, P, j, x);
        x[0] = 0.0;
        v_XL_YR_2D[j] = CalcPkBasis(2, P, j, x);
    }

    for (iy = 0; iy < Ny; iy++)
    {
        for (ix = 0; ix < Nx; ix++)
        {
            io = Index2D(ix, iy, Nx, Ny);
            elem_e[io].ndof = N_dof_loc;
            elem_e[io].udof = NewVector<double>(N_dof_loc);
            elem_u[io].vcoe = NewMatrix<double>(N_dof_loc, N_dof_loc);
        }
    }

    N_dof_glo = N_elem * N_dof_loc;
    U_tn = NewVector<double>(N_dof_glo);
    U_tn1 = NewVector<double>(N_dof_glo);
    U_RHS = NewVector<double>(N_dof_glo);
}

void Solver::Simulate()
{
    int iter;
    double t, iter_time, cpu_time;

    Project();

    cpu_time = 0.;
    for (iter = 0, t = 0., ht = 0.; t < tmax; iter++, t += ht)
    {
        gtime.Start();

        // set time step
        if (SetTimeStep(t))
        {
            TrackBack();
            ReconstUShape();
            ReconstUTest();
            Clipping();
        }

        // evolute solution
        Step();

        iter_time = gtime.Stop();
        cpu_time += iter_time;

        std::cout << std::fixed;
        std::cout << "num of time steps: " << iter + 1 << ';'
                  << "  time step: " << std::setprecision(8) << ht << ';'
                  << "  curr time: " << std::setprecision(8) << t + ht << ';'
                  << "  iter time: " << std::setprecision(2) << iter_time * 1e3 << " ms" << ';'
                  << std::setprecision(8) << '\n';
    }
    std::cout << '\n';
    std::cout << "iteration terminate at time t = " << t << '\n';
    std::cout << "total cpu time = " << std::setprecision(2) << cpu_time * 1e3 << " ms" << '\n';
}

void Solver::Refine()
{
    hx *= 2;
    hy *= 2;
}

void Solver::CalcResidual()
{
    const int P_GL_new = 10;
    const int N_GL_new = P_GL_new * P_GL_new;
    int ix, iy, j, k, io;
    double x[2];
    double x_GL_new[N_GL_new];
    double y_GL_new[N_GL_new];
    double w_GL_new[N_GL_new];
    double v_GL_new[N_GL_new][N_dof_loc];
    double e_point;

    // preparation
    GetGaussLegendre(2, P_GL_new, x_GL_new, y_GL_new, w_GL_new);
    for (k = 0; k < N_GL_new; k++)
    {
        x[0] = x_GL_new[k], x[1] = y_GL_new[k];
        for (j = 0; j < N_dof_loc; j++)
        {
            v_GL_new[k][j] = CalcPkBasis(2, P, j, x);
        }
    }

    // calculate residuals
    e_L1 = e_L2 = e_Linf = 0.0;
    for (iy = 0; iy < Ny; iy++)
    {
        for (ix = 0; ix < Nx; ix++)
        {
            io = Index2D(ix, iy, Nx, Ny);
            for (k = 0; k < N_GL_new; k++)
            {
                x[0] = x_grid[ix] + hx * x_GL_new[k];
                x[1] = y_grid[iy] + hy * y_GL_new[k];
                e_point = 0.0;
                for (j = 0; j < N_dof_loc; j++)
                {
                    e_point += elem_e[io].udof[j] * v_GL_new[k][j];
                }

                e_point -= ExactData(tmax, x[0], x[1]);
                e_point = fabs(e_point);
                e_L1 += e_point * w_GL_new[k] * hx * hy;
                e_L2 += e_point * e_point * w_GL_new[k] * hx * hy;
                e_Linf = fmax(e_point, e_Linf);
            }
        }
    }
    e_L2 = sqrt(e_L2);
}

void Solver::GetResidual(const int nerr, double *e)
{
    if (nerr > 0)
    {
        e[0] = e_L1;
        if (nerr > 1)
        {
            e[1] = e_L2;
            if (nerr > 2)
            {
                e[2] = e_Linf;
            }
        }
    }
}

void Solver::OutputResidual(const char *file)
{
    std::ofstream out;
    out.open(file);
    out << std::fixed << std::setprecision(14);
    out << e_L1 << '\t' << e_L2 << '\t' << e_Linf << '\n';
    out.close();
}

void Solver::OutputSolution2D(const char *file)
{
    std::ofstream out;
    out.open(file);
    out << std::right << std::fixed << std::setprecision(8);
    int ix, iy, io, j;
    double xc, yc, uc, uc_exact;
    for (iy = 0; iy < Ny; iy++)
    {
        for (ix = 0; ix < Nx; ix++)
        {
            xc = x_grid[ix] + hx / 2.;
            yc = y_grid[iy] + hy / 2.;
            io = Index2D(ix, iy, Nx, Ny);
            uc = 0.;
            for (j = 0; j < N_dof_loc; j++)
            {
                uc += elem_e[io].udof[j] * v_XC_YC_2D[j];
            }

            uc_exact = ExactData(tmax, xc, yc);
            out << std::setw(13) << xc << '\t';
            out << std::setw(13) << yc << '\t';
            out << std::setw(13) << uc << '\t';
            out << std::setw(13) << uc_exact << '\n';
        }
    }
    out.close();
}

void Solver::PrintNodeE(const char *file)
{
    std::ofstream out;
    out.open(file);
    out << std::fixed << std::setprecision(8);
    for (int i = 0; i < N_node; i++)
    {
        Print(node_e[i], out);
    }
    out.close();
}

void Solver::PrintEdgeE(const char *file)
{
    std::ofstream out;
    out.open(file);
    for (int i = 0; i < N_edge; i++)
    {
        Print(edge_e[i], out);
    }
    out.close();
}

void Solver::PrintElemE(const char *file)
{
    std::ofstream out;
    out.open(file);
    for (int i = 0; i < N_elem; i++)
    {
        out << "[ElemE " << i << "]\n";
        Print(elem_e[i], out);
    }
    out.close();
}

void Solver::PrintNodeU(const char *file)
{
    std::ofstream out;
    out.open(file);
    out << std::fixed << std::setprecision(8);
    for (int i = 0; i < N_node; i++)
    {
        Print(node_u[i], out);
    }
    out.close();
}

void Solver::PrintEdgeU(const char *file)
{
    std::ofstream out;
    out.open(file);
    for (int i = 0; i < N_edge; i++)
    {
        Print(edge_u[i], out);
    }
    out.close();
}

void Solver::PrintElemU(const char *file)
{
    std::ofstream out;
    out.open(file);
    for (int i = 0; i < N_elem; i++)
    {
        out << "[ElemU " << i << "]\n";
        Print(elem_u[i], out);
    }
    out.close();
}

void Solver::Destroy()
{
    DeleteVector<double>(x_grid, Nx + 1);
    DeleteVector<double>(y_grid, Ny + 1);
    DeleteVector<double>(x_GL_1D, N_GL_1D);
    DeleteVector<double>(w_GL_1D, N_GL_1D);
    DeleteVector<double>(x_GL_2D, N_GL_2D);
    DeleteVector<double>(y_GL_2D, N_GL_2D);
    DeleteVector<double>(w_GL_2D, N_GL_2D);
    DeleteMatrix<double>(v_GL_2D, N_dof_loc, N_GL_2D);
    DeleteMatrix<double>(v_XL_GL_2D, N_dof_loc, N_GL_2D);
    DeleteMatrix<double>(v_XR_GL_2D, N_dof_loc, N_GL_2D);
    DeleteMatrix<double>(v_GL_YL_2D, N_dof_loc, N_GL_2D);
    DeleteMatrix<double>(v_GL_YR_2D, N_dof_loc, N_GL_2D);
    DeleteVector<double>(v_XC_YC_2D, N_dof_loc);
    DeleteVector<double>(v_XL_YL_2D, N_dof_loc);
    DeleteVector<double>(v_XR_YL_2D, N_dof_loc);
    DeleteVector<double>(v_XL_YR_2D, N_dof_loc);
    DeleteVector<double>(v_XR_YR_2D, N_dof_loc);
    DeleteVector<double>(U_tn, N_dof_glo);
    DeleteVector<double>(U_tn1, N_dof_glo);
    DeleteVector<double>(U_RHS, N_dof_glo);

    for (int i = 0; i < N_elem; i++)
    {
        DeleteVector<double>(elem_e[i].udof, N_dof_loc);
        DeleteMatrix<double>(elem_u[i].vcoe, N_dof_loc, N_dof_loc);
    }

    DeleteVector<NodeE>(node_e, N_node);
    DeleteVector<NodeU>(node_u, N_node);
    DeleteVector<EdgeE>(edge_e, N_edge);
    DeleteVector<EdgeU>(edge_u, N_edge);
    DeleteVector<ElemE>(elem_e, N_elem);
    DeleteVector<ElemU>(elem_u, N_elem);
}

double Solver::InitData(double x, double y)
{
    return sin(2. * M_PI * (x + y));
}

double Solver::ExactData(double t, double x, double y)
{
    return sin(2. * M_PI * (x + y - 2. * t));
}

void Solver::Project()
{
    // declarations
    int ix, iy, io, j, k;
    double a[kMAX_N_DOF];
    double f_GL_2D_k, x_GL_2D_k, y_GL_2D_k;

    // initialization: reference data
    for (j = 0; j < N_dof_loc; j++)
    {
        a[j] = 0.;
        for (k = 0; k < N_GL_2D; k++)
        {
            a[j] += w_GL_2D[k] * v_GL_2D[j][k] * v_GL_2D[j][k];
        }
    }

    // calculate DoFs for each element
    for (iy = 0; iy < Ny; iy++)
    {
        for (ix = 0; ix < Nx; ix++)
        {
            io = Index2D(ix, iy, Nx, Ny);
            for (j = 0; j < N_dof_loc; j++)
            {
                elem_e[io].udof[j] = 0.;
                for (k = 0; k < N_GL_2D; k++)
                {
                    x_GL_2D_k = x_grid[ix] + hx * x_GL_2D[k];
                    y_GL_2D_k = y_grid[iy] + hy * y_GL_2D[k];
                    f_GL_2D_k = InitData(x_GL_2D_k, y_GL_2D_k);
                    elem_e[io].udof[j] += w_GL_2D[k] * f_GL_2D_k * v_GL_2D[j][k];
                }
                elem_e[io].udof[j] /= a[j];
            }
        }
    }
}

int Solver::SetTimeStep(double t)
{
    double ht_old = ht;
    ht = CFL * fmax(hx, hy);
    ht = fmin(tmax - t, ht);
    return !FuzzyZero(ht - ht_old);
}

void Solver::TrackBack()
{
    // track NodeE to NodeU
    int i, jx, jy, k;
    for (i = 0; i < N_node; i++)
    {
        node_u[i].pos.x = node_e[i].pos.x - ax * ht;
        node_u[i].pos.y = node_e[i].pos.y - ay * ht;
        // periodic BC
        node_u[i].pos.x = Circulate(node_u[i].pos.x, xmin, xmax);
        node_u[i].pos.y = Circulate(node_u[i].pos.y, ymin, ymax);

        for (jx = 0; jx < Nx + 1; jx++)
        {
            if (node_u[i].pos.x < x_grid[jx])
            {
                break;
            }
        }

        for (jy = 0; jy < Ny + 1; jy++)
        {
            if (node_u[i].pos.y < y_grid[jy])
            {
                break;
            }
        }

        node_u[i].id = Index2D(jx - 1, jy - 1, Nx, Ny);
        if (node_u[i].id < 0 || node_u[i].id >= Nx * Ny)
        {
            node_u[i].id = -1;
        }
    }

    // track EdgeE to EdgeU
    for (i = 0; i < N_edge; i++)
    {
        edge_u[i].beg = edge_e[i].beg;
        edge_u[i].end = edge_e[i].end;
    }

    // track ElemE to ElemU
    for (i = 0; i < N_elem; i++)
    {
        elem_u[i].nvert = elem_e[i].nvert;
        elem_u[i].nedge = elem_e[i].nedge;
        for (k = 0; k < elem_u[i].nvert; k++)
        {
            elem_u[i].vert_buf[k] = elem_e[i].vert_buf[k];
        }
        for (k = 0; k < elem_u[i].nedge; k++)
        {
            elem_u[i].edge_buf[k] = elem_e[i].edge_buf[k];
        }
    }
}

void Solver::ReconstUShape()
{
    MCM_WARNING("ReconstUShape() is undefined!");
}

void Solver::ReconstUTest()
{
    int i, j, k;
    double v[4];
    double A[kMAX_N_DOF][kMAX_N_DOF];
    for (i = 0; i < N_elem; i++)
    {
        // solve reconstruction problem
        for (k = 0; k < N_dof_loc; k++)
        {
            A[0][k] = v_XL_YL_2D[k];
        }
        for (k = 0; k < N_dof_loc; k++)
        {
            A[1][k] = v_XR_YL_2D[k];
        }
        for (k = 0; k < N_dof_loc; k++)
        {
            A[2][k] = v_XR_YR_2D[k];
        }
        for (k = 0; k < N_dof_loc; k++)
        {
            A[3][k] = v_XL_YR_2D[k];
        }

        for (j = 0; j < N_dof_loc; j++)
        {
            v[0] = v_XL_YL_2D[j];
            v[1] = v_XR_YL_2D[j];
            v[2] = v_XR_YR_2D[j];
            v[3] = v_XL_YR_2D[j];
            PkReconst(P, (const double **)A, v, elem_u[i].vcoe[j]);
        }
    }
}

void Solver::Clipping()
{
    // find outer segments
    FindOuterSegments();

    // find inner segments
    FindInnerSegments();
}

void Solver::FindOuterSegments()
{
    const int M = Nx * (Ny + 1);
    int i, j, k[2], l[2], a, b;
    double x, y, xe[2], ye[2];
    POI p, q, r;
    p.type = 2;
    q.type = 2;
    r.type = 1;
    // horizontal edges
    for (i = 0; i < M; i++)
    {
        k[0] = edge_u[i].beg;
        k[1] = edge_u[i].end;
        l[0] = node_u[k[0]].id;
        l[1] = node_u[k[1]].id;
        xe[0] = node_u[k[0]].pos.x;
        xe[1] = node_u[k[1]].pos.x;
        y = node_u[k[0]].pos.y;

        edge_u[i].nsego = 2;
        a = (k[0] < k[1]);
        b = (l[0] < l[1]);

        // interior
        if (!(a ^ b))
        {
            for (j = 0; j < Nx + 1; j++)
            {
                x = x_grid[j];
                if ((x < xe[0]) ^ (x < xe[1]))
                {
                    p.id = l[0];
                    p.pos.x = xe[0];
                    p.pos.y = y;
                    r.id = l[0];
                    r.pos.x = x;
                    r.pos.y = y;
                    edge_u[i].sego_buf[0].beg = p;
                    edge_u[i].sego_buf[0].end = r;

                    r.id = l[1];
                    r.pos.x = x;
                    r.pos.y = y;
                    q.id = l[1];
                    q.pos.x = xe[1];
                    q.pos.y = y;
                    edge_u[i].sego_buf[1].beg = r;
                    edge_u[i].sego_buf[1].end = q;
                    break;
                }
            }
        }
        // boundary
        else
        {
            j = (a == 1) ? Nx : 0;
            x = x_grid[j];
            p.id = l[0];
            p.pos.x = xe[0];
            p.pos.y = y;
            r.id = l[0];
            r.pos.x = x;
            r.pos.y = y;
            edge_u[i].sego_buf[0].beg = p;
            edge_u[i].sego_buf[0].end = r;

            j = (a == 1) ? 0 : Nx;
            x = x_grid[j];
            r.id = l[1];
            r.pos.x = x;
            r.pos.y = y;
            q.id = l[1];
            q.pos.x = xe[1];
            q.pos.y = y;
            edge_u[i].sego_buf[1].beg = r;
            edge_u[i].sego_buf[1].end = q;
        }
    }

    // vertical edges
    for (i = M; i < N_edge; i++)
    {
        k[0] = edge_u[i].beg;
        k[1] = edge_u[i].end;
        l[0] = node_u[k[0]].id;
        l[1] = node_u[k[1]].id;
        ye[0] = node_u[k[0]].pos.y;
        ye[1] = node_u[k[1]].pos.y;
        x = node_u[k[0]].pos.x;

        edge_u[i].nsego = 2;
        a = (k[0] < k[1]);
        b = (l[0] < l[1]);

        // interior
        if (!(a ^ b))
        {
            for (j = 0; j < Ny + 1; j++)
            {
                y = y_grid[j];
                if ((y < ye[0]) ^ (y < ye[1]))
                {
                    p.type = 2;
                    p.id = l[0];
                    p.pos.x = x;
                    p.pos.y = ye[0];
                    r.type = 1;
                    r.id = l[0];
                    r.pos.x = x;
                    r.pos.y = y;
                    edge_u[i].sego_buf[0].beg = p;
                    edge_u[i].sego_buf[0].end = r;

                    r.type = 1;
                    r.id = l[1];
                    r.pos.x = x;
                    r.pos.y = y;
                    q.type = 2;
                    q.id = l[1];
                    q.pos.x = x;
                    q.pos.y = ye[1];
                    edge_u[i].sego_buf[1].beg = r;
                    edge_u[i].sego_buf[1].end = q;
                    break;
                }
            }
        }
        // boundary
        else
        {
            j = (a == 1) ? Ny : 0;
            y = y_grid[j];
            p.type = 2;
            p.id = l[0];
            p.pos.x = x;
            p.pos.y = ye[0];
            r.type = 1;
            r.id = l[0];
            r.pos.x = x;
            r.pos.y = y;
            edge_u[i].sego_buf[0].beg = p;
            edge_u[i].sego_buf[0].end = r;

            j = (a == 1) ? 0 : Nx;
            y = y_grid[j];
            r.type = 1;
            r.id = l[1];
            r.pos.x = x;
            r.pos.y = y;
            q.type = 2;
            q.id = l[1];
            q.pos.x = x;
            q.pos.y = ye[1];
            edge_u[i].sego_buf[1].beg = r;
            edge_u[i].sego_buf[1].end = q;
        }
    }
}

void Solver::FindInnerSegments()
{
    int i, j[4], k[4], l[4];
    double xe[2], ye[2];
    POI p, q, r;
    p.type = 1;
    q.type = 1;
    r.type = 0;
    for (i = 0; i < N_elem; i++)
    {
        j[0] = elem_u[i].vert_buf[0];
        j[1] = elem_u[i].vert_buf[1];
        j[2] = elem_u[i].vert_buf[2];
        j[3] = elem_u[i].vert_buf[3];
        k[0] = node_u[j[0]].id;
        k[1] = node_u[j[1]].id;
        k[2] = node_u[j[2]].id;
        k[3] = node_u[j[3]].id;
        l[0] = elem_u[i].edge_buf[0];
        l[1] = elem_u[i].edge_buf[1];
        l[2] = elem_u[i].edge_buf[2];
        l[3] = elem_u[i].edge_buf[3];
        elem_u[i].nsegi = 8;

        // left bottom
        xe[0] = edge_u[l[0]].sego_buf[0].beg.pos.x;
        xe[1] = edge_u[l[0]].sego_buf[0].end.pos.x;
        ye[0] = edge_u[l[3]].sego_buf[1].end.pos.y;
        ye[1] = edge_u[l[3]].sego_buf[1].beg.pos.y;
        p.id = k[0];
        r.id = k[0];
        q.id = k[0];
        p.pos.x = xe[1];
        p.pos.y = ye[0];
        r.pos.x = xe[1];
        r.pos.y = ye[1];
        q.pos.x = xe[0];
        q.pos.y = ye[1];
        elem_u[i].segi_buf[0].beg = p;
        elem_u[i].segi_buf[0].end = r;
        elem_u[i].segi_buf[1].beg = r;
        elem_u[i].segi_buf[1].end = q;

        // right bottom
        xe[0] = edge_u[l[0]].sego_buf[1].beg.pos.x;
        xe[1] = edge_u[l[0]].sego_buf[1].end.pos.x;
        ye[0] = edge_u[l[1]].sego_buf[0].beg.pos.y;
        ye[1] = edge_u[l[1]].sego_buf[0].end.pos.y;
        p.id = k[1];
        r.id = k[1];
        q.id = k[1];
        p.pos.x = xe[1];
        p.pos.y = ye[1];
        r.pos.x = xe[0];
        r.pos.y = ye[1];
        q.pos.x = xe[0];
        q.pos.y = ye[0];
        elem_u[i].segi_buf[2].beg = p;
        elem_u[i].segi_buf[2].end = r;
        elem_u[i].segi_buf[3].beg = r;
        elem_u[i].segi_buf[3].end = q;

        // right top
        xe[0] = edge_u[l[1]].sego_buf[1].end.pos.x;
        xe[1] = edge_u[l[1]].sego_buf[1].beg.pos.x;
        ye[0] = edge_u[l[2]].sego_buf[0].beg.pos.y;
        ye[1] = edge_u[l[2]].sego_buf[0].end.pos.y;
        p.id = k[2];
        r.id = k[2];
        q.id = k[2];
        p.pos.x = xe[0];
        p.pos.y = ye[1];
        r.pos.x = xe[0];
        r.pos.y = ye[0];
        q.pos.x = xe[1];
        q.pos.y = ye[0];
        elem_u[i].segi_buf[4].beg = p;
        elem_u[i].segi_buf[4].end = r;
        elem_u[i].segi_buf[5].beg = r;
        elem_u[i].segi_buf[5].end = q;

        // left top
        xe[0] = edge_u[l[2]].sego_buf[1].end.pos.x;
        xe[1] = edge_u[l[2]].sego_buf[1].beg.pos.x;
        ye[0] = edge_u[l[3]].sego_buf[0].beg.pos.y;
        ye[1] = edge_u[l[3]].sego_buf[0].end.pos.y;
        p.id = k[3];
        r.id = k[3];
        q.id = k[3];
        p.pos.x = xe[0];
        p.pos.y = ye[0];
        r.pos.x = xe[1];
        r.pos.y = ye[0];
        q.pos.x = xe[1];
        q.pos.y = ye[1];
        elem_u[i].segi_buf[6].beg = p;
        elem_u[i].segi_buf[6].end = r;
        elem_u[i].segi_buf[7].beg = r;
        elem_u[i].segi_buf[7].end = q;
    }
}

void Solver::Assemble()
{
    MCM_WARNING("Assemble() is undefined!");
}

void Solver::Step()
{
    SetZero(N_dof_glo, U_tn);
    SetZero(N_dof_glo, U_RHS);
    SetZero(N_dof_glo, U_tn1);

    CopyL2G(N_dof_glo, U_tn);

    Assemble();
    CopyG2G(N_dof_glo, U_RHS, U_tn1);

    CopyG2L(N_dof_glo, U_tn1);
}

void Solver::SetZero(const int size, double *v)
{
    std::fill(v, v + size, 0.);
}

void Solver::CopyL2G(const int size, double *v)
{
    for (int i = 0; i < N_elem; i++)
    {
        for (int j = 0; j < N_dof_loc; j++)
        {
            v[Index2D(i, j, N_elem, N_dof_loc)] = elem_e[i].udof[j];
        }
    }
}

void Solver::CopyG2L(const int size, const double *v)
{
    for (int i = 0; i < N_elem; i++)
    {
        for (int j = 0; j < N_dof_loc; j++)
        {
            elem_e[i].udof[j] = v[Index2D(i, j, N_elem, N_dof_loc)];
        }
    }
}

void Solver::CopyG2G(const int size, const double *v1, double *v2)
{
    std::copy(v1, v1 + size, v2);
}

void Solver::Print(NodeE &node, std::ostream &out)
{
    out << ' ' << std::setw(12) << node.pos.x
        << ' ' << std::setw(12) << node.pos.y;
    out << '\n';
}

void Solver::Print(EdgeE &edge, std::ostream &out)
{
    out << ' ' << std::setw(6) << edge.beg
        << ' ' << std::setw(6) << edge.end;
    out << '\n';
}

void Solver::Print(ElemE &elem, std::ostream &out)
{
    int i;
    out << "vertices:\n";
    for (i = 0; i < 4; i++)
    {
        out << ' ' << std::setw(6) << elem.vert_buf[i];
    }
    out << '\n';
    out << "edges:\n";
    for (i = 0; i < 4; i++)
    {
        out << ' ' << std::setw(6) << elem.edge_buf[i];
    }
    out << '\n';
    out << "DoFs:\n";
    out << std::fixed;
    for (i = 0; i < elem.ndof; i++)
    {
        out << std::setw(12) << std::setprecision(8) << elem.udof[i] << '\n';
    }
    out << '\n';
}

void Solver::Print(NodeU &node, std::ostream &out)
{
    out << ' ' << std::setw(6) << node.id
        << ' ' << std::setw(12) << node.pos.x
        << ' ' << std::setw(12) << node.pos.y;
    out << '\n';
}

void Solver::Print(EdgeU &edge, std::ostream &out)
{
    out << ' ' << std::setw(6) << edge.beg
        << ' ' << std::setw(6) << edge.end
        << ' ' << std::setw(6) << edge.nsego
        << ' ' << std::setw(6) << edge.sego_buf[0].beg.pos.x
        << ' ' << std::setw(6) << edge.sego_buf[0].beg.pos.y
        << ' ' << std::setw(6) << edge.sego_buf[0].end.pos.x
        << ' ' << std::setw(6) << edge.sego_buf[0].end.pos.y
        << ' ' << std::setw(6) << edge.sego_buf[1].beg.pos.x
        << ' ' << std::setw(6) << edge.sego_buf[1].beg.pos.y
        << ' ' << std::setw(6) << edge.sego_buf[1].end.pos.x
        << ' ' << std::setw(6) << edge.sego_buf[1].end.pos.y;
    out << '\n';
}

void Solver::Print(ElemU &elem, std::ostream &out)
{
    int i;
    out << "vertices:\n";
    for (i = 0; i < 4; i++)
    {
        out <<  ' ' << std::setw(6) << elem.vert_buf[i];
    }
    out << '\n';
    out << "edges:\n";
    for (i = 0; i < 4; i++)
    {
        out << ' ' << std::setw(6) << elem.edge_buf[i];
    }
    out << '\n';
    out << "inner segments:\n";
    for (i = 0; i < elem.nsegi; i++)
    {
        out << ' ' << std::setw(6) << elem.segi_buf[i].beg.pos.x
            << ' ' << std::setw(6) << elem.segi_buf[i].beg.pos.y
            << ' ' << std::setw(6) << elem.segi_buf[i].end.pos.x
            << ' ' << std::setw(6) << elem.segi_buf[i].end.pos.y;
        out << '\n';
    }
    out << '\n';
}