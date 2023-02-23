#include "solver.h"

using mcm::DeleteMatrix;
using mcm::DeleteVector;
using mcm::Index2D;
using mcm::NewMatrix;
using mcm::NewVector;
using mcm::kernels::FuzzyEQ;
using mcm::kernels::FuzzyGE;
using mcm::kernels::FuzzyGT;
using mcm::kernels::FuzzyLE;
using mcm::kernels::FuzzyLT;
using mcm::kernels::FuzzyZero;
using mcm::kernels::RealGE;
using mcm::kernels::RealGT;
using mcm::kernels::RealLE;
using mcm::kernels::RealLT;
using mcm::kernels::Square;

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

void Solver::ReadOptions(SolverOptions *opt)
{
    xmin = opt->xmin;
    xmax = opt->xmax;
    ymin = opt->ymin;
    ymax = opt->ymax;
    tmax = opt->tmax;
    ax = opt->ax;
    ay = opt->ay;
    Nx = opt->Nx;
    Ny = opt->Ny;
    CFL = opt->CFL;
    order = opt->order;
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
            elem_e[io].p0.x = node_e[elem_e[io].vert_buf[0]].pos.x;
            elem_e[io].p0.y = node_e[elem_e[io].vert_buf[0]].pos.y;
            elem_e[io].p1.x = node_e[elem_e[io].vert_buf[2]].pos.x;
            elem_e[io].p1.y = node_e[elem_e[io].vert_buf[2]].pos.y;
        }
    }
}

void Solver::ReadDG()
{
    int ix, iy, io, j, k, l;
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

    mat_v_u_2D = NewMatrix<double>(N_dof_loc, N_dof_loc);
    for (j = 0; j < N_dof_loc; j++) // i : test index
    {
        for (k = 0; k < N_dof_loc; k++) // j : trial index
        {
            mat_v_u_2D[j][k] = 0.;
            for (l = 0; l < N_GL_2D; l++)
            {
                mat_v_u_2D[j][k] += v_GL_2D[j][l] * w_GL_2D[l] * v_GL_2D[k][l];
            }
        }
    }

    // DoFs
    for (iy = 0; iy < Ny; iy++)
    {
        for (ix = 0; ix < Nx; ix++)
        {
            io = Index2D(ix, iy, Nx, Ny);
            elem_e[io].ndof = N_dof_loc;
            elem_e[io].udof = NewVector<double>(N_dof_loc);
        }
    }

    N_dof_glo = N_elem * N_dof_loc;
    U_tn = NewVector<double>(N_dof_glo);
    U_tn1 = NewVector<double>(N_dof_glo);
    U_RHS = NewVector<double>(N_dof_glo);

#ifdef DEBUG
    using mcm::PrintMatrix;
    using mcm::PrintVector;
    PrintVector(x_GL_2D, N_GL_2D, "x_GL_2D = ");
    PrintVector(y_GL_2D, N_GL_2D, "y_GL_2D = ");
    PrintVector(w_GL_2D, N_GL_2D, "w_GL_2D = ");
    PrintMatrix(v_GL_2D, N_dof_loc, N_GL_2D, "v_GL_2D = ");
    PrintVector(v_XC_YC_2D, N_dof_loc, "v_XC_YC_2D = ");
    PrintMatrix(mat_v_u_2D, N_dof_loc, N_dof_loc, "mat_v_u_2D = ");
#endif
}

void Solver::Simulate()
{
    int iter;
    double t, iter_time, cpu_time;

    Project();

    cpu_time = 0.;
    for (iter = 0, t = 0., ht = 0.; FuzzyLT(t, tmax); iter++, t += ht)
    {
        gtime.Start();

        // set time step
        if (SetTimeStep(t))
        {
            TrackBack();
            Clipping();
            SetQuadRule();
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
        // node_u[i].pos.x = Circulate(node_u[i].pos.x, xmin, xmax);
        // node_u[i].pos.y = Circulate(node_u[i].pos.y, ymin, ymax);

        for (jx = 0; jx < Nx + 1; jx++)
        {
            if (FuzzyLT(node_u[i].pos.x, x_grid[jx]))
            {
                break;
            }
        }

        for (jy = 0; jy < Ny + 1; jy++)
        {
            if (FuzzyLT(node_u[i].pos.y, y_grid[jy]))
            {
                break;
            }
        }

        node_u[i].par = Index2D(jx - 1, jy - 1, Nx, Ny);
        if (node_u[i].par < 0 || node_u[i].par >= Nx * Ny)
        {
            node_u[i].par = -1;
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
        elem_u[i].p0.x = elem_e[i].p0.x - ax * ht;
        elem_u[i].p0.y = elem_e[i].p0.y - ay * ht;
        elem_u[i].p1.x = elem_e[i].p1.x - ax * ht;
        elem_u[i].p1.y = elem_e[i].p1.y - ay * ht;
    }
}

void Solver::Clipping()
{
    int i;
    for (i = 0; i < N_elem; i++)
    {
        ClipElemU(&elem_u[i]);
        ClipElemE(&elem_e[i], &elem_u[i]);
        FinalizeClip(i, &elem_e[i], &elem_u[i]);
    }
}

void Solver::ClipElemU(ElemU *eu)
{
    double x0, y0, x1, y1, rx, ry, x, y, dx, dy;
    x0 = eu->p0.x;
    y0 = eu->p0.y;
    x1 = eu->p1.x;
    y1 = eu->p1.y;
    // x1_u-xmin = k*h+r, -h<r<h
    // r < 0, x = k*h-h = x1_u-xmin-r-h
    // r > 0, x = k*h = x1_u-xmin-r
    // r == 0, x = x1_u-xmin
    rx = fmod(x1, hx);
    if (rx < 0)
    {
        x = x1 - xmin - rx - hx;
    }
    else if (rx > 0)
    {
        x = x1 - xmin - rx;
    }
    else
    {
        x = x1 - xmin;
    }
    ry = fmod(y1, hy);
    if (ry < 0)
    {
        y = y1 - ymin - ry - hy;
    }
    else if (ry > 0)
    {
        y = y1 - ymin - ry;
    }
    else
    {
        y = y1 - ymin;
    }
    eu->poi.x = (x - x0) / hx;
    eu->poi.y = (y - y0) / hy;

    // add SubElems
    eu->nsub = 0;
    // left bottom SubElem
    dx = fabs(x - x0);
    dy = fabs(y - y0);
    if (!(FuzzyZero(dx) || FuzzyZero(dy)))
    {
        eu->sub_buf[eu->nsub].p0.x = x0;
        eu->sub_buf[eu->nsub].p0.y = y0;
        eu->sub_buf[eu->nsub].p1.x = x;
        eu->sub_buf[eu->nsub].p1.y = y;
        eu->nsub++;
    }
    // right bottom SubElem
    dx = fabs(x1 - x);
    dy = fabs(y - y0);
    if (!(FuzzyZero(dx) || FuzzyZero(dy)))
    {
        eu->sub_buf[eu->nsub].p0.x = x;
        eu->sub_buf[eu->nsub].p0.y = y0;
        eu->sub_buf[eu->nsub].p1.x = x1;
        eu->sub_buf[eu->nsub].p1.y = y;
        eu->nsub++;
    }
    // right top SubElem
    dx = fabs(x1 - x);
    dy = fabs(y1 - y);
    if (!(FuzzyZero(dx) || FuzzyZero(dy)))
    {
        eu->sub_buf[eu->nsub].p0.x = x;
        eu->sub_buf[eu->nsub].p0.y = y;
        eu->sub_buf[eu->nsub].p1.x = x1;
        eu->sub_buf[eu->nsub].p1.y = y1;
        eu->nsub++;
    }
    // left top SubElem
    dx = fabs(x - x0);
    dy = fabs(y1 - y);
    if (!(FuzzyZero(dx) || FuzzyZero(dy)))
    {
        eu->sub_buf[eu->nsub].p0.x = x0;
        eu->sub_buf[eu->nsub].p0.y = y;
        eu->sub_buf[eu->nsub].p1.x = x;
        eu->sub_buf[eu->nsub].p1.y = y1;
        eu->nsub++;
    }
}

void Solver::ClipElemE(ElemE *ee, ElemU *eu)
{
    double x0, y0, x1, y1, x, y, dx, dy;
    x0 = ee->p0.x;
    y0 = ee->p0.y;
    x1 = ee->p1.x;
    y1 = ee->p1.y;
    ee->poi.x = eu->poi.x;
    ee->poi.y = eu->poi.y;
    x = ee->p0.x + ee->poi.x * hx;
    y = ee->p0.y + ee->poi.y * hy;

    // add SubElems
    ee->nsub = 0;
    // left bottom SubElem
    dx = fabs(x - x0);
    dy = fabs(y - y0);
    if (!(FuzzyZero(dx) || FuzzyZero(dy)))
    {
        ee->sub_buf[ee->nsub].p0.x = x0;
        ee->sub_buf[ee->nsub].p0.y = y0;
        ee->sub_buf[ee->nsub].p1.x = x;
        ee->sub_buf[ee->nsub].p1.y = y;
        ee->nsub++;
    }
    // right bottom SubElem
    dx = fabs(x1 - x);
    dy = fabs(y - y0);
    if (!(FuzzyZero(dx) || FuzzyZero(dy)))
    {
        ee->sub_buf[ee->nsub].p0.x = x;
        ee->sub_buf[ee->nsub].p0.y = y0;
        ee->sub_buf[ee->nsub].p1.x = x1;
        ee->sub_buf[ee->nsub].p1.y = y;
        ee->nsub++;
    }
    // right top SubElem
    dx = fabs(x1 - x);
    dy = fabs(y1 - y);
    if (!(FuzzyZero(dx) || FuzzyZero(dy)))
    {
        ee->sub_buf[ee->nsub].p0.x = x;
        ee->sub_buf[ee->nsub].p0.y = y;
        ee->sub_buf[ee->nsub].p1.x = x1;
        ee->sub_buf[ee->nsub].p1.y = y1;
        ee->nsub++;
    }
    // left top SubElem
    dx = fabs(x - x0);
    dy = fabs(y1 - y);
    if (!(FuzzyZero(dx) || FuzzyZero(dy)))
    {
        ee->sub_buf[ee->nsub].p0.x = x0;
        ee->sub_buf[ee->nsub].p0.y = y;
        ee->sub_buf[ee->nsub].p1.x = x;
        ee->sub_buf[ee->nsub].p1.y = y1;
        ee->nsub++;
    }
}

void Solver::FinalizeClip(int par_e, ElemE *ee, ElemU *eu)
{
    int i, jx, jy;
    double cx, cy;
    for (i = 0; i < ee->nsub; i++)
    {
        ee->sub_buf[i].par = par_e;
    }
    for (i = 0; i < eu->nsub; i++)
    {
        cx = (eu->sub_buf[i].p0.x + eu->sub_buf[i].p1.x) / 2.;
        cy = (eu->sub_buf[i].p0.y + eu->sub_buf[i].p1.y) / 2.;
        cx = Circulate(cx, xmin, xmax);
        cy = Circulate(cy, ymin, ymax);
        for (jx = 0; jx < Nx + 1; jx++)
        {
            if (FuzzyLT(cx, x_grid[jx]))
            {
                break;
            }
        }

        for (jy = 0; jy < Ny + 1; jy++)
        {
            if (FuzzyLT(cy, y_grid[jy]))
            {
                break;
            }
        }

        eu->sub_buf[i].par = Index2D(jx - 1, jy - 1, Nx, Ny);
    }
}

void Solver::SetQuadRule()
{
    int i, j;
    for (i = 0; i < N_elem; i++)
    {
        for (j = 0; j < 4; j++)
        {
            MakeSubQuadRuleE(&elem_e[i], &elem_e[i].sub_buf[j]);
            MakeSubQuadRuleU(&elem_u[i], &elem_u[i].sub_buf[j]);
        }
    }
}

void Solver::MakeSubQuadRuleE(ElemE *elem, SubElem *sub)
{
    int i;
    double dx, dy, dxy;
    dx = fabs(sub->p1.x - sub->p0.x);
    dy = fabs(sub->p1.y - sub->p0.y);
    dxy = dx * dy;

    for (i = 0; i < N_GL_2D; i++)
    {
        sub->qp[i].x = sub->p0.x + dx * x_GL_2D[i];
        sub->qp[i].y = sub->p0.y + dy * y_GL_2D[i];
        sub->qp[i].w = w_GL_2D[i] * dxy;
    }
}

void Solver::MakeSubQuadRuleU(ElemU *elem, SubElem *sub)
{
    int i;
    double dx, dy, dxy;
    dx = fabs(sub->p1.x - sub->p0.x);
    dy = fabs(sub->p1.y - sub->p0.y);
    dxy = dx * dy;

    for (i = 0; i < N_GL_2D; i++)
    {
        sub->qp[i].x = sub->p0.x + dx * x_GL_2D[i];
        sub->qp[i].y = sub->p0.y + dy * y_GL_2D[i];
        sub->qp[i].w = w_GL_2D[i] * dxy;

        // periodic BC
        sub->qp[i].x = Circulate(sub->qp[i].x, xmin, xmax);
        sub->qp[i].y = Circulate(sub->qp[i].y, ymin, ymax);
    }
}

void Solver::Assemble()
{
    int i, j, k, l, m, jg;
    double u, v, udof, mass;
    QuadPoint *qe, *qu;
    for (i = 0; i < N_elem; i++)
    {
        for (j = 0; j < N_dof_loc; j++)
        {
            jg = Index2D(i, j, N_elem, N_dof_loc);
            mass = hx * hy * mat_v_u_2D[j][j];
            for (l = 0; l < elem_e[i].nsub; l++)
            {
                for (m = 0; m < N_GL_2D; m++)
                {
                    qe = &elem_e[i].sub_buf[l].qp[m];
                    qu = &elem_u[i].sub_buf[l].qp[m];
                    u = 0.;
                    v = CalcElemEPkBasis(j, &elem_e[i], qe->xy);
                    for (k = 0; k < N_dof_loc; k++)
                    {
                        udof = elem_e[elem_u[i].sub_buf[l].par].udof[k];
                        u += udof * CalcElemUPkBasis(k, &elem_u[i], &elem_u[i].sub_buf[l], qu->xy);
                    }
                    // std::cout << std::fixed << std::setprecision(8);
                    // std::cout << i << ' ' << j << ' ' << u << ' ' << v << ' ' << qe->w << '\n';
                    U_RHS[jg] += qe->w * u * v / mass;
                }
            }
        }
    }
}

double Solver::CalcElemEPkBasis(const int k, ElemE *elem, double *x)
{
    bool outside = (FuzzyLT(x[0], elem->p0.x)) || (FuzzyGT(x[0], elem->p1.x)) ||
                   (FuzzyLT(x[1], elem->p0.y)) || (FuzzyGT(x[1], elem->p1.y));
    // bool interface = (FuzzyEQ(x[0], elem->x0)) || (FuzzyEQ(x[0], elem->x1)) ||
    //                  (FuzzyEQ(x[1], elem->y0)) || (FuzzyEQ(x[1], elem->y1));
    if (outside)
    {
        return 0.;
    }

    double xx[2] = {(x[0] - elem->p0.x) / hx, (x[1] - elem->p0.y) / hy};

    return CalcPkBasis(2, P, k, xx);
}

double Solver::CalcElemUPkBasis(const int k, ElemU *elem, SubElem *sub, double *x)
{
    return CalcElemEPkBasis(k, &elem_e[sub->par], x);
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

void Solver::CopyG2L(int size, const double *v)
{
    for (int i = 0; i < N_elem; i++)
    {
        for (int j = 0; j < N_dof_loc; j++)
        {
            elem_e[i].udof[j] = v[Index2D(i, j, N_elem, N_dof_loc)];
        }
    }
}

void Solver::CopyG2G(int size, const double *v1, double *v2)
{
    std::copy(v1, v1 + size, v2);
}
