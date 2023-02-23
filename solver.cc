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
    node_e_buf = NewVector<NodeE>(N_node);
    node_u_buf = NewVector<NodeU>(N_node);
    edge_e_buf = NewVector<EdgeE>(N_edge);
    edge_u_buf = NewVector<EdgeU>(N_edge);
    elem_e_buf = NewVector<ElemE>(N_elem);
    elem_u_buf = NewVector<ElemU>(N_elem);

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
            node_e_buf[io].pos = {x_grid[ix], y_grid[iy]};
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
                edge_e_buf[io].beg = id[0];
                edge_e_buf[io].end = id[1];
            }
            else
            {
                edge_e_buf[io].beg = id[1];
                edge_e_buf[io].end = id[0];
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
                edge_e_buf[io].beg = id[1];
                edge_e_buf[io].end = id[0];
            }
            else
            {
                edge_e_buf[io].beg = id[0];
                edge_e_buf[io].end = id[1];
            }
        }
    }

    // make ElemE buffer
    for (iy = 0; iy < Ny; iy++)
    {
        for (ix = 0; ix < Nx; ix++)
        {
            io = Index2D(ix, iy, Nx, Ny);
            elem_e_buf[io].nvert = 4;
            elem_e_buf[io].nedge = 4;
            elem_e_buf[io].ndof = N_dof_loc;
            elem_e_buf[io].vert_buf[0] = Index2D(ix, iy, Gx, Gy);
            elem_e_buf[io].vert_buf[1] = Index2D(ix + 1, iy, Gx, Gy);
            elem_e_buf[io].vert_buf[2] = Index2D(ix + 1, iy + 1, Gx, Gy);
            elem_e_buf[io].vert_buf[3] = Index2D(ix, iy + 1, Gx, Gy);
            elem_e_buf[io].edge_buf[0] = Index2D(ix, iy, Nx, Gy);
            elem_e_buf[io].edge_buf[1] = Nx * Gy + Index2D(ix + 1, iy, Gx, Ny);
            elem_e_buf[io].edge_buf[2] = Index2D(ix, iy + 1, Nx, Gy);
            elem_e_buf[io].edge_buf[3] = Nx * Gy + Index2D(ix, iy, Gx, Ny);
            elem_e_buf[io].p0.x = node_e_buf[elem_e_buf[io].vert_buf[0]].pos.x;
            elem_e_buf[io].p0.y = node_e_buf[elem_e_buf[io].vert_buf[0]].pos.y;
            elem_e_buf[io].p1.x = node_e_buf[elem_e_buf[io].vert_buf[2]].pos.x;
            elem_e_buf[io].p1.y = node_e_buf[elem_e_buf[io].vert_buf[2]].pos.y;
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
            elem_e_buf[io].ndof = N_dof_loc;
            elem_e_buf[io].udof = NewVector<double>(N_dof_loc);
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
                    e_point += elem_e_buf[io].udof[j] * v_GL_new[k][j];
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
                uc += elem_e_buf[io].udof[j] * v_XC_YC_2D[j];
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
        DeleteVector<double>(elem_e_buf[i].udof, N_dof_loc);
    }

    DeleteVector<NodeE>(node_e_buf, N_node);
    DeleteVector<NodeU>(node_u_buf, N_node);
    DeleteVector<EdgeE>(edge_e_buf, N_edge);
    DeleteVector<EdgeU>(edge_u_buf, N_edge);
    DeleteVector<ElemE>(elem_e_buf, N_elem);
    DeleteVector<ElemU>(elem_u_buf, N_elem);
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
    double a[MAX_DOF];
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
                elem_e_buf[io].udof[j] = 0.;
                for (k = 0; k < N_GL_2D; k++)
                {
                    x_GL_2D_k = x_grid[ix] + hx * x_GL_2D[k];
                    y_GL_2D_k = y_grid[iy] + hy * y_GL_2D[k];
                    f_GL_2D_k = InitData(x_GL_2D_k, y_GL_2D_k);
                    elem_e_buf[io].udof[j] += w_GL_2D[k] * f_GL_2D_k * v_GL_2D[j][k];
                }
                elem_e_buf[io].udof[j] /= a[j];
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
        node_u_buf[i].pos.x = node_e_buf[i].pos.x - ax * ht;
        node_u_buf[i].pos.y = node_e_buf[i].pos.y - ay * ht;

        // periodic BC
        // node_u[i].pos.x = Circulate(node_u[i].pos.x, xmin, xmax);
        // node_u[i].pos.y = Circulate(node_u[i].pos.y, ymin, ymax);

        for (jx = 0; jx < Nx + 1; jx++)
        {
            if (FuzzyLT(node_u_buf[i].pos.x, x_grid[jx]))
            {
                break;
            }
        }

        for (jy = 0; jy < Ny + 1; jy++)
        {
            if (FuzzyLT(node_u_buf[i].pos.y, y_grid[jy]))
            {
                break;
            }
        }

        node_u_buf[i].par = Index2D(jx - 1, jy - 1, Nx, Ny);
        if (node_u_buf[i].par < 0 || node_u_buf[i].par >= Nx * Ny)
        {
            node_u_buf[i].par = -1;
        }
    }

    // track EdgeE to EdgeU
    for (i = 0; i < N_edge; i++)
    {
        edge_u_buf[i].beg = edge_e_buf[i].beg;
        edge_u_buf[i].end = edge_e_buf[i].end;
    }

    // track ElemE to ElemU
    for (i = 0; i < N_elem; i++)
    {
        elem_u_buf[i].nvert = elem_e_buf[i].nvert;
        elem_u_buf[i].nedge = elem_e_buf[i].nedge;
        for (k = 0; k < elem_u_buf[i].nvert; k++)
        {
            elem_u_buf[i].vert_buf[k] = elem_e_buf[i].vert_buf[k];
        }
        for (k = 0; k < elem_u_buf[i].nedge; k++)
        {
            elem_u_buf[i].edge_buf[k] = elem_e_buf[i].edge_buf[k];
        }
        elem_u_buf[i].p0.x = elem_e_buf[i].p0.x - ax * ht;
        elem_u_buf[i].p0.y = elem_e_buf[i].p0.y - ay * ht;
        elem_u_buf[i].p1.x = elem_e_buf[i].p1.x - ax * ht;
        elem_u_buf[i].p1.y = elem_e_buf[i].p1.y - ay * ht;
    }
}

void Solver::Clipping()
{
    int i;
    ClipREF();
    for (i = 0; i < N_elem; i++)
    {
        // clip ElemU and set QuadRuleU
        ClipElem(&elem_u_buf[i]);
        // clip ElemE and set QuadRuleE
        ClipElem(&elem_e_buf[i]);
    }
#ifdef FAST_SLDG
    MakeQuadRuleREF(); // fast computing
#endif
}

void Solver::ClipREF()
{
    double x0, y0, x1, y1, rx, ry, x, y, dx, dy;
    x0 = ref_buf.p0.x - ax * ht / hx;
    y0 = ref_buf.p0.y - ay * ht / hy;
    x1 = ref_buf.p1.x - ax * ht / hx;
    y1 = ref_buf.p1.y - ay * ht / hy;
    // x0 = i*hx-a*ht, x1 = (i+1)*hx-a*ht;
    // xx0 = -a*ht/hx, xx1 = 1-a*ht/hx;
    // 1 = k+r, -1<r<1
    // r < 0, xx = k-1 = -r
    // r > 0, xx = k = 1-r
    // r = 0, xx = 1
    // we have x0 < x = x0+xx*hx < x1.
    rx = fmod(x1, 1.);
    x = (rx < 0) ? -rx : 1. - rx;
    ry = fmod(y1, 1.);
    y = (ry < 0) ? -ry : 1. - ry;

    ref_buf.poi.x = x;
    ref_buf.poi.y = y;

    // add SubElems
    ref_buf.nsub = 0;
    x0 = ref_buf.p0.x;
    y0 = ref_buf.p0.y;
    x1 = ref_buf.p1.x;
    y1 = ref_buf.p1.y;
    // left bottom SubElem
    dx = fabs(x - x0);
    dy = fabs(y - y0);
    if (!(FuzzyZero(dx) || FuzzyZero(dy)))
    {
        ref_buf.sub_buf[ref_buf.nsub].p0.x = x0;
        ref_buf.sub_buf[ref_buf.nsub].p0.y = y0;
        ref_buf.sub_buf[ref_buf.nsub].p1.x = x;
        ref_buf.sub_buf[ref_buf.nsub].p1.y = y;
        ref_buf.nsub++;
    }
    // right bottom SubElem
    dx = fabs(x1 - x);
    dy = fabs(y - y0);
    if (!(FuzzyZero(dx) || FuzzyZero(dy)))
    {
        ref_buf.sub_buf[ref_buf.nsub].p0.x = x;
        ref_buf.sub_buf[ref_buf.nsub].p0.y = y0;
        ref_buf.sub_buf[ref_buf.nsub].p1.x = x1;
        ref_buf.sub_buf[ref_buf.nsub].p1.y = y;
        ref_buf.nsub++;
    }
    // right top SubElem
    dx = fabs(x1 - x);
    dy = fabs(y1 - y);
    if (!(FuzzyZero(dx) || FuzzyZero(dy)))
    {
        ref_buf.sub_buf[ref_buf.nsub].p0.x = x;
        ref_buf.sub_buf[ref_buf.nsub].p0.y = y;
        ref_buf.sub_buf[ref_buf.nsub].p1.x = x1;
        ref_buf.sub_buf[ref_buf.nsub].p1.y = y1;
        ref_buf.nsub++;
    }
    // left top SubElem
    dx = fabs(x - x0);
    dy = fabs(y1 - y);
    if (!(FuzzyZero(dx) || FuzzyZero(dy)))
    {
        ref_buf.sub_buf[ref_buf.nsub].p0.x = x0;
        ref_buf.sub_buf[ref_buf.nsub].p0.y = y;
        ref_buf.sub_buf[ref_buf.nsub].p1.x = x;
        ref_buf.sub_buf[ref_buf.nsub].p1.y = y1;
        ref_buf.nsub++;
    }
}

void Solver::ClipElem(ElemU *elem_u)
{
    int i, jx, jy, k;
    double x0, y0, x1, y1, cx, cy, dx, dy, dxy;
    // find parents
    for (i = 0; i < ref_buf.nsub; i++)
    {
        x0 = elem_u->p0.x + ref_buf.sub_buf[i].p0.x * hx;
        y0 = elem_u->p0.y + ref_buf.sub_buf[i].p0.y * hy;
        x1 = elem_u->p0.x + ref_buf.sub_buf[i].p1.x * hx;
        y1 = elem_u->p0.y + ref_buf.sub_buf[i].p1.y * hy;
        dx = fabs(x1 - x0);
        dy = fabs(y1 - y0);
        dxy = dx * dy;
        cx = (x0 + x1) / 2.;
        cy = (y0 + y1) / 2.;
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

        elem_u->qr.par_buf[i] = Index2D(jx - 1, jy - 1, Nx, Ny);

        for (k = 0; k < N_GL_2D; k++)
        {
            // This is for periodic BC
            elem_u->qr.qp_buf[i][k].x = Circulate(x0 + dx * x_GL_2D[k], xmin, xmax);
            elem_u->qr.qp_buf[i][k].y = Circulate(y0 + dy * y_GL_2D[k], xmin, xmax);

            // This is for other BC
            // elem_u->qr.qp_buf[i][k].x = x0 + dx * x_GL_2D[k];
            // elem_u->qr.qp_buf[i][k].y = y0 + dy * y_GL_2D[k];
            elem_u->qr.qp_buf[i][k].w = w_GL_2D[k] * dxy;
        }
    }
}

void Solver::ClipElem(ElemE *elem_e)
{
    int i, j;
    double x0, y0, x1, y1, dx, dy, dxy;
    for (i = 0; i < ref_buf.nsub; i++)
    {
        x0 = elem_e->p0.x + ref_buf.sub_buf[i].p0.x * hx;
        y0 = elem_e->p0.y + ref_buf.sub_buf[i].p0.y * hy;
        x1 = elem_e->p0.x + ref_buf.sub_buf[i].p1.x * hx;
        y1 = elem_e->p0.y + ref_buf.sub_buf[i].p1.y * hy;
        dx = fabs(x1 - x0);
        dy = fabs(y1 - y0);
        dxy = dx * dy;
        for (j = 0; j < N_GL_2D; j++)
        {
            elem_e->qr.qp_buf[i][j].x = x0 + dx * x_GL_2D[j];
            elem_e->qr.qp_buf[i][j].y = y0 + dy * y_GL_2D[j];
            elem_e->qr.qp_buf[i][j].w = w_GL_2D[j] * dxy;
        }
    }
}

void Solver::MakeQuadRuleREF()
{
    int i, j, k;
    double x0, y0;
    ElemE *elem_e = &elem_e_buf[0];
    ElemU *elem_u = &elem_u_buf[0];
    for (i = 0; i < ref_buf.nsub; i++)
    {
        // reference QuadRule for ElemE
        x0 = elem_e->p0.x;
        y0 = elem_e->p0.y;
        for (j = 0; j < N_GL_2D; j++)
        {
            ref_buf.qp_e_buf[i][j].x = (elem_e->qr.qp_buf[i][j].x - x0) / hx;
            ref_buf.qp_e_buf[i][j].y = (elem_e->qr.qp_buf[i][j].y - y0) / hy;
            ref_buf.qp_e_buf[i][j].w = elem_e->qr.qp_buf[i][j].w / (hx * hy);
        }

        // reference QuadRule for ElemU
        x0 = elem_e[elem_u->qr.par_buf[i]].p0.x;
        y0 = elem_e[elem_u->qr.par_buf[i]].p0.x;
        for (j = 0; j < N_GL_2D; j++)
        {
            ref_buf.qp_u_buf[i][j].x = (elem_u->qr.qp_buf[i][j].x - x0) / hx;
            ref_buf.qp_u_buf[i][j].y = (elem_u->qr.qp_buf[i][j].y - y0) / hy;
            ref_buf.qp_u_buf[i][j].w = elem_u->qr.qp_buf[i][j].w / (hx * hy);
        }

        // reference basis values for ElemE
        for (j = 0; j < N_dof_loc; j++)
        {
            for (k = 0; k < N_GL_2D; k++)
            {
                ref_buf.v_qp_e_buf[i][j][k] = CalcPkBasis(2, P, j, ref_buf.qp_e_buf[i][k].xy);
            }
        }

        // reference basis values for ElemU
        for (j = 0; j < N_dof_loc; j++)
        {
            for (k = 0; k < N_GL_2D; k++)
            {
                ref_buf.v_qp_u_buf[i][j][k] = CalcPkBasis(2, P, j, ref_buf.qp_u_buf[i][k].xy);
            }
        }
    }
}

void Solver::Assemble()
{
    int i, j, k, l, m, jg;
    double u, v, udof, mass;
    for (i = 0; i < N_elem; i++)
    {
        for (j = 0; j < N_dof_loc; j++)
        {
            jg = Index2D(i, j, N_elem, N_dof_loc);
            mass = hx * hy * mat_v_u_2D[j][j];
            for (l = 0; l < ref_buf.nsub; l++)
            {
                for (m = 0; m < N_GL_2D; m++)
                {
#ifndef FAST_SLDG
                    QuadPoint *qe, *qu;
                    qe = &elem_e_buf[i].qr.qp_buf[l][m];
                    qu = &elem_u_buf[i].qr.qp_buf[l][m];
                    u = 0.;
                    v = CalcElemEPkBasis(j, qe->xy, &elem_e_buf[i]);
                    for (k = 0; k < N_dof_loc; k++)
                    {
                        udof = elem_e_buf[elem_u_buf[i].qr.par_buf[l]].udof[k];
                        u += udof * CalcElemUPkBasis(k, qu->xy, &elem_u_buf[i], elem_u_buf[i].qr.par_buf[l]);
                    }
                    U_RHS[jg] += qe->w * u * v / mass;
#else
                    u = 0.;
                    v = ref_buf.v_qp_e_buf[l][j][m];
                    for (k = 0; k < N_dof_loc; k++)
                    {
                        udof = elem_e_buf[elem_u_buf[i].qr.par_buf[l]].udof[k];
                        u += udof * ref_buf.v_qp_u_buf[l][k][m];
                    }
                    u *= hx * hy;
                    U_RHS[jg] += ref_buf.qp_e_buf[l][m].w * u * v / mass;
#endif
                    // std::cout << std::fixed << std::setprecision(8);
                    // std::cout << i << ' ' << j << ' ' << u << ' ' << v << ' ' << qe->w << '\n';
                }
            }
        }
    }
}

double Solver::CalcElemEPkBasis(int k, double *x, ElemE *elem)
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

double Solver::CalcElemUPkBasis(int k, double *x, ElemU *elem, int sub_par)
{
    return CalcElemEPkBasis(k, x, &elem_e_buf[sub_par]);
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
            v[Index2D(i, j, N_elem, N_dof_loc)] = elem_e_buf[i].udof[j];
        }
    }
}

void Solver::CopyG2L(int size, const double *v)
{
    for (int i = 0; i < N_elem; i++)
    {
        for (int j = 0; j < N_dof_loc; j++)
        {
            elem_e_buf[i].udof[j] = v[Index2D(i, j, N_elem, N_dof_loc)];
        }
    }
}

void Solver::CopyG2G(int size, const double *v1, double *v2)
{
    std::copy(v1, v1 + size, v2);
}
