#include "solver.h"

void Solver::Setup()
{
    nx = 40;
    ny = 40;
    nghost = 5;
    nmoment = 6;
    time_final = 1.5;
    cfl = 2.;
    iexample = 1;
    irk = 3;
    relative = 1;
    idebug = 0;
    if (iexample == 1)
    {
        xleft = -M_PI;
        xright = M_PI;
        ybottom = -M_PI;
        ytop = M_PI;
    }
    else if (iexample == 2)
    {
        xleft = -M_PI;
        xright = M_PI;
        ybottom = -M_PI;
        ytop = M_PI;
    }
    else if (iexample == 3)
    {
        xleft = -M_PI * 2.;
        xright = M_PI * 2.;
        ybottom = -M_PI * 2.;
        ytop = M_PI * 2.;
    }

    dx = (xright - xleft) / (double)nx;
    dy = (ytop - ybottom) / (double)ny;

    std::cout << "iexample   = " << iexample << '\n';
    std::cout << "xleft      = " << xleft << '\n';
    std::cout << "xright     = " << xright << '\n';
    std::cout << "ybottom    = " << ybottom << '\n';
    std::cout << "ytop       = " << ytop << '\n';
    std::cout << "time_final = " << time_final << '\n';
    std::cout << "nx         = " << nx << '\n';
    std::cout << "ny         = " << ny << '\n';
    std::cout << "cfl        = " << cfl << '\n';
    std::cout << "dx         = " << dx << '\n';
    std::cout << "dy         = " << dy << '\n';
    std::cout << "irk        = " << irk << '\n';
}

void Solver::Allocate()
{
    xgrid = NewVector<double>(nx + 1 + 2 * nghost);
    ygrid = NewVector<double>(ny + 1 + 2 * nghost);

    x = NewVector<double>(nx + 2 * nghost);
    y = NewVector<double>(ny + 2 * nghost);

    vertex = NewMatrix<NodeE>(nx + 1, ny + 1);
    vertex_star = NewMatrix<NodeU>(nx + 1, ny + 1);

    nodex = NewMatrix<NodeE>(nx, ny + 1);
    nodey = NewMatrix<NodeE>(nx + 1, ny);
    nodec = NewMatrix<NodeE>(nx, ny);
    nodex_star = NewMatrix<NodeU>(nx, ny + 1);
    nodey_star = NewMatrix<NodeU>(nx + 1, ny);
    nodec_star = NewMatrix<NodeU>(nx, ny);

    face_lr = NewMatrix<Face>(nx, ny + 1);
    face_bt = NewMatrix<Face>(nx + 1, ny);

    element = NewMatrix<CellE>(nx + 2 * nghost, ny + 2 * nghost);
    element_star = NewMatrix<CellU>(nx, ny);

    for (int i = 0; i < nx + 2 * nghost; i++)
    {
        for (int j = 0; j < ny + 2 * nghost; j++)
        {
            element[i][j].umodal = NewVector<double>(nmoment);
        }
    }
    umod_t = NewTensor3D<double>(nx, ny, nmoment);
    com_mass = NewMatrix<double>(nx, ny);
}

void Solver::Destroy()
{
    DeleteVector<double>(xgrid, nx + 1 + 2 * nghost);
    DeleteVector<double>(ygrid, nx + 1 + 2 * nghost);

    DeleteVector<double>(x, nx + 2 * nghost);
    DeleteVector<double>(y, ny + 2 * nghost);

    DeleteMatrix<NodeE>(vertex, nx + 1, ny + 1);
    DeleteMatrix<NodeU>(vertex_star, nx + 1, ny + 1);

    DeleteMatrix<NodeE>(nodex, nx, ny + 1);
    DeleteMatrix<NodeE>(nodey, nx + 1, ny);
    DeleteMatrix<NodeE>(nodec, nx, ny);
    DeleteMatrix<NodeU>(nodex_star, nx, ny + 1);
    DeleteMatrix<NodeU>(nodey_star, nx + 1, ny);
    DeleteMatrix<NodeU>(nodec_star, nx, ny);

    DeleteMatrix<Face>(face_lr, nx, ny + 1);
    DeleteMatrix<Face>(face_bt, nx + 1, ny);

    for (int i = 0; i < nx + 2 * nghost; i++)
    {
        for (int j = 0; j < ny + 2 * nghost; j++)
        {
            DeleteVector<double>(element[i][j].umodal, nmoment);
        }
    }
    DeleteMatrix<CellE>(element, nx + 2 * nghost, ny + 2 * nghost);
    DeleteMatrix<CellU>(element_star, nx, ny);

    DeleteTensor3D<double>(umod_t, nx, ny, nmoment);
    DeleteMatrix<double>(com_mass, nx, ny);
}

void Solver::Init()
{
    int i, j, k, lx, ly;

    t = 0.;

    for (i = -nghost; i < nx + nghost; i++)
    {
        x[nghost + i] = xleft + ((double)i + 0.5) * dx;
    }
    for (j = -nghost; j < ny + nghost; j++)
    {
        y[nghost + j] = ybottom + ((double)j + 0.5) * dy;
    }
    for (i = -nghost; i < nx + 1 + nghost; i++)
    {
        xgrid[nghost + i] = xleft + (double)i * dx;
    }
    for (j = -nghost; j < ny + 1 + nghost; j++)
    {
        ygrid[nghost + j] = ybottom + (double)j * dy;
    }
    for (i = 0; i < nx + 1; i++)
    {
        for (j = 0; j < ny + 1; j++)
        {
            vertex[i][j].pos.x = xleft + (double)i * dx;
            vertex[i][j].pos.y = ybottom + (double)j * dy;
        }
    }
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny + 1; j++)
        {
            nodex[i][j].pos.x = xleft + ((double)i + 0.5) * dx;
            nodex[i][j].pos.y = ybottom + (double)j * dy;
        }
    }
    for (i = 0; i < nx + 1; i++)
    {
        for (j = 0; j < ny; j++)
        {
            nodey[i][j].pos.x = xleft + (double)i * dx;
            nodey[i][j].pos.y = ybottom + ((double)j + 0.5) * dy;
        }
    }
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            nodec[i][j].pos.x = xleft + ((double)i + 0.5) * dx;
            nodec[i][j].pos.y = ybottom + ((double)j + 0.5) * dy;
        }
    }
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            element[i][j].vertex1 = &vertex[i][j];
            element[i][j].vertex2 = &vertex[i + 1][j];
            element[i][j].vertex3 = &vertex[i + 1][j + 1];
            element[i][j].vertex4 = &vertex[i][j + 1];
            if (nmoment == 6)
            {
                element[i][j].vertex5 = &nodex[i][j];
                element[i][j].vertex6 = &nodey[i + 1][j];
                element[i][j].vertex7 = &nodex[i][j + 1];
                element[i][j].vertex8 = &nodey[i][j];
                element[i][j].vertex9 = &nodec[i][j];
            }
        }
    }
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            for (k = 0; k < nmoment; k++)
            {
                double utmp = 0.;
                for (lx = 0; lx < 6; lx++)
                {
                    for (ly = 0; ly < 6; ly++)
                    {
                        utmp += ExactData(x[i] + xg[lx] * dx, y[j] + xg[ly] * dy, 0.) * Phi(k, xg[lx], xg[ly]) * wg[lx] * wg[ly];
                    }
                }
                element[i][j].umodal[k] = utmp * ai[k];
            }
        }
    }
}

void Solver::Debug()
{
    std::ofstream out("elem.dat");
    int i, j;
    double x, y;
    out << std::fixed << std::setprecision(8);
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            out << "[Vertex (" << i << "," << j << ") ]" << '\n';
            x = element[i][j].vertex1->pos.x;
            y = element[i][j].vertex1->pos.y;
            out << x << ' ' << y << '\n';
            x = element[i][j].vertex2->pos.x;
            y = element[i][j].vertex2->pos.y;
            out << x << ' ' << y << '\n';
            x = element[i][j].vertex3->pos.x;
            y = element[i][j].vertex3->pos.y;
            out << x << ' ' << y << '\n';
            x = element[i][j].vertex4->pos.x;
            y = element[i][j].vertex4->pos.y;
            out << x << ' ' << y << '\n';
        }
    }
    out.close();
}

void Solver::Output(const char *file)
{
    std::ofstream out;
    out.open(file);

    int i, j;
    double xc, yc, uc, uce;
    out << std::fixed << std::setprecision(8);
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            xc = nodec[i][j].pos.x;
            yc = nodec[i][j].pos.y;
            uc = CalcPoly(element[i][j].umodal, xc, xc, dx, yc, yc, dy, 2);
            uce = ExactData(xc, yc, t);
            out << std::setw(13) << xc << ' '
                << std::setw(13) << yc << ' '
                << std::setw(13) << uc << ' '
                << std::setw(13) << uce << '\n';
        }
    }

    out.close();
}