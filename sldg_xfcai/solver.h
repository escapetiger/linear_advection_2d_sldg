#ifndef _SOLVER_H_
#define _SOLVER_H_

#include "globals.h"

const double gau2[2][2] = {{-sqrt(3.) / 6., 0.5},
                           {sqrt(3.) / 6., 0.5}};
const double gau3[3][2] = {{-sqrt(15.) / 10., 5. / 18.},
                           {0., 4. / 9.},
                           {sqrt(15.) / 10., 5. / 18.}};

const double ai[6] = {1., 12., 12., 180., 144., 180.};

const double xg[6] = {
    -0.466234757101576013906150777246997304567,
    -0.330604693233132256830699797509952673503,
    -0.119309593041598454315250860840355967709,
    0.119309593041598454315250860840355967709,
    0.330604693233132256830699797509952673503,
    0.466234757101576013906150777246997304567};

const double wg[6] = {
    1.71324492379170345040296142172733e-1 / 2.,
    3.60761573048138607569833513837716e-1 / 2.,
    4.67913934572691047389870343989551e-1 / 2.,
    4.67913934572691047389870343989551e-1 / 2.,
    3.60761573048138607569833513837716e-1 / 2.,
    1.71324492379170345040296142172733e-1 / 2.,
};

class Solver
{
public:
    double xleft = -M_PI;
    double xright = M_PI;
    double ybottom = -M_PI;
    double ytop = M_PI;
    int nx = 4;
    int ny = 4;
    int nghost = 5;
    int nmoment = 6;
    double time_final = 1.5;
    double cfl = 2.;
    int iexample = 1;
    int irk = 3;
    int relative = 1;
    int idebug = 0;
    double dx = (xright - xleft) / (double)nx;
    double dy = (ytop - ybottom) / (double)ny;

    double t;
    double *xgrid;
    double *ygrid;
    double *x;
    double *y;
    NodeE **vertex;
    NodeU **vertex_star;
    NodeE **nodex;
    NodeE **nodey;
    NodeE **nodec;
    NodeU **nodex_star;
    NodeU **nodey_star;
    NodeU **nodec_star;

    Face **face_lr;
    Face **face_bt;

    CellE **element;
    CellU **element_star;

    double ***umod_t;
    double **com_mass;

    double VelX(double x, double y, double t)
    {
        switch (iexample)
        {
        case 1:
            return 1.;
        case 2:
            return -cos(x / 2.) * cos(x / 2.) * sin(y) * cos(M_PI * t / time_final) * M_PI;
        case 3:
            return -y;
        }
        return 0.;
    }
    double VelY(double x, double y, double t)
    {
        switch (iexample)
        {
        case 1:
            return 1.;
        case 2:
            return cos(y / 2.) * cos(y / 2.) * sin(x) * cos(M_PI * t / time_final) * M_PI;
        case 3:
            return x;
        }
        return 0.;
    }
    double ExactData(double x, double y, double t)
    {
        double rb, rb0;
        switch (iexample)
        {
        case 1:
            return sin(x + y - 2. * t);
        case 2:
            rb = sqrt(pow(x - 0.15 / 0.5 * M_PI, 2.) + pow(y, 2.));
            rb0 = 0.3 * M_PI;
            if (rb < rb0)
            {
                return rb0 * pow(cos(M_PI * rb / (2. * rb0)), 6.);
            }
            else
            {
                return 0.;
            }
        case 3:
            return exp(-(x * x + y * y));
        }
        return 0.;
    }

    double Phi(int k, double x, double y)
    {
        switch (k)
        {
        case 0:
            return 1.;
        case 1:
            return x;
        case 2:
            return y;
        case 3:
            return x * x - 1. / 12.;
        case 4:
            return x * y;
        case 5:
            return y * y - 1. / 12.;
        }
        return 0.;
    }
    double CalcPoly(double *a00, double x16, double xc16, double dx16, double y16, double yc16, double dy16, int k)
    {
        switch (k)
        {
        case 0:
            return a00[0];
        case 1:
            return a00[0] + a00[1] * (x16 - xc16) / dx16 + a00[2] * (y16 - yc16) / dy16;
        case 2:
            return a00[0] + a00[1] * (x16 - xc16) / dx16 + a00[2] * (y16 - yc16) / dy16 + a00[3] * (pow((x16 - xc16) / dx16, 2.) - 1. / 12.) + a00[4] * ((x16 - xc16) / dx16 * (y16 - yc16) / dy16) + a00[5] * (pow((y16 - yc16) / dy16, 2.) - 1. / 12.);
        }
        return 0.;
    }

    void Setup();
    void Allocate();
    void Init();
    // void Simulate();
    void Output(const char *file);
    void Debug();
    void Destroy();
};

#endif