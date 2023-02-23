#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double Circulate(double x, double xL, double xR)
{
    // z = Circulate(x, xL, xR), z in [xL, xR]
    // y = x - xL = k*P + r,
    // where P = xR - xL and r = (x - xL) % (xR - xL) in (-P,P);

    double P = xR - xL;
    double y = x - xL;
    double r = fmod(y, P);

    // 0 < y < P, z = x;
    if (y > 0. && y < P)
    {
        return x;
    }
    // y <= 0, r == 0, z = xL;
    if (y <= 0. && fabs(r) < 1e-8)
    {
        return xL;
    }
    // y <= 0, r < 0, z = xR + r;
    if (y <= 0. && r < 0.)
    {
        return xR + r;
    }
    // y >= P, r == 0, z = xR;
    if (y >= P && fabs(r) < 1e-8)
    {
        return xR;
    }
    // y >= P, r > 0, z = xL + r;
    if (y >= P && r > 0.)
    {
        return xL + r;
    }
    // otherwise, do nothing
    return x;
}

int Circulate(int x, int xL, int xR)
{
    // z = Circulate(x, xL, xR), z in [xL, xR)
    // y = x - xL = k*P + r,
    // where P = xR - xL and r = (x - xL) % (xR - xL) in (-P,P);
    int P = xR - xL;
    int y = x - xL;
    int r = y % P;

    // 0 <= y < P, z = x
    if (y >= 0 && y < P)
    {
        return x;
    }
    // y < 0, r == 0, z = xL
    if (y < 0 && r == 0)
    {
        return xL;
    }
    // y < 0, r < 0, z = xR + r
    if (y < 0 && r < 0)
    {
        return xR + r;
    }
    // y >= P, r == 0, z = xL
    if (y >= P && r == 0)
    {
        return xL;
    }
    // y >= P, r > 0, z = xL + r
    if (y >= P && r > 0)
    {
        return xL + r;
    }
    // otherwise, do nothing
    return x;
}

int main(int argc, char const *argv[])
{
    double xmin = 0.;
    double xmax = 1.;
    double x0, x1, hx, x0_u, x1_u, x_u, ax, ht, r;

    ax = 100000;
    ht = -.0001;
    hx = .0001;
    x0 = 1 * hx;
    x1 = 2 * hx;
    hx = x1 - x0;
    x0_u = x0 - ax * ht;
    x1_u = x1 - ax * ht;
    printf("x0 = %.4f\n", x0);
    printf("x1 = %.4f\n", x1);
    printf("x0_u = %.4f\n", x0_u);
    printf("x1_u = %.4f\n", x1_u);

    r = fmod(x1_u, hx);
    // x1_u-xmin = k*h+r, -h<r<h
    // r < 0, x = k*h-h = x1_u-xmin-r-h
    // r > 0, x = k*h = x1_u-xmin-r
    // r == 0, x = x1_u-xmin
    if (r < 0)
    {
        x_u = x1_u - xmin - r - hx;
    }
    else if (r > 0)
    {
        x_u = x1_u - xmin - r;
    }
    else
    {
        x_u = x1_u - xmin;
    }

    printf("r = %.4f\n", r);
    printf("x = %.4f\n", x_u);
    printf("check upper bound: %d\n", x_u <= x1_u);
    printf("check lower bound: %d\n", x_u >= x0_u);

    return 0;
}
