#include <iostream>
#include <cmath>

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
    double xL = -1.;
    double xR = 1.;
    double x = -1.5;
    std::cout << Circulate(x, xL, xR) << std::endl;
    x = .5;
    std::cout << Circulate(x, xL, xR) << std::endl;
    x = 1.5;
    std::cout << Circulate(x, xL, xR) << std::endl;
    x = -1.;
    std::cout << Circulate(x, xL, xR) << std::endl;
    x = 1.;
    std::cout << Circulate(x, xL, xR) << std::endl;
    x = -3.;
    std::cout << Circulate(x, xL, xR) << std::endl;
    x = 3.;
    std::cout << Circulate(x, xL, xR) << std::endl;

    int ia = 0;
    int ib = 10;
    int i = 2;
    std::cout << Circulate(i, ia, ib) << std::endl;
    i = -2;
    std::cout << Circulate(i, ia, ib) << std::endl;
    i = 12;
    std::cout << Circulate(i, ia, ib) << std::endl;
    i = 0;
    std::cout << Circulate(i, ia, ib) << std::endl;
    i = 10;
    std::cout << Circulate(i, ia, ib) << std::endl;
    i = 20;
    std::cout << Circulate(i, ia, ib) << std::endl;
    i = -10;
    std::cout << Circulate(i, ia, ib) << std::endl;

    return 0;
}
