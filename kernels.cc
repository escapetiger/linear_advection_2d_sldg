#include "kernels.h"

void SetQuadPoint(const int i, const double a, const double v,
                  double *x, double *w)
{
    x[i] = a;
    w[i] = v;
}

void GetGaussLegendre(const int n, double *x, double *w)
{
    double a, v;
    switch (n)
    {
    case 1:
        SetQuadPoint(0, 0.5, 1.0, x, w);
        return;
    case 2:
        a = (1.0 - sqrt(1.0 / 3.0)) / 2.0;
        v = 0.5;
        SetQuadPoint(0, a, v, x, w);
        SetQuadPoint(1, 1.0 - a, v, x, w);
        return;
    case 3:
        a = (1.0 - sqrt(3.0 / 5.0)) / 2.0;
        v = 5.0 / 18.0;
        SetQuadPoint(0, a, v, x, w);
        a = 0.5;
        v = 4.0 / 9.0;
        SetQuadPoint(1, a, v, x, w);
        a = (1.0 + sqrt(3.0 / 5.0)) / 2.0;
        v = 5.0 / 18.0;
        SetQuadPoint(2, a, v, x, w);
        return;
    case 4:
        a = (1.0 - sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0))) / 2.0;
        v = (18.0 - sqrt(30.0)) / 72.0;
        SetQuadPoint(0, a, v, x, w);
        a = (1.0 - sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0))) / 2.0;
        v = (18.0 + sqrt(30.0)) / 72.0;
        SetQuadPoint(1, a, v, x, w);
        a = (1.0 + sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0))) / 2.0;
        v = (18.0 + sqrt(30.0)) / 72.0;
        SetQuadPoint(2, a, v, x, w);
        a = (1.0 + sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0))) / 2.0;
        v = (18.0 - sqrt(30.0)) / 72.0;
        SetQuadPoint(3, a, v, x, w);
        return;
    case 5:
        a = (1.0 - sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0)) / 3.0) / 2.0;
        v = (322.0 - 13.0 * sqrt(70.0)) / 1800.0;
        SetQuadPoint(0, a, v, x, w);
        a = (1.0 - sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0)) / 3.0) / 2.0;
        v = (322.0 + 13.0 * sqrt(70.0)) / 1800.0;
        SetQuadPoint(1, a, v, x, w);
        a = 0.5;
        v = 64.0 / 225.0;
        SetQuadPoint(2, a, v, x, w);
        a = (1.0 + sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0)) / 3.0) / 2.0;
        v = (322.0 + 13.0 * sqrt(70.0)) / 1800.0;
        SetQuadPoint(3, a, v, x, w);
        a = (1.0 + sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0)) / 3.0) / 2.0;
        v = (322.0 - 13.0 * sqrt(70.0)) / 1800.0;
        SetQuadPoint(4, a, v, x, w);
        return;
    case 6:
        a = 0.033765242898424;
        v = 0.085662246189585;
        SetQuadPoint(0, a, v, x, w);
        a = 0.169395306766868;
        v = 0.180380786524069;
        SetQuadPoint(1, a, v, x, w);
        a = 0.380690406958402;
        v = 0.233956967286345;
        SetQuadPoint(2, a, v, x, w);
        a = 0.619309593041598;
        v = 0.233956967286345;
        SetQuadPoint(3, a, v, x, w);
        a = 0.830604693233132;
        v = 0.180380786524069;
        SetQuadPoint(4, a, v, x, w);
        a = 0.966234757101576;
        v = 0.085662246189585;
        SetQuadPoint(5, a, v, x, w);
        return;
    }

    const int m = (n + 1.0) / 2.0;

    for (int i = 1; i <= m; i++)
    {
        double z = cos(M_PI * (i - 0.25) / (n + 0.5));
        double pp, p1, dz, xi = 0.0;
        bool done = false;
        while (1)
        {
            double p2 = 1.0;
            p1 = z;
            for (int j = 2; j <= n; j++)
            {
                double p3 = p2;
                p2 = p1;
                p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j;
            }
            // p1 is Legendre polynomial

            pp = n * (z * p1 - p2) / (z * z - 1);
            if (done)
            {
                break;
            }

            dz = p1 / pp;
            if (fabs(dz) < 1e-16)
            {
                done = true;
                // map the new point (z-dz) to (0,1):
                xi = ((1 - z) + dz) / 2; // (1 - (z - dz))/2 has bad round-off
                                         // continue the computation: get pp at the new point, then exit
            }
            // update: z = z - dz
            z -= dz;
        }

        double weight = 1.0 / (4.0 * xi * (1.0 - xi) * pp * pp);
        SetQuadPoint(i - 1, xi, weight, x, w);
        SetQuadPoint(n - i, 1.0 - xi, weight, x, w);
    }
}

void GetGaussLegendre(const int dim, const int n, double *x, double *y, double *w)
{
    using mcm::Index2D;

    int i, j, o;
    double x_GL_1D[20];
    double w_GL_1D[20];
    GetGaussLegendre(n, x_GL_1D, w_GL_1D);

    for (j = 0; j < n; j++)
    {
        for (i = 0; i < n; i++)
        {
            o = Index2D(i, j, n);
            x[o] = x_GL_1D[i];
            y[o] = x_GL_1D[j];
            w[o] = w_GL_1D[i] * w_GL_1D[j];
        }
    }
}

double Legendre1D(const int n, const int m, const double x)
{
    const double xi = 2.0 * x - 1.0;
    int n_row = m + 1;
    int n_col = n + 1;
    double s[n_row * n_col];

    // Store L_0, L_1, ..., L_{n} at the first row
    s[0] = 1.0;
    s[1] = xi;
    for (int j = 2; j < n_col; j++)
    {
        s[j] = ((2 * j - 1) * xi * s[j - 1] - (j - 1) * s[j - 2]) / j;
    }

    // Compute derivatives
    for (int i = 1; i < n_row; i++)
    {
        s[i * n_col] = 0.0;
        s[1 + i * n_col] = (1 == i) ? 1.0 : 0.0;
        for (int j = 2; j < n_col; j++)
        {
            // L_n^{(m)}(x) = (2n-1)L_{n-1}^{(m-1)}(x) + L_{n-2}^{(m)}(x)
            s[j + i * n_col] = (2 * j - 1) * s[j - 1 + (i - 1) * n_col] + s[j - 2 + i * n_col];
        }
    }

    return s[n_row * n_col - 1] * pow(2.0, (double)m);
}

double CalcQkBasis(const int dim, const int k, const int *n, const double *x)
{
    MCM_CONTRACT_VAR(k);
    MCM_ASSERT(k < 3, "it is recommended to use k < 3");
    double y = 1.0;
    for (int i = 0; i < dim; i++)
    {
        y *= Legendre1D(n[i], 0, x[i]);
    }
    return y;
}

double CalcQkBasis(const int dim, const int k, const int j, const int *n, const double *x)
{
    MCM_CONTRACT_VAR(k);
    MCM_ASSERT(k < 3, "it is recommended to use k < 3");
    int m;
    double y = 1.0;
    for (int i = 0; i < dim; i++)
    {
        m = (i == j) ? 1 : 0;
        y *= Legendre1D(n[i], m, x[i]);
    }
    return y;
}

double CalcQkBasis(const int dim, const int k, const int n, const double *x)
{
    MCM_CONTRACT_VAR(dim);
    MCM_CONTRACT_VAR(k);
    MCM_ASSERT(dim == 2, "dim must be equal to 2");
    MCM_ASSERT(k < 3, "k must be less than 3");
    double xx[2] = {x[0] - .5, x[1] - .5};
    switch (n)
    {
    case 0:
        return 1.;
    case 1:
        return xx[0];
    case 2:
        return xx[0] * xx[0] - 1. / 12.;
    case 3:
        return xx[1];
    case 4:
        return xx[0] * xx[1];
    case 5:
        return (xx[0] * xx[0] - 1. / 12.) * xx[1];
    case 6:
        return xx[1] * xx[1] - 1. / 12.;
    case 7:
        return xx[0] * (xx[1] * xx[1] - 1. / 12.);
    case 8:
        return (xx[0] * xx[0] - 1. / 12.) * (xx[1] * xx[1] - 1. / 12.);
    }
    return 0.;
}

double CalcPkBasis(const int dim, const int k, const int n, const double *x)
{
    MCM_CONTRACT_VAR(dim);
    MCM_CONTRACT_VAR(k);
    MCM_ASSERT(dim == 2, "dim must be equal to 2");
    MCM_ASSERT(k < 3, "k must be less than 3");
    double xx[2] = {x[0] - .5, x[1] - .5};
    switch (n)
    {
    case 0:
        return 1.;
    case 1:
        return xx[0];
    case 2:
        return xx[1];
    case 3:
        return xx[0] * xx[0] - 1. / 12.;
    case 4:
        return xx[0] * xx[1];
    case 5:
        return xx[1] * xx[1] - 1. / 12.;
    }
    return 0.;
}

double CalcWeightedSum(const int n, const double *w, const double *x)
{
    double y = 0.;
    for (int i = 0; i < n; i++)
    {
        y += w[i] * x[i];
    }
    return y;
}

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

void PkReconst(const int k, const double **Ap, const double *bp, double *xp)
{
    if (k == 0)
    {
        xp[0] = bp[0];
    }
    if (k == 1)
    {
        int i, j;
        Eigen::Matrix<double, 4, 3> A;
        for (i = 0; i < 4; i++)
        {
            for (j = 0; j < 3; j++)
            {
                A(i, j) = Ap[i][j];
            }
        }

        Eigen::Vector3d b;
        for (i = 0; i < 4; i++)
        {
            b[i] = bp[i];
        }

        // Eigen::Vector3d x = (A.transpose() * A).ldlt().solve(A.transpose() * b);

        // Eigen::Map<Eigen::Vector3d>(xp, 3) = x;
    }
}

bool CollinearAndIntersect(const int dim, const double *p, const double *q, const double *r)
{
    if (q[0] <= std::max(p[0], r[0]) && q[0] >= std::min(p[0], r[0]) &&
        q[1] <= std::max(p[1], r[1]) && q[1] >= std::min(p[1], r[1]))
        return true;

    return false;
}

int Orientation(const int dim, const double *p, const double *q, const double *r)
{
    // See https://www.geeksforgeeks.org/orientation-3-ordered-points/
    // for details of below formula.
    double val = (q[1] - p[1]) * (r[0] - q[0]) -
                 (q[0] - p[0]) * (r[1] - q[1]);

    if (mcm::kernels::FuzzyZero(val))
        return 0; // collinear

    return (val > 0.) ? 1 : 2; // clock or counterclock wise
}

// /** Check whether two Segment2D objects @a s1 and @a s2 intersect. */
// bool CheckIntersect(const int dim, const Segment2D &s1, const Segment2D &s2)
// {
//     // Find the four orientations needed for general and
//     // special cases
//     int o1 = Orientation(dim, s1.p[0], s1.p[1], s2.p[0]);
//     int o2 = Orientation(dim, s1.p[0], s1.p[1], s2.p[1]);
//     int o3 = Orientation(dim, s2.p[0], s2.p[1], s1.p[0]);
//     int o4 = Orientation(dim, s2.p[0], s2.p[1], s1.p[1]);

//     // General case
//     if (o1 != o2 && o3 != o4)
//         return true;

//     // Special Cases
//     // p1, q1 and p2 are collinear and p2 lies on segment p1q1
//     if (o1 == 0 && CollinearAndIntersect(dim, s1.p[0], s2.p[0], s1.p[1]))
//         return true;

//     // p1, q1 and q2 are collinear and q2 lies on segment p1q1
//     if (o2 == 0 && CollinearAndIntersect(dim, s1.p[0], s2.p[1], s1.p[1]))
//         return true;

//     // p2, q2 and p1 are collinear and p1 lies on segment p2q2
//     if (o3 == 0 && CollinearAndIntersect(dim, s2.p[0], s1.p[0], s2.p[1]))
//         return true;

//     // p2, q2 and q1 are collinear and q1 lies on segment p2q2
//     if (o4 == 0 && CollinearAndIntersect(dim, s2.p[0], s1.p[1], s2.p[1]))
//         return true;

//     return false; // Doesn't fall in any of the above cases
// }

// /** Find intersection @p of two Segment2D objects @a s1 and @a s2.
//  *  This method must be called after CheckIntersect(s1, s2). */
// void FindIntersection(const int dim, const Segment &s1, const Segment &s2, double *p)
// {
//     // Line AB represented as a1x + b1y = c1
//     double a1 = s1.p[1][1] - s1.p[0][1];
//     double b1 = s1.p[0][0] - s1.p[1][0];
//     double c1 = a1 * (s1.p[0][0]) + b1 * (s1.p[0][1]);

//     // Line CD represented as a2x + b2y = c2
//     double a2 = s2.p[1][1] - s2.p[0][1];
//     double b2 = s2.p[0][0] - s2.p[1][0];
//     double c2 = a2 * (s2.p[0][0]) + b2 * (s2.p[0][1]);

//     double determinant = a1 * b2 - a2 * b1;

//     p[0] = (b2 * c1 - b1 * c2) / determinant;
//     p[1] = (a1 * c2 - a2 * c1) / determinant;
// }