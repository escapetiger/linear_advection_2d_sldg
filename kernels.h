#ifndef SLDG_KERNELS_H
#define SLDG_KERNELS_H

#include "globals.h"

/** Quickly set quadrature points. */
void SetQuadPoint(const int i, const double a, const double v,
                  double *x, double *w);

/** Get 1D Gauss-Legendre quadrature points. */
void GetGaussLegendre(const int n, double *x, double *w);

/** Get multi-dimensional Gauss-Legendre quadrature points. */
void GetGaussLegendre(const int dim, const int n, double *x, double *y, double *w);

/** Calculate (derivative) values of 1D Legendre polynomials. */
double Legendre1D(const int n, const int m, const double x);

/** Calculate values of multi-dimensional Qk basis. */
double CalcQkBasis(const int dim, const int k, const int *n, const double *x);

/** Calculate derivative values of multi-dimensional Legendre polynomials. */
double CalcQkBasis(const int dim, const int k, const int j, const int *n, const double *x);

/** Calculate value of multi-dimensional Pk basis. */
double CalcPkBasis(const int dim, const int k, const int n, const double *x);

/** Calculate weighted summation. */
double CalcWeightedSum(const int n, const double *w, const double *x);

/** Circulate `double` number in [xL,xR]. */
double Circulate(double x, double xL, double xR);

/** Circulate `int` number in [xL,xR). */
int Circulate(int x, int xL, int xR);

/** The reconstruction problem is equivalent to solve Ax = b.
 * P1-Reconstion:
 *   A  = [1 x1 y1;
 *         1 x2 y2;
 *         1 x3 y3;
 *         1 x4 y4;];
 *   b  = [v1 v2 v3 v4]';
 *   A2 = [1  sx  sy;
 *        sx sx2 sxy;
 *        sy sxy sy2];
 *   b2 = [sv1 sv2 sv3 sv4]';
 * Q1-Reconstion (reduced system!!):
 *   A = [1 x1 y1 x1*y1;
 *        1 x2 y2 x2*y2;
 *        1 x3 y3 x3*y3;
 *        1 x4 y4 x4*y4];
 *   b = [v1 v2 v3 v4];
 *   x = [1 0 0 0] or [0 1 0 0] or [0 0 1 0] or [0 0 0 1];
 *
 * SUMMARY:
 * For Pk-reconstruction, a least square problem shall be solved.
 * For Qk-reconstruction, we need redefine the problem.
 */
void PkReconst(const int k, const double **Ap, const double *bp, double *xp);


#endif