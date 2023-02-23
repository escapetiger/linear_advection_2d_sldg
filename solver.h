#ifndef SLDG_SOLVER_H
#define SLDG_SOLVER_H

#include "globals.h"
#include "kernels.h"

static constexpr int kMAX_N_DOF = 25;
static constexpr int kMAX_N_GL_1D = 6;
static constexpr int kMAX_N_GL_2D = 36;

/** Solve u_t + a*u_x + b*u_y = 0 */

typedef struct
{
    double xmin = 0.; ///< minimum value of x
    double xmax = 1.; ///< maximum value of x
    double ymin = 0.; ///< minimum value of y
    double ymax = 1.; ///< maximum value of y
    double tmax = 0.; ///< final simulation time
    double CFL = .8;  ///< CFL number
    double ax = 1.;   ///< wave speed in x
    double ay = 1.;   ///< wave speed in y
    int Nx;           ///< size of x-partition
    int Ny;           ///< size of y-partition
    int Nt;           ///< size of t-partition
    int order;        ///< approximation order of numerical scheme
} SolverOptions;

class Solver
{
public:
    // model setups
    double xmin = 0.; ///< minimum value of x
    double xmax = 1.; ///< maximum value of x
    double ymin = 0.; ///< minimum value of y
    double ymax = 1.; ///< maximum value of y
    double tmax = 0.; ///< final simulation time
    double CFL = .8;  ///< CFL number
    double ax = 1.;   ///< wave speed in x
    double ay = 1.;   ///< wave speed in y
    int Nx;           ///< size of x-partition
    int Ny;           ///< size of y-partition
    int Nt;           ///< size of t-partition
    int order;        ///< approximation order of numerical scheme

    // mesh information
    double *x_grid; ///< x-partition
    double *y_grid; ///< y-partition
    double hx;      ///< ...
    double hy;      ///< ...
    double ht;      ///< ...

    // DG information in the reference space
    int P;               ///< degrees of DG polynomial
    int N_dof_loc;       ///< number of local degrees of freedom
    int N_GL_1D;         ///< number of 1D Gauss-Legendre points
    int N_GL_2D;         ///< number of 2D Gauss-Legendre points
    double *x_GL_1D;     ///< 1D Gauss-Legendre points on (0,1)
    double *w_GL_1D;     ///< 1D Gauss-Legendre weights on (0,1)
    double *x_GL_2D;     ///< 2D Gauss-Legendre points on (0,1)x(0,1)
    double *y_GL_2D;     ///< 2D Gauss-Legendre points on (0,1)x(0,1)
    double *w_GL_2D;     ///< 2D Gauss-Legendre weights on (0,1)x(0,1)
    double **v_GL_2D;    ///< 2D basis values at x_GL_2D
    double **v_XL_GL_2D; ///< 2D basis values at (0, x_GL_1D)
    double **v_XR_GL_2D; ///< 2D basis values at (1, x_GL_1D)
    double **v_GL_YL_2D; ///< 2D basis values at (x_GL_1D, 0)
    double **v_GL_YR_2D; ///< 2D basis values at (x_GL_1D, 1)
    double *v_XC_YC_2D;  ///< 2D basis values at (0.5, 0.5)
    double *v_XL_YL_2D;  ///< 2D basis values at (0.0, 0.0)
    double *v_XR_YL_2D;  ///< 2D basis values at (1.0, 0.0)
    double *v_XL_YR_2D;  ///< 2D basis values at (0.0, 1.0)
    double *v_XR_YR_2D;  ///< 2D basis values at (1.0, 1.0)
    double **mat_v_u_2D; ///< DG matrix: v*u

    // local data structures
    int N_node;    ///< number of nodes
    int N_edge;    ///< number of edges
    int N_elem;    ///< number of elements
    NodeE *node_e; ///< global Eulerian node buffer
    NodeU *node_u; ///< global Upstream node buffer
    EdgeE *edge_e; ///< global Eulerian edge buffer
    EdgeU *edge_u; ///< global Upstream edge buffer
    ElemE *elem_e; ///< global Eulerian element buffer
    ElemU *elem_u; ///< global Upstream element buffer

    // global algebraic elements: vectors and matrices
    int N_dof_glo; ///< number of global degrees freedom
    double *U_tn;  ///< global DoFs at time t_n
    double *U_tn1; ///< global DoFs at time t_{n+1}
    double *U_RHS; ///< RHS terms

    // residuals
    double e_L1;   ///< residual in L1 norm
    double e_L2;   ///< residual in L2 norm
    double e_Linf; ///< residual in Linf norm

    Solver() = default;
    ~Solver() = default;

    void ReadOptions(int argc, char const *argv[]);
    void ReadOptions(SolverOptions *opt);
    void ReadMesh();
    void ReadDG();
    void Simulate();
    void CalcResidual();
    void GetResidual(const int nerr, double *e);
    void OutputResidual(const char *file);
    void OutputSolution2D(const char *file);
    void Destroy();

private:
    double InitData(double x, double y);
    double ExactData(double t, double x, double y);
    void Project();
    int SetTimeStep(double t);
    void TrackBack();
    void Assemble();
    void Clipping();
    void ClipElemU(ElemU *elem);
    void ClipElemE(ElemE *elem, ElemU *eu);
    void FinalizeClip(int par_e, ElemE *elem, ElemU *eu);
    void SetQuadRule();
    void MakeSubQuadRuleE(ElemE *elem, SubElem *sub);
    void MakeSubQuadRuleU(ElemU *elem, SubElem *sub);
    double CalcElemEPkBasis(const int k, ElemE *elem, double *x);
    double CalcElemUPkBasis(const int k, ElemU *elem, SubElem *sub, double *x);
    void Step();
    void SetZero(int size, double *v);
    void CopyL2G(int size, double *v);
    void CopyG2L(int size, const double *v);
    void CopyG2G(int size, const double *v1, double *v2);
};

#endif