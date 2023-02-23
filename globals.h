#ifndef SLDG_GLOBALS_H
#define SLDG_GLOBALS_H

#include "mcm/mcm.h"

using mcm::CommandLineReader;
using mcm::gerr;
using mcm::gout;
using mcm::gtime;
using mcm::mesh::Graph;

// constexpr
static constexpr int MAX_P = 5;
static constexpr int MAX_DOF = 25;
static constexpr int MAX_GL_1D = 5;
static constexpr int MAX_GL_2D = 25;
static constexpr int MAX_SUB = 4;
static constexpr int MAX_VERT = 4;
static constexpr int MAX_EDGE = 4;

// 2D `double` vector
typedef union
{
    struct
    {
        double x, y;
    };
    double xy[2];
} Vec2d;

// 2D `int` vector
typedef union
{
    struct
    {
        int i, j;
    };
    int ij[2];
} Vec2i;

typedef struct
{
    union
    {
        struct
        {
            double x;
            double y;
        };
        double xy[2];
    };
    double w;
} QuadPoint;

// eulerian node
typedef struct
{
    Vec2d pos; ///< position
} NodeE;

// upstream node
typedef struct
{
    int par;   ///< index of Eulerian element
    Vec2d pos; ///< position
} NodeU;

// eulerian edge
typedef union
{
    struct
    {
        int beg; // begin index
        int end; // end index
    };
    int node[2]; // index of nodes in NodeE buffer
} EdgeE;

// upstream edge
typedef union
{
    struct
    {
        int beg; // begin index
        int end; // end index
    };
    int node[2]; // index of nodes in NodeU buffer
} EdgeU;

// sub-element consists of Segments
typedef struct
{
    Vec2d p0; ///< left bottom
    Vec2d p1; ///< right top
} SubElem;

// reference element
typedef struct REF
{
    const Vec2d p0 = {0., 0.};                      ///< left bottom
    const Vec2d p1 = {1., 1.};                      ///< right top
    int nsub;                                       ///< number of SubElems
    SubElem sub_buf[MAX_SUB];                       ///< local SubElem buffer
    Vec2d poi;                                      ///< local coordinate of inner POI, in REF coordinate
    QuadPoint qp_e_buf[MAX_SUB][MAX_GL_2D];         ///< reference quadrature points for ElemE
    QuadPoint qp_u_buf[MAX_SUB][MAX_GL_2D];         ///< reference quadrature points for ElemU
    double v_qp_e_buf[MAX_SUB][MAX_DOF][MAX_GL_2D]; ///< precomputed basis values at quadrature points
    double v_qp_u_buf[MAX_SUB][MAX_DOF][MAX_GL_2D]; ///< precomputed basis values at quadrature points
} ElemR;

// quadrature rule in ElemE
typedef struct
{
    QuadPoint qp_buf[MAX_SUB][MAX_GL_2D]; ///< quadrature points, in parent ElemE's REF coordinate
} QuadRuleE;

// quadrature rule in ElemU
typedef struct
{
    int par_buf[MAX_SUB];                 ///< parent Eulerian element index for each SubElem
    QuadPoint qp_buf[MAX_SUB][MAX_GL_2D]; ///< quadrature points, in parent ElemE's REF coordinate
} QuadRuleU;

// eulerian element
typedef struct
{
    // standard DG members
    int nvert;              ///< number of vertices
    int nedge;              ///< number of edges
    int ndof;               ///< number of DoFs
    int vert_buf[MAX_VERT]; ///< local vertex buffer, storing NodeE's index
    int edge_buf[MAX_EDGE]; ///< local edge buffer, storing EdgeE's index
    double *udof;           ///< trial DoFs
    Vec2d p0;               ///< left bottom
    Vec2d p1;               ///< right top
    QuadRuleE qr;           ///< quadrature rule
} ElemE;

// upstream element
typedef struct
{
    // standard DG members
    int nvert;              ///< number of vertices
    int nedge;              ///< number of edges
    int vert_buf[MAX_VERT]; ///< vertex buffer, storing NodeU's index
    int edge_buf[MAX_EDGE]; ///< edge buffer, storing EdgeU's index
    Vec2d p0;               ///< left bottom
    Vec2d p1;               ///< right top
    QuadRuleU qr;           ///< quadrature rule
} ElemU;

#endif