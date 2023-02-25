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
static constexpr int MAX_SUB_ELEM = 4;
static constexpr int MAX_SUB_EDGE = 2;
static constexpr int MAX_VERT = 4;
static constexpr int MAX_EDGE = 4;
static constexpr int NUM_LIMIT = 2;

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

// sub-element
typedef struct
{
    Vec2d p0; ///< left bottom
    Vec2d p1; ///< right top
} SubElem;

// sub-edge
typedef struct
{
    Vec2d p0; ///< begin
    Vec2d p1; ///< end
} SubEdge;

// reference element
typedef struct ELEM_REF
{
    const Vec2d p0 = {0., 0.};                           ///< left bottom
    const Vec2d p1 = {1., 1.};                           ///< right top
    int nsub;                                            ///< number of SubElems
    SubElem sub_buf[MAX_SUB_ELEM];                       ///< local SubElem buffer
    Vec2d poi;                                           ///< local coordinate of inner POI, in REF coordinate
    QuadPoint qp_e_buf[MAX_SUB_ELEM][MAX_GL_2D];         ///< reference quadrature points for ElemE
    QuadPoint qp_u_buf[MAX_SUB_ELEM][MAX_GL_2D];         ///< reference quadrature points for ElemU
    double v_qp_e_buf[MAX_SUB_ELEM][MAX_DOF][MAX_GL_2D]; ///< precomputed basis values at qp_e
    double v_qp_u_buf[MAX_SUB_ELEM][MAX_DOF][MAX_GL_2D]; ///< precomputed basis values at qp_u
} ElemR;

// reference bottom edge
typedef struct EDGE_REF_B
{
    const Vec2d p0 = {0., 0.};                                      ///< left
    const Vec2d p1 = {1., 0.};                                      ///< right
    int nsub;                                                       ///< number of SubEdges
    SubEdge sub_buf[MAX_SUB_EDGE];                                  ///< local SubEdge buffer
    QuadPoint qp_e_buf[MAX_SUB_EDGE][MAX_GL_1D];                    ///< reference quadrature points for EdgeE
    QuadPoint qp_u_buf[MAX_SUB_EDGE][MAX_GL_1D];                    ///< reference quadrature points for EdgeU
    double v_qp_e_buf[MAX_SUB_EDGE][MAX_DOF][MAX_GL_1D][NUM_LIMIT]; ///< precomputed basis values at qp_e
    double v_qp_u_buf[MAX_SUB_EDGE][MAX_DOF][MAX_GL_1D][NUM_LIMIT]; ///< precomputed basis values at qp_u
} EdgeR_B;

// reference right edge
typedef struct EDGE_REF_R
{
    const Vec2d p0 = {1., 0.};                                      ///< bottom
    const Vec2d p1 = {1., 1.};                                      ///< top
    int nsub;                                                       ///< number of SubEdges
    SubEdge sub_buf[MAX_SUB_EDGE];                                  ///< local SubEdge buffer
    QuadPoint qp_e_buf[MAX_SUB_EDGE][MAX_GL_1D];                    ///< reference quadrature points for EdgeE
    QuadPoint qp_u_buf[MAX_SUB_EDGE][MAX_GL_1D];                    ///< reference quadrature points for EdgeU
    double v_qp_e_buf[MAX_SUB_EDGE][MAX_DOF][MAX_GL_1D][NUM_LIMIT]; ///< precomputed basis values at qp_e
    double v_qp_u_buf[MAX_SUB_EDGE][MAX_DOF][MAX_GL_1D][NUM_LIMIT]; ///< precomputed basis values at qp_u
} EdgeR_R;

// reference top edge
typedef struct EDGE_REF_T
{
    const Vec2d p0 = {0., 1.};                                      ///< left
    const Vec2d p1 = {1., 1.};                                      ///< right
    int nsub;                                                       ///< number of SubEdges
    SubEdge sub_buf[MAX_SUB_EDGE];                                  ///< local SubEdge buffer
    QuadPoint qp_e_buf[MAX_SUB_EDGE][MAX_GL_1D];                    ///< reference quadrature points for EdgeE
    QuadPoint qp_u_buf[MAX_SUB_EDGE][MAX_GL_1D];                    ///< reference quadrature points for EdgeU
    double v_qp_e_buf[MAX_SUB_EDGE][MAX_DOF][MAX_GL_1D][NUM_LIMIT]; ///< precomputed basis values at qp_e
    double v_qp_u_buf[MAX_SUB_EDGE][MAX_DOF][MAX_GL_1D][NUM_LIMIT]; ///< precomputed basis values at qp_u
} EdgeR_T;

// reference edge
typedef struct EDGE_REF
{
    const Vec2d p0 = {0., 0.};                                      ///< bottom
    const Vec2d p1 = {0., 1.};                                      ///< top
    int nsub;                                                       ///< number of SubEdges
    SubEdge sub_buf[MAX_SUB_EDGE];                                  ///< local SubEdge buffer
    QuadPoint qp_e_buf[MAX_SUB_EDGE][MAX_GL_1D];                    ///< reference quadrature points for EdgeE
    QuadPoint qp_u_buf[MAX_SUB_EDGE][MAX_GL_1D];                    ///< reference quadrature points for EdgeU
    double v_qp_e_buf[MAX_SUB_EDGE][MAX_DOF][MAX_GL_1D][NUM_LIMIT]; ///< precomputed basis values at qp_e
    double v_qp_u_buf[MAX_SUB_EDGE][MAX_DOF][MAX_GL_1D][NUM_LIMIT]; ///< precomputed basis values at qp_u
} EdgeR_L;

// quadrature rule in ElemE
typedef struct
{
    QuadPoint qp_buf[MAX_SUB_ELEM][MAX_GL_2D]; ///< quadrature points, in parent ElemE's REF coordinate
} QuadRuleElemE;

// quadrature rule in ElemU
typedef struct
{
    int par_buf[MAX_SUB_ELEM];                 ///< parent Eulerian element index for each SubElem
    QuadPoint qp_buf[MAX_SUB_ELEM][MAX_GL_2D]; ///< quadrature points, in parent ElemE's REF coordinate
} QuadRuleElemU;

// quadrature rule in EdgeE
typedef struct
{
    QuadPoint qp_buf[MAX_SUB_EDGE][MAX_GL_1D]; ///< quadrature points, in parent EdgeE's REF coordinate
} QuadRuleEdgeE;

// quadrature rule in EdgeU
typedef struct
{
    int par_buf[MAX_SUB_EDGE];                 ///< parent Eulerian element index for each SubEdge
    QuadPoint qp_buf[MAX_SUB_EDGE][MAX_GL_1D]; ///< quadrature points, in parent EdgeU's REF coordinate
} QuadRuleEdgeU;

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
    QuadRuleElemE qr;       ///< quadrature rule
    QuadRuleEdgeE qr_b;     ///< quadrature rule at bottom EdgeE
    QuadRuleEdgeE qr_r;     ///< quadrature rule at right EdgeE
    QuadRuleEdgeE qr_t;     ///< quadrature rule at top EdgeE
    QuadRuleEdgeE qr_l;     ///< quadrature rule at left EdgeE
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
    QuadRuleElemU qr;       ///< quadrature rule
    QuadRuleEdgeU qr_b;     ///< quadrature rule at bottom EdgeU
    QuadRuleEdgeU qr_r;     ///< quadrature rule at right EdgeU
    QuadRuleEdgeU qr_t;     ///< quadrature rule at top EdgeU
    QuadRuleEdgeU qr_l;     ///< quadrature rule at left EdgeU
} ElemU;

#endif