#ifndef SLDG_GLOBALS_H
#define SLDG_GLOBALS_H

#include "mcm/mcm.h"

using mcm::CommandLineReader;
using mcm::mesh::Graph;
using mcm::gerr;
using mcm::gout;
using mcm::gtime;

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

// intersection point
typedef struct
{
    Vec2d pos; // position
    int type;  // 0 : Edge-Edge; 1 : Grid-Edge; 2 : Grid-Grid;
    int par;   // index of Eulerian element
} POI;

// segment consists of POIs
typedef struct
{
    union
    {
        struct
        {
            int beg; // begin index of POI
            int end; // end index of POI
        };
        int poi[2]; // index of POIs in local buffer
    };
} Segment;

// sub-element consists of Segments
typedef struct
{
    int par;               ///< index of Eulerian element
    int nseg;              ///< number of Segments
    int seg[4];            ///< indices of Segments
    double x0, y0, x1, y1; ///< fast check bound
    QuadPoint qp[4];       ///< quadrature points
} SubElem;

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
    int nvert;             ///< number of vertices
    int nedge;             ///< number of edges
    int ndof;              ///< number of DoFs
    int npoi;              ///< number of POIs
    int nseg;              ///< number of segments
    int nsub;              ///< number of SubElemEs
    int vert_buf[4];       ///< local vertex buffer, storing NodeE's index
    int edge_buf[4];       ///< local edge buffer, storing EdgeE's index
    POI poi_buf[32];       ///< local POIs buffer
    Segment seg_buf[16];   ///< local Segment buffer
    SubElem sub_buf[4];    ///< local SubElemEs buffer
    double *udof;          ///< trial DoFs
    double x0, y0, x1, y1; ///< fast check bound
} ElemE;

// upstream element
typedef struct
{
    int nvert;           ///< number of vertices
    int nedge;           ///< number of edges
    int npoi;            ///< number of POIs
    int nseg;            ///< number of segments
    int nsub;            ///< number of SubElemUs
    int vert_buf[4];     ///< vertex buffer, storing NodeU's index
    int edge_buf[4];     ///< edge buffer, storing EdgeU's index
    POI poi_buf[32];     ///< local POIs buffer
    Segment seg_buf[16]; ///< local Segment buffer
    SubElem sub_buf[4];  ///< local SubElemUs buffer
} ElemU;

#endif