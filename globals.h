#ifndef SLDG_GLOBALS_H
#define SLDG_GLOBALS_H

#include "mcm/mcm.h"

using mcm::CommandLineReader;
using mcm::gerr;
using mcm::gout;
using mcm::gtime;
using mcm::mesh::Graph;

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
    int par;          ///< index of Eulerian element
    Vec2d p0, p1;     ///< p0: left bottom; p1: right top
    QuadPoint qp[25]; ///< quadrature points
} SubElem;

// eulerian element
typedef struct
{
    int nvert;          ///< number of vertices
    int nedge;          ///< number of edges
    int ndof;           ///< number of DoFs
    int nsub;           ///< number of SubElemEs
    int vert_buf[4];    ///< local vertex buffer, storing NodeE's index
    int edge_buf[4];    ///< local edge buffer, storing EdgeE's index
    SubElem sub_buf[4]; ///< local SubElem buffer
    double *udof;       ///< trial DoFs
    Vec2d p0;           ///< left bottom
    Vec2d p1;           ///< right top
    Vec2d poi;          ///< local position of inner POI
} ElemE;

// upstream element
typedef struct
{
    int nvert;          ///< number of vertices
    int nedge;          ///< number of edges
    int nsub;           ///< number of SubElemUs
    int vert_buf[4];    ///< vertex buffer, storing NodeU's index
    int edge_buf[4];    ///< edge buffer, storing EdgeU's index
    SubElem sub_buf[4]; ///< local SubElem buffer
    Vec2d p0;           ///< left bottom
    Vec2d p1;           ///< right top
    Vec2d poi;          ///< local position of inner POI
} ElemU;

#endif