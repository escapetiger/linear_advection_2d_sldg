#ifndef SLDG_GLOBALS_H
#define SLDG_GLOBALS_H

#include "mcm/mcm.h"

using mcm::CommandLineReader;
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

// eulerian node
typedef struct
{
    Vec2d pos; ///< position
} NodeE;

// upstream node
typedef struct
{
    int id;    ///< index of Eulerian element
    Vec2d pos; ///< position
} NodeU;

// intersection point
typedef struct
{
    Vec2d pos; // position
    int type;  // 0 : Grid-Grid; 1 : Grid-ElemU; 2 : ElemU-ElemU
    int id;    // index of Eulerian element
} POI;

// eulerian edge
typedef union
{
    struct
    {
        int beg; // begin index
        int end; // end index
    };
    int id[2]; // index of nodes in NodeE buffer
} EdgeE;

// upstream segment
typedef struct
{
    union
    {
        struct
        {
            POI beg;
            POI end;
        };
        POI poi[2];
    };
    double c_ab[6];
} SegmentU;

// upstream edge
typedef struct
{
    union
    {
        struct
        {
            int beg; // begin index
            int end; // end index
        };
        int id[2]; // index of nodes in NodeU buffer
    };
    int nsego;
    SegmentU sego_buf[2];
} EdgeU;

// eulerian element
typedef struct
{
    int nvert;       ///< number of vertices
    int nedge;       ///< number of edges
    int ndof;        ///< number of DoFs
    int vert_buf[9]; ///< local vertex buffer, storing NodeE's index
    int edge_buf[4]; ///< local edge buffer, storing EdgeE's index
    double *udof;    ///< trial DoFs
} ElemE;

// upstream element
typedef struct
{
    int nvert;            ///< number of vertices
    int nedge;            ///< number of edges
    int nsegi;            ///< number of inner segments
    int vert_buf[9];      ///< vertex buffer, storing NodeU's index
    int edge_buf[4];      ///< edge buffer, storing EdgeU's index
    SegmentU segi_buf[8]; ///< local inner segments buffer
    double **vcoe;        ///< coefficients of reconstructed test basis 
} ElemU;

#endif