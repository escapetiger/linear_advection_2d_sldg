#ifndef _GLOBALS_H_
#define _GLOBALS_H_

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstring>

/** Allocate the Vector data. */
template <typename T>
T *NewVector(size_t size)
{
    T *v = new T[size];
    return v;
}

/** Allocate the row-major Matrix data. */
template <typename T>
T **NewMatrix(size_t rows, size_t cols)
{
    T **m = new T *[rows];
    for (size_t i = 0; i < rows; i++)
    {
        m[i] = new T[cols];
    }
    return m;
}

/** Allocate the row-major Tensor data. */
template <typename T>
T ***NewTensor3D(size_t n1, size_t n2, size_t n3)
{
    T ***p = new T **[n1];
    for (size_t i = 0; i < n1; i++)
    {
        p[i] = new T *[n2];
        for (size_t j = 0; j < n2; j++)
        {
            p[i][j] = new T[n3];
        }
    }
    return p;
}

/** Delete the Vector data. */
template <typename T>
void DeleteVector(T *p, size_t size)
{
    delete[] p;
    p = nullptr;
}

/** Delete the row-major Matrix data. */
template <typename T>
void DeleteMatrix(T **p, size_t rows, size_t cols)
{
    for (size_t i = 0; i < rows; i++)
    {
        delete[] (p[i]);
    }
    delete []p;
    p = nullptr;
}

/** Delete the row-major Matrix data. */
template <typename T>
void DeleteTensor3D(T ***p, size_t n1, size_t n2, size_t n3)
{
    for (size_t i = 0; i < n1; i++)
    {
        for (size_t j = 0; j < n2; j++)
        {
            delete[] (p[i][j]);
        }
        delete[] (p[i]);
    }
    delete[] p;
    p = nullptr;
}

/** Print the Vector data. */
template <typename T>
void PrintVector(T *v, size_t size, const char *msg = "",
                 std::ostream &out = std::cout)
{
    out << std::fixed << std::setprecision(8);
    out << msg << "Vector<" << size << ">\n";
    for (size_t i = 0; i < size; i++)
    {
        out << std::setw(13) << v[i] << '\n';
    }
    out << '\n';
}

/** Print the column-major Matrix data. */
template <typename T>
void PrintMatrix(T **m, size_t rows, size_t cols, const char *msg = "",
                 std::ostream &out = std::cout)
{
    out << std::fixed << std::setprecision(8);
    out << msg << "Matrix<" << rows << ',' << cols << ">\n";
    for (auto i = 0u; i < rows; i++)
    {
        for (auto j = 0u; j < cols; j++)
        {
            out << std::setw(13) << m[i][j] << "\t";
        }
        out << '\n';
    }
}

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

// eulerian segment
typedef struct
{
    NodeE node[2];
} SegmentE;

// upstream node
typedef struct
{
    Vec2i id;  ///< index of Eulerian element
    Vec2d pos; ///< position
} NodeU;

// intersection point
typedef struct
{
    Vec2d pos;
    int ixy_type;
    // note that on grid, index is i+1/2, so idx2=2i+1
    int idx2;
    int idy2;
    int igrid;
    int id_xy;
} POI;

typedef struct
{
    POI porigin;
    POI pend;
    Vec2i id;
    double c_ab[6];
} SegmentU;

typedef struct
{
    NodeU porigin;
    NodeU pend;
    POI poi[11];
    int nsub_outer;
    SegmentU segment_outer[10];
} Face;

typedef struct
{
    NodeE *vertex1, *vertex2, *vertex3, *vertex4;
    NodeE *vertex5, *vertex6, *vertex7, *vertex8;
    NodeE *vertex9;
    double *umodal;
} CellE;

typedef struct
{
    NodeU *vertex1, *vertex2, *vertex3, *vertex4;
    NodeU *vertex5, *vertex6, *vertex7, *vertex8;
    NodeU *vertex9;
    Face *edge1_lr, *edge2_bt, *edge3_lr, *edge4_bt;
    int nsub_inner;
    SegmentU segment_inner[32];
} CellU;

#endif