# SLDG_LINEAR_ADVECTION_2D

DESCRIPTION : C++ implementation of SLDG method for 2D linear advection equation.

VERSION     : 1.0

REFERENCES  : 
[1] X. Cai, W. Guo, J.-M. Qiu, A High Order Conservative Semi-Lagrangian Discontinuous Galerkin Method for Two-Dimensional Transport Simulations, J Sci Comput. 73 (2017) 514–542. https://doi.org/10.1007/s10915-017-0554-0.

---

2023/2/24 

Error order tests for order = 1, 2, 3, 4 are passed. When order >= 5, "stack smashing detected".

BUGs :
1. When ElemU overlaps ElemE, the program may fail to capture the correct result.
2. The clipping algorithm is too complicate, due to complex local data structures, such as POI, Segment, SubElem.

IDEAs for optimizing algorithm:
1. We shall simplify the local data structures. For uniform rectangular mesh, SubElem can be determined only by (x0,y0), (x1,y1) and par.
