# SLDG_LINEAR_ADVECTION_2D

DESCRIPTION : C++ implementation of SLDG method for 2D linear advection equation.

VERSION     : 3.0

REFERENCES  : 
[1] X. Cai, W. Guo, J.-M. Qiu, A High Order Conservative Semi-Lagrangian Discontinuous Galerkin Method for Two-Dimensional Transport Simulations, J Sci Comput. 73 (2017) 514–542. https://doi.org/10.1007/s10915-017-0554-0.

---

## update at 2023/2/25

Universal fast implementation for uniform mesh, 10 times faster than VERSION 2.0.

**Benchmark:** 

Error order tests for order = 1, 2, 3, 4, 5 are passed. 
Only periodic boundary conditions are considered.

**Future work:**

Extend this idea to 3D case.

