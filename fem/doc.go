// Package fem implements the finite element method (FEM) for elliptic
// boundary value problems in one and two spatial dimensions.
//
// The package is written entirely with the Go standard library and provides a
// self-contained toolkit that covers the full FEM pipeline:
//
//   - Dense and sparse linear algebra (LU factorisation, conjugate gradients).
//   - Numerical quadrature on intervals (Gauss–Legendre) and on triangles
//     (symmetric Dunavant/Strang rules up to degree five).
//   - Lagrange P1 and P2 basis functions on the reference interval and the
//     reference triangle, together with their gradients.
//   - Interval and triangular mesh generation, boundary extraction and uniform
//     mesh refinement.
//   - Element stiffness, mass and load matrices and their assembly into global
//     sparse operators.
//   - Application of Dirichlet, Neumann and Robin boundary conditions.
//   - Ready-made solvers for the Poisson equation, reaction–diffusion problems
//     and two–dimensional linear elasticity.
//   - L2 and H1 (semi)norms and the corresponding error estimators.
//
// Barycentric coordinates on a triangle with vertices v1, v2, v3 are ordered
// (L1, L2, L3) with Li = 1 at vi. P2 nodes are ordered as the three vertices
// followed by the midpoints of edges (v2,v3), (v1,v3) and (v1,v2).
//
// Randomness, where required, always uses a caller supplied *math/rand.Rand so
// that results are reproducible.
package fem
