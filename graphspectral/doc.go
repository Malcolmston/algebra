// Package graphspectral implements spectral graph theory using only the Go
// standard library.
//
// The package centres on two data types. [Matrix] is a small dense real matrix
// with the linear-algebra operations required for spectral analysis: addition,
// multiplication, transposition, determinants, linear solves, inverses,
// symmetric eigendecomposition (via the cyclic Jacobi method) and the leading
// eigenpair by power iteration. [Graph] is a finite undirected, optionally
// weighted, graph stored as a symmetric adjacency matrix; it produces the
// adjacency, degree, (combinatorial) Laplacian, normalized Laplacian,
// random-walk Laplacian, signless Laplacian, transition and incidence matrices,
// and answers structural questions such as connectivity, bipartiteness,
// regularity and density.
//
// On top of these it provides the classical spectral quantities: the adjacency
// and Laplacian spectra, the algebraic connectivity and the Fiedler vector,
// spectral bisection and k-way spectral clustering, the adjacency and Laplacian
// spectral gaps, the graph energy and Laplacian energy, the Estrada index, the
// number of spanning trees via the Matrix-Tree theorem, the Kirchhoff index and
// effective resistances through the Laplacian pseudoinverse, and the two Cheeger
// bounds relating the algebraic connectivity to the isoperimetric (Cheeger)
// constant.
//
// A family of centrality measures is included: degree, eigenvector, Katz,
// closeness and harmonic centralities, together with PageRank (both uniform and
// personalized) computed by power iteration on the Google matrix. Standard
// convenient graph constructors (path, cycle, complete, star, wheel, complete
// bipartite and empty graphs) are provided, largely because their spectra are
// known in closed form and make excellent test fixtures.
//
// Every routine is deterministic and depends only on packages math, sort,
// errors, fmt and strings. Vectors are plain []float64 slices; matrices are
// dense and are never assumed sparse.
package graphspectral
