// Package tilings provides tools for planar symmetry and tilings: the 17
// wallpaper groups and 7 frieze groups (their generators, point groups and
// orbits), planar isometries (rotations, reflections, translations and glide
// reflections) together with the finite point groups they build, the regular
// and semiregular (Archimedean) tilings described by vertex configurations,
// substitution tilings (Penrose via Robinson triangles, the chair/L-tromino
// rep-tile and the pinwheel triangle), and two-dimensional lattice and
// orbifold utilities.
//
// The geometric core is the [Isometry] type, an affine map x -> M x + t whose
// linear part M is orthogonal (determinant +/-1). Every planar isometry is one
// of five kinds - identity, translation, rotation, reflection or glide
// reflection - and [Isometry.Classify] recovers the kind together with the
// invariant data (rotation centre and angle, mirror axis, glide vector). Finite
// collections of origin-fixing isometries form the cyclic groups [CyclicGroup]
// and dihedral groups [DihedralGroup] used as planar point groups.
//
// Wallpaper and frieze groups are represented by their translation lattice and
// a finite set of coset representatives (the quotient by the translation
// subgroup, which is the point group). The representatives are obtained by
// closing a small set of generators under composition, reducing translation
// parts modulo the lattice; [WallpaperGroup.Orbit] and [FriezeGroup.Orbit] then
// expand a seed point into a finite patch of its symmetry orbit.
//
// Regular and Archimedean tilings are handled combinatorially through
// [VertexConfiguration], the cyclic sequence of polygon sizes meeting at a
// vertex; the interior-angle sum of a legal planar vertex is exactly 2*pi.
// Substitution tilings are modelled as sets of tiles that a substitution rule
// deflates into smaller similar tiles: [DeflateRobinson] realises the Penrose
// P2 and P3 tilings, [Chair] realises the chair (L-tromino) rep-tile, and
// [PinwheelTriangle] realises Conway and Radin's pinwheel.
//
// All computation is performed in float64 using only the Go standard library
// and is deterministic. Any randomness is driven by a caller-supplied
// *math/rand.Rand so that results are reproducible. Predicates that must
// tolerate floating-point noise take or use an explicit epsilon; the package
// default is [Epsilon].
package tilings
