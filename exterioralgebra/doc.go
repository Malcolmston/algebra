// Package exterioralgebra implements the exterior (Grassmann) algebra Λ(Rⁿ)
// together with the calculus of differential forms on Rⁿ.
//
// The central value type is [Form], an element of the exterior algebra over an
// n-dimensional real vector space. A Form is a (possibly mixed-grade)
// multivector: a real linear combination of basis blades e_{i₁}∧…∧e_{i_k},
// with 0 ≤ i₁ < … < i_k < n. Each basis blade is encoded compactly as a
// bitmask whose set bits are the participating indices, so the grade of a
// blade is the population count of its mask. Coefficients are stored sparsely,
// so only the nonzero terms occupy memory.
//
// On Forms the package provides:
//
//   - construction ([New], [Scalar], [Basis1], [BasisBlade], [FromVector],
//     [VolumeForm], [One]) and inspection ([Form.Coeff], [Form.Terms],
//     [Form.Grade], [Form.Grades], [Form.IsHomogeneous], [Form.String]);
//   - the graded-vector-space operations ([Form.Add], [Form.Sub], [Form.Neg],
//     [Form.Scale], [Form.GradeProject], [Form.ScalarPart]);
//   - the exterior (wedge) product ([Wedge], [Form.Wedge], [WedgeAll],
//     [Form.WedgePow]);
//   - the interior product / contraction ([InteriorProduct],
//     [Form.LeftContract], [Form.RightContract]) and the Euclidean inner
//     product and norm ([InnerProduct], [Form.Norm]);
//   - the Hodge star for the Euclidean and for a general diagonal metric
//     ([HodgeStar], [HodgeStarMetric], [InverseHodgeStar]), verified against
//     the de Rham identity ** = (−1)^{k(n−k)} id.
//
// For the differential calculus, [Poly] is a dense-free multivariate real
// polynomial and [PForm] is a differential form whose coefficients are
// polynomials. Because polynomial partial derivatives are computed exactly,
// the exterior derivative [PForm.ExteriorDerivative] satisfies d² = 0 exactly,
// the pullback [PForm.Pullback] along a polynomial map commutes with d, and
// the Hodge–de Rham Laplacian [PForm.HodgeLaplacian] reproduces the componentwise
// scalar Laplacian (up to the geometer's sign). A numerical exterior derivative
// ([NumExteriorDerivative]) is also provided for arbitrary smooth coefficient
// fields.
//
// All computation is deterministic and uses only the Go standard library.
// Operations that cannot express failure in their signature (for example
// [Form.Coeff]) return a neutral value on misuse; the others return one of the
// sentinel errors [ErrDim], [ErrIndex] or [ErrGrade].
package exterioralgebra
