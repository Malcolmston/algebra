package physics

import "math"

// physicsDefaultStep is the finite-difference step used by the numerical
// vector-calculus operators when the caller passes a non-positive step. It is
// small enough for second-order central differences to be accurate yet large
// enough to stay clear of float64 cancellation for typical inputs.
const physicsDefaultStep = 1e-6

// physicsStep returns h when it is positive and [physicsDefaultStep] otherwise.
// Centralising the fallback keeps every operator's step handling identical and
// deterministic.
func physicsStep(h float64) float64 {
	if h > 0 {
		return h
	}
	return physicsDefaultStep
}

// Unit returns the unit vector (length 1) pointing in the direction of a. If a
// is the zero vector its direction is undefined and the zero vector is returned
// unchanged. Unit is the allocation-free companion of [Vec3.Normalize] and, like
// all Vec3 arithmetic, returns its result by value without touching the heap.
func (a Vec3) Unit() Vec3 {
	n := a.Norm()
	if n == 0 {
		return a
	}
	inv := 1 / n
	return Vec3{a.X * inv, a.Y * inv, a.Z * inv}
}

// AngleBetween returns the unsigned angle in radians between a and b, in the
// range [0, π]. If either operand is the zero vector the angle is undefined and
// 0 is returned. The cosine is clamped to [-1, 1] before inversion so that
// rounding cannot push [math.Acos] out of its domain.
func (a Vec3) AngleBetween(b Vec3) float64 {
	na, nb := a.Norm(), b.Norm()
	if na == 0 || nb == 0 {
		return 0
	}
	return math.Acos(physicsClampUnit(a.Dot(b) / (na * nb)))
}

// AddInto stores the component-wise sum a + b into *dst. It is the in-place form
// of [Vec3.Add]: hot loops that repeatedly accumulate into an existing vector
// can use it to make their intent explicit while performing, like the value
// form, no heap allocation. dst may alias a or b.
func AddInto(dst *Vec3, a, b Vec3) {
	dst.X = a.X + b.X
	dst.Y = a.Y + b.Y
	dst.Z = a.Z + b.Z
}

// ScaleInto stores a scaled by the scalar s into *dst. It is the in-place form
// of [Vec3.Scale] and, like [AddInto], allocates nothing. dst may alias a.
func ScaleInto(dst *Vec3, a Vec3, s float64) {
	dst.X = a.X * s
	dst.Y = a.Y * s
	dst.Z = a.Z * s
}

// SumSlice returns the component-wise sum of every vector in v. It makes a
// single pass over the slice accumulating into scalar registers, allocating no
// intermediate vectors, and returns the zero vector for an empty or nil slice.
func SumSlice(v []Vec3) Vec3 {
	var sx, sy, sz float64
	for i := range v {
		sx += v[i].X
		sy += v[i].Y
		sz += v[i].Z
	}
	return Vec3{sx, sy, sz}
}

// Gradient returns a numerical approximation of the gradient ∇f of the scalar
// field f at the point p, using second-order central finite differences with
// step h. Each component is (f(p+h·eᵢ) − f(p−h·eᵢ)) / (2h) along the
// corresponding axis. If h ≤ 0 a sensible default step is used. Central
// differences make the result deterministic and second-order accurate.
func Gradient(f func(Vec3) float64, p Vec3, h float64) Vec3 {
	h = physicsStep(h)
	inv := 1 / (2 * h)
	dx := (f(Vec3{p.X + h, p.Y, p.Z}) - f(Vec3{p.X - h, p.Y, p.Z})) * inv
	dy := (f(Vec3{p.X, p.Y + h, p.Z}) - f(Vec3{p.X, p.Y - h, p.Z})) * inv
	dz := (f(Vec3{p.X, p.Y, p.Z + h}) - f(Vec3{p.X, p.Y, p.Z - h})) * inv
	return Vec3{dx, dy, dz}
}

// Divergence returns a numerical approximation of the divergence ∇·f of the
// vector field f at the point p, using second-order central finite differences
// with step h. It sums the diagonal partial derivatives
// ∂fₓ/∂x + ∂f_y/∂y + ∂f_z/∂z. If h ≤ 0 a sensible default step is used.
func Divergence(f func(Vec3) Vec3, p Vec3, h float64) float64 {
	h = physicsStep(h)
	inv := 1 / (2 * h)
	dfxdx := (f(Vec3{p.X + h, p.Y, p.Z}).X - f(Vec3{p.X - h, p.Y, p.Z}).X) * inv
	dfydy := (f(Vec3{p.X, p.Y + h, p.Z}).Y - f(Vec3{p.X, p.Y - h, p.Z}).Y) * inv
	dfzdz := (f(Vec3{p.X, p.Y, p.Z + h}).Z - f(Vec3{p.X, p.Y, p.Z - h}).Z) * inv
	return dfxdx + dfydy + dfzdz
}

// Curl returns a numerical approximation of the curl ∇×f of the vector field f
// at the point p, using second-order central finite differences with step h.
// The components are (∂f_z/∂y − ∂f_y/∂z, ∂fₓ/∂z − ∂f_z/∂x, ∂f_y/∂x − ∂fₓ/∂y).
// If h ≤ 0 a sensible default step is used.
func Curl(f func(Vec3) Vec3, p Vec3, h float64) Vec3 {
	h = physicsStep(h)
	inv := 1 / (2 * h)

	fxp := f(Vec3{p.X + h, p.Y, p.Z})
	fxm := f(Vec3{p.X - h, p.Y, p.Z})
	fyp := f(Vec3{p.X, p.Y + h, p.Z})
	fym := f(Vec3{p.X, p.Y - h, p.Z})
	fzp := f(Vec3{p.X, p.Y, p.Z + h})
	fzm := f(Vec3{p.X, p.Y, p.Z - h})

	dfzdy := (fyp.Z - fym.Z) * inv
	dfydz := (fzp.Y - fzm.Y) * inv
	dfxdz := (fzp.X - fzm.X) * inv
	dfzdx := (fxp.Z - fxm.Z) * inv
	dfydx := (fxp.Y - fxm.Y) * inv
	dfxdy := (fyp.X - fym.X) * inv

	return Vec3{
		dfzdy - dfydz,
		dfxdz - dfzdx,
		dfydx - dfxdy,
	}
}

// Laplacian returns a numerical approximation of the Laplacian ∇²f of the scalar
// field f at the point p, using the second-order central second-difference
// stencil with step h. It sums the unmixed second partial derivatives
// ∂²f/∂x² + ∂²f/∂y² + ∂²f/∂z², each evaluated as
// (f(p+h·eᵢ) − 2f(p) + f(p−h·eᵢ)) / h². If h ≤ 0 a sensible default step is used.
func Laplacian(f func(Vec3) float64, p Vec3, h float64) float64 {
	h = physicsStep(h)
	inv := 1 / (h * h)
	c := f(p)
	dxx := (f(Vec3{p.X + h, p.Y, p.Z}) - 2*c + f(Vec3{p.X - h, p.Y, p.Z})) * inv
	dyy := (f(Vec3{p.X, p.Y + h, p.Z}) - 2*c + f(Vec3{p.X, p.Y - h, p.Z})) * inv
	dzz := (f(Vec3{p.X, p.Y, p.Z + h}) - 2*c + f(Vec3{p.X, p.Y, p.Z - h})) * inv
	return dxx + dyy + dzz
}
