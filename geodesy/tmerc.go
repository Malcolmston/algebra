package geodesy

import "math"

// TransverseMercatorForward projects geodetic latitude/longitude (degrees) to a
// transverse Mercator grid on the given ellipsoid, using the Krüger n-series
// (accurate to sub-millimetre within a few degrees of the central meridian).
// lon0 is the central meridian (degrees), k0 the central scale factor, and
// falseEasting/falseNorthing the grid origin offsets (metres). It returns the
// easting and northing in metres.
func TransverseMercatorForward(lat, lon, lon0, k0, falseEasting, falseNorthing float64, e Ellipsoid) (easting, northing float64) {
	φ := rad(lat)
	λ := rad(lon - lon0)

	ecc := e.FirstEccentricity()
	n := e.ThirdFlattening()
	n2 := n * n
	n3 := n2 * n
	n4 := n3 * n
	n5 := n4 * n
	n6 := n5 * n

	cosλ := math.Cos(λ)
	sinλ := math.Sin(λ)

	τ := math.Tan(φ)
	σ := math.Sinh(ecc * math.Atanh(ecc*τ/math.Sqrt(1+τ*τ)))
	τʹ := τ*math.Sqrt(1+σ*σ) - σ*math.Sqrt(1+τ*τ)

	ξʹ := math.Atan2(τʹ, cosλ)
	ηʹ := math.Asinh(sinλ / math.Sqrt(τʹ*τʹ+cosλ*cosλ))

	A := e.A / (1 + n) * (1 + n2/4 + n4/64 + n6/256)

	α := [7]float64{
		0,
		1.0/2*n - 2.0/3*n2 + 5.0/16*n3 + 41.0/180*n4 - 127.0/288*n5 + 7891.0/37800*n6,
		13.0/48*n2 - 3.0/5*n3 + 557.0/1440*n4 + 281.0/630*n5 - 1983433.0/1935360*n6,
		61.0/240*n3 - 103.0/140*n4 + 15061.0/26880*n5 + 167603.0/181440*n6,
		49561.0/161280*n4 - 179.0/168*n5 + 6601661.0/7257600*n6,
		34729.0/80640*n5 - 3418889.0/1995840*n6,
		212378941.0 / 319334400 * n6,
	}

	ξ := ξʹ
	η := ηʹ
	for j := 1; j <= 6; j++ {
		ξ += α[j] * math.Sin(2*float64(j)*ξʹ) * math.Cosh(2*float64(j)*ηʹ)
		η += α[j] * math.Cos(2*float64(j)*ξʹ) * math.Sinh(2*float64(j)*ηʹ)
	}

	easting = k0*A*η + falseEasting
	northing = k0*A*ξ + falseNorthing
	return easting, northing
}

// TransverseMercatorInverse recovers geodetic latitude/longitude (degrees) from
// a transverse Mercator grid coordinate (easting, northing in metres) on the
// given ellipsoid, inverting TransverseMercatorForward with the same
// parameters.
func TransverseMercatorInverse(easting, northing, lon0, k0, falseEasting, falseNorthing float64, e Ellipsoid) (lat, lon float64) {
	ecc := e.FirstEccentricity()
	n := e.ThirdFlattening()
	n2 := n * n
	n3 := n2 * n
	n4 := n3 * n
	n5 := n4 * n
	n6 := n5 * n

	A := e.A / (1 + n) * (1 + n2/4 + n4/64 + n6/256)

	η := (easting - falseEasting) / (k0 * A)
	ξ := (northing - falseNorthing) / (k0 * A)

	β := [7]float64{
		0,
		1.0/2*n - 2.0/3*n2 + 37.0/96*n3 - 1.0/360*n4 - 81.0/512*n5 + 96199.0/604800*n6,
		1.0/48*n2 + 1.0/15*n3 - 437.0/1440*n4 + 46.0/105*n5 - 1118711.0/3870720*n6,
		17.0/480*n3 - 37.0/840*n4 - 209.0/4480*n5 + 5569.0/90720*n6,
		4397.0/161280*n4 - 11.0/504*n5 - 830251.0/7257600*n6,
		4583.0/161280*n5 - 108847.0/3991680*n6,
		20648693.0 / 638668800 * n6,
	}

	ξʹ := ξ
	ηʹ := η
	for j := 1; j <= 6; j++ {
		ξʹ -= β[j] * math.Sin(2*float64(j)*ξ) * math.Cosh(2*float64(j)*η)
		ηʹ -= β[j] * math.Cos(2*float64(j)*ξ) * math.Sinh(2*float64(j)*η)
	}

	sinhηʹ := math.Sinh(ηʹ)
	sinξʹ := math.Sin(ξʹ)
	cosξʹ := math.Cos(ξʹ)
	τʹ := sinξʹ / math.Sqrt(sinhηʹ*sinhηʹ+cosξʹ*cosξʹ)

	// Solve for τ from τʹ by Newton's method (Karney).
	τ := τʹ
	for i := 0; i < 100; i++ {
		σ := math.Sinh(ecc * math.Atanh(ecc*τ/math.Sqrt(1+τ*τ)))
		τi := τ*math.Sqrt(1+σ*σ) - σ*math.Sqrt(1+τ*τ)
		dτ := (τʹ - τi) / math.Sqrt(1+τi*τi) *
			(1 + (1-ecc*ecc)*τ*τ) / ((1 - ecc*ecc) * math.Sqrt(1+τ*τ))
		τ += dτ
		if math.Abs(dτ) < 1e-14 {
			break
		}
	}

	φ := math.Atan(τ)
	λ := math.Atan2(sinhηʹ, cosξʹ)
	return deg(φ), lon0 + deg(λ)
}
