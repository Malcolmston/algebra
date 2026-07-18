package fractal

// Viewport is an axis-aligned rectangle in the complex plane, used to map
// integer pixel coordinates onto complex values. XMin/XMax bound the real axis
// and YMin/YMax bound the imaginary axis.
type Viewport struct {
	XMin, XMax, YMin, YMax float64
}

// NewViewport returns the square viewport centered at (centerRe, centerIm) with
// the given half-width (radius) along each axis.
func NewViewport(centerRe, centerIm, radius float64) Viewport {
	return Viewport{
		XMin: centerRe - radius,
		XMax: centerRe + radius,
		YMin: centerIm - radius,
		YMax: centerIm + radius,
	}
}

// SpanX returns the width of the viewport along the real axis.
func (v Viewport) SpanX() float64 { return v.XMax - v.XMin }

// SpanY returns the height of the viewport along the imaginary axis.
func (v Viewport) SpanY() float64 { return v.YMax - v.YMin }

// Center returns the complex number at the center of the viewport.
func (v Viewport) Center() complex128 {
	return complex((v.XMin+v.XMax)/2, (v.YMin+v.YMax)/2)
}

// Aspect returns the width-to-height ratio of the viewport. It returns 0 when
// the height is zero.
func (v Viewport) Aspect() float64 {
	h := v.SpanY()
	if h == 0 {
		return 0
	}
	return v.SpanX() / h
}

// PixelToComplex maps the pixel at column px, row py of a width×height raster to
// the complex value it samples. Column 0 maps to XMin and column width-1 maps
// to XMax; row 0 maps to YMax (top) and row height-1 maps to YMin (bottom), so
// increasing the row index moves downward as is conventional for images. When
// width or height is 1 the corresponding axis collapses to its minimum edge.
func (v Viewport) PixelToComplex(px, py, width, height int) complex128 {
	var re, im float64
	if width <= 1 {
		re = v.XMin
	} else {
		re = v.XMin + v.SpanX()*float64(px)/float64(width-1)
	}
	if height <= 1 {
		im = v.YMax
	} else {
		im = v.YMax - v.SpanY()*float64(py)/float64(height-1)
	}
	return complex(re, im)
}

// Zoom returns a new viewport with the same center whose span along each axis is
// multiplied by factor. A factor below 1 zooms in; above 1 zooms out.
func (v Viewport) Zoom(factor float64) Viewport {
	cx := (v.XMin + v.XMax) / 2
	cy := (v.YMin + v.YMax) / 2
	hx := v.SpanX() / 2 * factor
	hy := v.SpanY() / 2 * factor
	return Viewport{cx - hx, cx + hx, cy - hy, cy + hy}
}
