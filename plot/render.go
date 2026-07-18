package plot

import (
	"bytes"
	"image"
	"image/color"
	"image/png"
	"math"
	"os"
)

// transform maps data-space coordinates to pixel coordinates within the plot
// rectangle. The pixel y axis points down, so larger data y maps to smaller
// pixel y.
type transform struct {
	xmin, xmax, ymin, ymax float64
	px0, py0, px1, py1     float64 // plot rect: (px0,py0) top-left, (px1,py1) bottom-right
}

// X maps a data x coordinate to a pixel x coordinate.
func (t transform) X(x float64) float64 {
	return t.px0 + (x-t.xmin)/(t.xmax-t.xmin)*(t.px1-t.px0)
}

// Y maps a data y coordinate to a pixel y coordinate.
func (t transform) Y(y float64) float64 {
	return t.py1 - (y-t.ymin)/(t.ymax-t.ymin)*(t.py1-t.py0)
}

// canvas wraps an RGBA image with clipped drawing primitives.
type canvas struct {
	img *image.RGBA
}

// RenderPNG renders the figure and returns the encoded PNG bytes. The output is
// deterministic for a given figure.
func (f *Figure) RenderPNG() ([]byte, error) {
	img := f.raster()
	var buf bytes.Buffer
	enc := png.Encoder{CompressionLevel: png.DefaultCompression}
	if err := enc.Encode(&buf, img); err != nil {
		return nil, err
	}
	return buf.Bytes(), nil
}

// SavePNG renders the figure and writes the PNG to the file at path.
func (f *Figure) SavePNG(path string) error {
	b, err := f.RenderPNG()
	if err != nil {
		return err
	}
	return os.WriteFile(path, b, 0o644)
}

// raster draws the whole figure into a fresh RGBA image.
func (f *Figure) raster() *image.RGBA {
	w, h := f.Width, f.Height
	img := image.NewRGBA(image.Rect(0, 0, w, h))
	c := &canvas{img: img}
	c.fillRect(0, 0, w, h, f.Background)

	a := f.ax
	tr := transform{
		px0: marginLeft, py0: marginTop,
		px1: float64(w - marginRight), py1: float64(h - marginBottom),
	}
	tr.xmin, tr.xmax, tr.ymin, tr.ymax = a.dataRange()

	xticks, xstep := ticksFor(tr.xmin, tr.xmax, 8)
	yticks, ystep := ticksFor(tr.ymin, tr.ymax, 6)

	// Grid.
	if a.grid {
		for _, xv := range xticks {
			px := int(math.Round(tr.X(xv)))
			c.vLine(px, int(tr.py0), int(tr.py1), LightGray)
		}
		for _, yv := range yticks {
			py := int(math.Round(tr.Y(yv)))
			c.hLine(int(tr.px0), int(tr.px1), py, LightGray)
		}
	}

	// Series (clipped to the plot rectangle).
	clip := image.Rect(int(tr.px0), int(tr.py0), int(tr.px1)+1, int(tr.py1)+1)
	for _, s := range a.series {
		c.drawSeries(s, tr, clip)
	}

	// Plot border.
	c.rectOutline(int(tr.px0), int(tr.py0), int(tr.px1), int(tr.py1), Black)

	// Ticks and their labels.
	for _, xv := range xticks {
		px := int(math.Round(tr.X(xv)))
		c.vLine(px, int(tr.py1), int(tr.py1)+4, Black)
		lbl := formatTick(xv, xstep)
		c.text(px-textWidth(lbl, 1)/2, int(tr.py1)+7, lbl, Black, 1)
	}
	for _, yv := range yticks {
		py := int(math.Round(tr.Y(yv)))
		c.hLine(int(tr.px0)-4, int(tr.px0), py, Black)
		lbl := formatTick(yv, ystep)
		c.text(int(tr.px0)-6-textWidth(lbl, 1), py-textHeight(1)/2, lbl, Black, 1)
	}

	// Titles and axis labels.
	if a.title != "" {
		c.text((w-textWidth(a.title, 2))/2, 8, a.title, Black, 2)
	}
	if a.xlabel != "" {
		c.text((w-textWidth(a.xlabel, 1))/2, h-14, a.xlabel, Black, 1)
	}
	if a.ylabel != "" {
		cx := 12
		cy := (int(tr.py0) + int(tr.py1) + textWidth(a.ylabel, 1)) / 2
		c.textVert(cx, cy, a.ylabel, Black, 1)
	}

	if a.legend {
		c.drawLegend(a, w, tr)
	}
	return img
}

// --- primitives ---------------------------------------------------------

// set writes a single opaque pixel, ignoring out-of-bounds coordinates.
func (c *canvas) set(x, y int, col color.RGBA) {
	if x < 0 || y < 0 || x >= c.img.Rect.Dx() || y >= c.img.Rect.Dy() {
		return
	}
	c.img.SetRGBA(x, y, col)
}

// blend writes a pixel alpha-composited over the existing pixel (source-over).
func (c *canvas) blend(x, y int, col color.RGBA) {
	if x < 0 || y < 0 || x >= c.img.Rect.Dx() || y >= c.img.Rect.Dy() {
		return
	}
	if col.A == 0xff {
		c.img.SetRGBA(x, y, col)
		return
	}
	dst := c.img.RGBAAt(x, y)
	a := float64(col.A) / 255
	mix := func(s, d uint8) uint8 {
		return uint8(math.Round(float64(s)*a + float64(d)*(1-a)))
	}
	c.img.SetRGBA(x, y, color.RGBA{
		R: mix(col.R, dst.R), G: mix(col.G, dst.G), B: mix(col.B, dst.B), A: 0xff,
	})
}

// fillRect fills the rectangle [x0,x1) x [y0,y1) with col.
func (c *canvas) fillRect(x0, y0, x1, y1 int, col color.RGBA) {
	for y := y0; y < y1; y++ {
		for x := x0; x < x1; x++ {
			c.set(x, y, col)
		}
	}
}

// hLine draws a horizontal line from x0 to x1 (inclusive) at y.
func (c *canvas) hLine(x0, x1, y int, col color.RGBA) {
	if x1 < x0 {
		x0, x1 = x1, x0
	}
	for x := x0; x <= x1; x++ {
		c.set(x, y, col)
	}
}

// vLine draws a vertical line from y0 to y1 (inclusive) at x.
func (c *canvas) vLine(x, y0, y1 int, col color.RGBA) {
	if y1 < y0 {
		y0, y1 = y1, y0
	}
	for y := y0; y <= y1; y++ {
		c.set(x, y, col)
	}
}

// rectOutline strokes the border of the rectangle with corners (x0,y0),(x1,y1).
func (c *canvas) rectOutline(x0, y0, x1, y1 int, col color.RGBA) {
	c.hLine(x0, x1, y0, col)
	c.hLine(x0, x1, y1, col)
	c.vLine(x0, y0, y1, col)
	c.vLine(x1, y0, y1, col)
}

// line draws a line between two integer endpoints using Bresenham's algorithm,
// with a square pen of the given width, clipped to clip.
func (c *canvas) line(x0, y0, x1, y1, width int, col color.RGBA, clip image.Rectangle) {
	if width < 1 {
		width = 1
	}
	dx := abs(x1 - x0)
	dy := -abs(y1 - y0)
	sx := 1
	if x0 > x1 {
		sx = -1
	}
	sy := 1
	if y0 > y1 {
		sy = -1
	}
	err := dx + dy
	for {
		c.pen(x0, y0, width, col, clip)
		if x0 == x1 && y0 == y1 {
			break
		}
		e2 := 2 * err
		if e2 >= dy {
			err += dy
			x0 += sx
		}
		if e2 <= dx {
			err += dx
			y0 += sy
		}
	}
}

// pen stamps a filled square of the given width centred on (x,y), clipped.
func (c *canvas) pen(x, y, width int, col color.RGBA, clip image.Rectangle) {
	if width <= 1 {
		if inClip(x, y, clip) {
			c.blend(x, y, col)
		}
		return
	}
	half := width / 2
	for oy := -half; oy <= half; oy++ {
		for ox := -half; ox <= half; ox++ {
			px, py := x+ox, y+oy
			if inClip(px, py, clip) {
				c.blend(px, py, col)
			}
		}
	}
}

// disk fills a solid circle of the given radius centred on (x,y), clipped.
func (c *canvas) disk(x, y, r int, col color.RGBA, clip image.Rectangle) {
	if r < 1 {
		r = 1
	}
	for oy := -r; oy <= r; oy++ {
		for ox := -r; ox <= r; ox++ {
			if ox*ox+oy*oy <= r*r {
				px, py := x+ox, y+oy
				if inClip(px, py, clip) {
					c.blend(px, py, col)
				}
			}
		}
	}
}

// text draws s at the top-left pixel (x,y) using the built-in bitmap font at
// the given integer scale.
func (c *canvas) text(x, y int, s string, col color.RGBA, scale int) {
	if scale < 1 {
		scale = 1
	}
	cx := x
	for _, r := range s {
		g, ok := glyph5x7[r]
		if ok {
			for col2 := 0; col2 < glyphW; col2++ {
				bits := g[col2]
				for row := 0; row < glyphH; row++ {
					if bits&(1<<uint(row)) != 0 {
						c.fillScaled(cx+col2*scale, y+row*scale, scale, col)
					}
				}
			}
		}
		cx += glyphAdv * scale
	}
}

// textVert draws s rotated 90 degrees counter-clockwise, reading upward, with
// the first character's baseline near (x, y). Used for the y-axis label.
func (c *canvas) textVert(x, y int, s string, col color.RGBA, scale int) {
	if scale < 1 {
		scale = 1
	}
	cy := y
	for _, r := range s {
		g, ok := glyph5x7[r]
		if ok {
			for col2 := 0; col2 < glyphW; col2++ {
				bits := g[col2]
				for row := 0; row < glyphH; row++ {
					if bits&(1<<uint(row)) != 0 {
						// Rotate: glyph column -> screen y (up), glyph row -> screen x.
						c.fillScaled(x+row*scale, cy-col2*scale, scale, col)
					}
				}
			}
		}
		cy -= glyphAdv * scale
	}
}

// fillScaled fills a scale x scale block at (x,y) for scaled font rendering.
func (c *canvas) fillScaled(x, y, scale int, col color.RGBA) {
	for dy := 0; dy < scale; dy++ {
		for dx := 0; dx < scale; dx++ {
			c.set(x+dx, y+dy, col)
		}
	}
}

// --- series drawing -----------------------------------------------------

// drawSeries dispatches to the drawing routine for the concrete series type.
func (c *canvas) drawSeries(s series, tr transform, clip image.Rectangle) {
	switch v := s.(type) {
	case *LineSeries:
		c.drawPolyline(v.Xs, v.Ys, v.Width, v.Color, tr, clip)
	case *ScatterSeries:
		r := v.Size
		for i := range v.Xs {
			if i >= len(v.Ys) || isBad(v.Xs[i]) || isBad(v.Ys[i]) {
				continue
			}
			c.disk(int(math.Round(tr.X(v.Xs[i]))), int(math.Round(tr.Y(v.Ys[i]))), r, v.Color, clip)
		}
	case *BarSeries:
		for i := range v.Xs {
			if i >= len(v.Heights) || isBad(v.Xs[i]) || isBad(v.Heights[i]) {
				continue
			}
			c.drawBar(v.Xs[i]-v.Width/2, v.Xs[i]+v.Width/2, v.Baseline, v.Baseline+v.Heights[i], v.Color, tr, clip)
		}
	case *HistSeries:
		for i := range v.Counts {
			c.drawBar(v.Edges[i], v.Edges[i+1], 0, v.Counts[i], v.Color, tr, clip)
		}
	case *StepSeries:
		c.drawStep(v.Xs, v.Ys, v.Width, v.Color, tr, clip)
	case *FillSeries:
		c.drawFill(v, tr, clip)
	}
}

// drawPolyline strokes the connected finite runs of the point list.
func (c *canvas) drawPolyline(xs, ys []float64, width int, col color.RGBA, tr transform, clip image.Rectangle) {
	n := len(xs)
	if len(ys) < n {
		n = len(ys)
	}
	prevOK := false
	var px, py int
	for i := 0; i < n; i++ {
		if isBad(xs[i]) || isBad(ys[i]) {
			prevOK = false
			continue
		}
		x := int(math.Round(tr.X(xs[i])))
		y := int(math.Round(tr.Y(ys[i])))
		if prevOK {
			c.line(px, py, x, y, width, col, clip)
		}
		px, py, prevOK = x, y, true
	}
}

// drawStep strokes a piecewise-constant (post-step) path through the points.
func (c *canvas) drawStep(xs, ys []float64, width int, col color.RGBA, tr transform, clip image.Rectangle) {
	n := len(xs)
	if len(ys) < n {
		n = len(ys)
	}
	prevOK := false
	var px, py int
	for i := 0; i < n; i++ {
		if isBad(xs[i]) || isBad(ys[i]) {
			prevOK = false
			continue
		}
		x := int(math.Round(tr.X(xs[i])))
		y := int(math.Round(tr.Y(ys[i])))
		if prevOK {
			c.line(px, py, x, py, width, col, clip) // horizontal segment
			c.line(x, py, x, y, width, col, clip)   // vertical riser
		}
		px, py, prevOK = x, y, true
	}
}

// drawBar fills and outlines a single bar spanning data x [xl,xr], y [yb,yt].
func (c *canvas) drawBar(xl, xr, yb, yt float64, col color.RGBA, tr transform, clip image.Rectangle) {
	x0 := int(math.Round(tr.X(xl)))
	x1 := int(math.Round(tr.X(xr)))
	y0 := int(math.Round(tr.Y(yt)))
	y1 := int(math.Round(tr.Y(yb)))
	if x1 < x0 {
		x0, x1 = x1, x0
	}
	if y1 < y0 {
		y0, y1 = y1, y0
	}
	for y := y0; y <= y1; y++ {
		for x := x0; x <= x1; x++ {
			if inClip(x, y, clip) {
				c.blend(x, y, col)
			}
		}
	}
	edge := color.RGBA{R: col.R / 2, G: col.G / 2, B: col.B / 2, A: 0xff}
	for x := x0; x <= x1; x++ {
		if inClip(x, y0, clip) {
			c.set(x, y0, edge)
		}
	}
	c.line(x0, y0, x0, y1, 1, edge, clip)
	c.line(x1, y0, x1, y1, 1, edge, clip)
}

// drawFill shades the vertical band between the two curves of a FillSeries by
// filling each x column between y1 and y2.
func (c *canvas) drawFill(v *FillSeries, tr transform, clip image.Rectangle) {
	for i := 0; i+1 < len(v.Xs); i++ {
		if isBad(v.Xs[i]) || isBad(v.Xs[i+1]) ||
			isBad(v.Y1[i]) || isBad(v.Y2[i]) || isBad(v.Y1[i+1]) || isBad(v.Y2[i+1]) {
			continue
		}
		xa := int(math.Round(tr.X(v.Xs[i])))
		xb := int(math.Round(tr.X(v.Xs[i+1])))
		for x := xa; x <= xb; x++ {
			var t float64
			if xb != xa {
				t = float64(x-xa) / float64(xb-xa)
			}
			ylo := v.Y1[i] + t*(v.Y1[i+1]-v.Y1[i])
			yhi := v.Y2[i] + t*(v.Y2[i+1]-v.Y2[i])
			p0 := int(math.Round(tr.Y(ylo)))
			p1 := int(math.Round(tr.Y(yhi)))
			if p1 < p0 {
				p0, p1 = p1, p0
			}
			for y := p0; y <= p1; y++ {
				if inClip(x, y, clip) {
					c.blend(x, y, v.Color)
				}
			}
		}
	}
}

// drawLegend renders a boxed legend in the top-right corner listing every
// labelled series.
func (c *canvas) drawLegend(a *Axes, w int, tr transform) {
	type entry struct {
		label string
		col   color.RGBA
	}
	var entries []entry
	maxLbl := 0
	for _, s := range a.series {
		lbl := s.legendLabel()
		if lbl == "" {
			continue
		}
		entries = append(entries, entry{lbl, s.styleColor()})
		if tw := textWidth(lbl, 1); tw > maxLbl {
			maxLbl = tw
		}
	}
	if len(entries) == 0 {
		return
	}
	rowH := 14
	boxW := 12 + 18 + maxLbl
	boxH := rowH*len(entries) + 6
	x1 := int(tr.px1) - 6
	x0 := x1 - boxW
	y0 := int(tr.py0) + 6
	y1 := y0 + boxH
	c.fillRect(x0, y0, x1, y1, White)
	c.rectOutline(x0, y0, x1, y1, Black)
	for i, e := range entries {
		ly := y0 + 4 + i*rowH
		col := e.col
		col.A = 0xff
		c.fillRect(x0+6, ly+3, x0+22, ly+7, col)
		c.text(x0+26, ly, e.label, Black, 1)
	}
}

// --- small helpers ------------------------------------------------------

func abs(x int) int {
	if x < 0 {
		return -x
	}
	return x
}

// inClip reports whether (x,y) lies within the clip rectangle.
func inClip(x, y int, clip image.Rectangle) bool {
	return x >= clip.Min.X && x < clip.Max.X && y >= clip.Min.Y && y < clip.Max.Y
}
