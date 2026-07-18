package plot

import (
	"fmt"
	"image/color"
)

// colorRGBA is a short internal alias for color.RGBA used in signatures.
type colorRGBA = color.RGBA

// Named opaque colors used throughout the package. They mirror the default
// Matplotlib "tab10" style palette so that figures rendered by this package
// look familiar next to Matplotlib output.
var (
	// Blue is the first default series color (Matplotlib C0).
	Blue = color.RGBA{R: 0x1f, G: 0x77, B: 0xb4, A: 0xff}
	// Orange is the second default series color (Matplotlib C1).
	Orange = color.RGBA{R: 0xff, G: 0x7f, B: 0x0e, A: 0xff}
	// Green is the third default series color (Matplotlib C2).
	Green = color.RGBA{R: 0x2c, G: 0xa0, B: 0x2c, A: 0xff}
	// Red is the fourth default series color (Matplotlib C3).
	Red = color.RGBA{R: 0xd6, G: 0x27, B: 0x28, A: 0xff}
	// Purple is the fifth default series color (Matplotlib C4).
	Purple = color.RGBA{R: 0x94, G: 0x67, B: 0xbd, A: 0xff}
	// Brown is the sixth default series color (Matplotlib C5).
	Brown = color.RGBA{R: 0x8c, G: 0x56, B: 0x4b, A: 0xff}
	// Pink is the seventh default series color (Matplotlib C6).
	Pink = color.RGBA{R: 0xe3, G: 0x77, B: 0xc2, A: 0xff}
	// Gray is the eighth default series color (Matplotlib C7).
	Gray = color.RGBA{R: 0x7f, G: 0x7f, B: 0x7f, A: 0xff}
	// Olive is the ninth default series color (Matplotlib C8).
	Olive = color.RGBA{R: 0xbc, G: 0xbd, B: 0x22, A: 0xff}
	// Cyan is the tenth default series color (Matplotlib C9).
	Cyan = color.RGBA{R: 0x17, G: 0xbe, B: 0xcf, A: 0xff}

	// Black is fully opaque black, used for axes, text and borders.
	Black = color.RGBA{R: 0, G: 0, B: 0, A: 0xff}
	// White is fully opaque white, used as the default figure background.
	White = color.RGBA{R: 0xff, G: 0xff, B: 0xff, A: 0xff}
	// LightGray is the default grid line color.
	LightGray = color.RGBA{R: 0xcc, G: 0xcc, B: 0xcc, A: 0xff}
)

// DefaultColors is the ordered cycle of colors assigned to successive series
// when a caller does not set an explicit color. It matches the Matplotlib
// default property cycle.
var DefaultColors = []color.RGBA{
	Blue, Orange, Green, Red, Purple, Brown, Pink, Gray, Olive, Cyan,
}

// ColorCycle returns the color at position i in [DefaultColors], wrapping
// around modulo the palette length. It is the rule this package uses to assign
// automatic colors to new series.
func ColorCycle(i int) color.RGBA {
	n := len(DefaultColors)
	return DefaultColors[((i%n)+n)%n]
}

// hexColor formats c as a CSS "#rrggbb" string for use in SVG output. The
// alpha channel is emitted separately by callers via fill-opacity.
func hexColor(c color.RGBA) string {
	return fmt.Sprintf("#%02x%02x%02x", c.R, c.G, c.B)
}
