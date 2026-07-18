package discretemath

import "fmt"

// discretemathErrorf builds a formatted error. It centralizes error creation so
// that every exported routine reports failures in a consistent style.
func discretemathErrorf(format string, args ...any) error {
	return fmt.Errorf(format, args...)
}
