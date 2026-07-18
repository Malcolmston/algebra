package controltheory

// PIDController is a discrete-time proportional-integral-derivative controller
// with configurable gains. It maintains internal state (the accumulated
// integral and the previous error) so successive calls to Update implement the
// standard incremental control law.
type PIDController struct {
	// Kp is the proportional gain.
	Kp float64
	// Ki is the integral gain.
	Ki float64
	// Kd is the derivative gain.
	Kd float64

	integral    float64
	prevError   float64
	initialized bool
}

// NewPIDController returns a PIDController with the given proportional,
// integral, and derivative gains and zero internal state.
func NewPIDController(kp, ki, kd float64) *PIDController {
	return &PIDController{Kp: kp, Ki: ki, Kd: kd}
}

// Reset clears the accumulated integral term and derivative history so the
// controller behaves as if freshly constructed.
func (c *PIDController) Reset() {
	c.integral = 0
	c.prevError = 0
	c.initialized = false
}

// Update advances the controller by one time step of length dt given the
// current error (setpoint minus measurement) and returns the control output
//
//	u = Kp·e + Ki·∫e dt + Kd·de/dt
//
// The integral is accumulated with the rectangle rule and the derivative uses a
// backward difference. On the first call after construction or Reset the
// derivative term is taken as zero. It panics if dt is not positive.
func (c *PIDController) Update(errValue, dt float64) float64 {
	if dt <= 0 {
		panic("controltheory: PID time step dt must be positive")
	}
	c.integral += errValue * dt
	var deriv float64
	if c.initialized {
		deriv = (errValue - c.prevError) / dt
	}
	c.prevError = errValue
	c.initialized = true
	return c.Kp*errValue + c.Ki*c.integral + c.Kd*deriv
}

// Integral returns the current accumulated integral of the error.
func (c *PIDController) Integral() float64 {
	return c.integral
}

// TransferFunction returns the ideal continuous-time transfer function of the
// controller, C(s) = Kp + Ki/s + Kd·s = (Kd·s^2 + Kp·s + Ki) / s.
func (c *PIDController) TransferFunction() TransferFunction {
	return TransferFunction{
		Num: Poly{c.Ki, c.Kp, c.Kd},
		Den: Poly{0, 1},
	}
}

// ZieglerNicholsPID returns PID gains tuned by the classical Ziegler-Nichols
// ultimate-sensitivity rule from the ultimate gain ku (the proportional gain at
// which the loop sustains oscillation) and the ultimate period pu (the
// oscillation period, in seconds). The returned controller uses
// Kp = 0.6·ku, Ki = Kp/(0.5·pu), Kd = Kp·0.125·pu.
func ZieglerNicholsPID(ku, pu float64) *PIDController {
	kp := 0.6 * ku
	ti := 0.5 * pu
	td := 0.125 * pu
	return NewPIDController(kp, kp/ti, kp*td)
}
