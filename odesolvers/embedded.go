package odesolvers

// HeunEulerTableau returns the embedded Heun-Euler 2(1) pair, the simplest
// adaptive Runge-Kutta method.
func HeunEulerTableau() *ButcherTableau {
	return &ButcherTableau{
		A:         [][]float64{{0, 0}, {1, 0}},
		B:         []float64{0.5, 0.5},
		BStar:     []float64{1, 0},
		C:         []float64{0, 1},
		Order:     2,
		OrderStar: 1,
		Name:      "Heun-Euler 2(1)",
	}
}

// BogackiShampineTableau returns the embedded Bogacki-Shampine 3(2) pair used
// by MATLAB's ode23. It is FSAL (first-same-as-last).
func BogackiShampineTableau() *ButcherTableau {
	return &ButcherTableau{
		A: [][]float64{
			{0, 0, 0, 0},
			{1.0 / 2.0, 0, 0, 0},
			{0, 3.0 / 4.0, 0, 0},
			{2.0 / 9.0, 1.0 / 3.0, 4.0 / 9.0, 0},
		},
		B:         []float64{2.0 / 9.0, 1.0 / 3.0, 4.0 / 9.0, 0},
		BStar:     []float64{7.0 / 24.0, 1.0 / 4.0, 1.0 / 3.0, 1.0 / 8.0},
		C:         []float64{0, 1.0 / 2.0, 3.0 / 4.0, 1},
		Order:     3,
		OrderStar: 2,
		Name:      "Bogacki-Shampine 3(2)",
	}
}

// FehlbergTableau returns the classical Runge-Kutta-Fehlberg RKF45 embedded
// pair. B is the fifth-order solution and BStar the fourth-order companion.
func FehlbergTableau() *ButcherTableau {
	return &ButcherTableau{
		A: [][]float64{
			{0, 0, 0, 0, 0, 0},
			{1.0 / 4.0, 0, 0, 0, 0, 0},
			{3.0 / 32.0, 9.0 / 32.0, 0, 0, 0, 0},
			{1932.0 / 2197.0, -7200.0 / 2197.0, 7296.0 / 2197.0, 0, 0, 0},
			{439.0 / 216.0, -8.0, 3680.0 / 513.0, -845.0 / 4104.0, 0, 0},
			{-8.0 / 27.0, 2.0, -3544.0 / 2565.0, 1859.0 / 4104.0, -11.0 / 40.0, 0},
		},
		B:         []float64{16.0 / 135.0, 0, 6656.0 / 12825.0, 28561.0 / 56430.0, -9.0 / 50.0, 2.0 / 55.0},
		BStar:     []float64{25.0 / 216.0, 0, 1408.0 / 2565.0, 2197.0 / 4104.0, -1.0 / 5.0, 0},
		C:         []float64{0, 1.0 / 4.0, 3.0 / 8.0, 12.0 / 13.0, 1, 1.0 / 2.0},
		Order:     5,
		OrderStar: 4,
		Name:      "RKF45",
	}
}

// CashKarpTableau returns the embedded Cash-Karp 5(4) pair.
func CashKarpTableau() *ButcherTableau {
	return &ButcherTableau{
		A: [][]float64{
			{0, 0, 0, 0, 0, 0},
			{1.0 / 5.0, 0, 0, 0, 0, 0},
			{3.0 / 40.0, 9.0 / 40.0, 0, 0, 0, 0},
			{3.0 / 10.0, -9.0 / 10.0, 6.0 / 5.0, 0, 0, 0},
			{-11.0 / 54.0, 5.0 / 2.0, -70.0 / 27.0, 35.0 / 27.0, 0, 0},
			{1631.0 / 55296.0, 175.0 / 512.0, 575.0 / 13824.0, 44275.0 / 110592.0, 253.0 / 4096.0, 0},
		},
		B:         []float64{37.0 / 378.0, 0, 250.0 / 621.0, 125.0 / 594.0, 0, 512.0 / 1771.0},
		BStar:     []float64{2825.0 / 27648.0, 0, 18575.0 / 48384.0, 13525.0 / 55296.0, 277.0 / 14336.0, 1.0 / 4.0},
		C:         []float64{0, 1.0 / 5.0, 3.0 / 10.0, 3.0 / 5.0, 1, 7.0 / 8.0},
		Order:     5,
		OrderStar: 4,
		Name:      "Cash-Karp",
	}
}

// DormandPrinceTableau returns the embedded Dormand-Prince DOPRI5 5(4) pair
// used by MATLAB's ode45. B is the fifth-order solution.
func DormandPrinceTableau() *ButcherTableau {
	return &ButcherTableau{
		A: [][]float64{
			{0, 0, 0, 0, 0, 0, 0},
			{1.0 / 5.0, 0, 0, 0, 0, 0, 0},
			{3.0 / 40.0, 9.0 / 40.0, 0, 0, 0, 0, 0},
			{44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0, 0, 0, 0},
			{19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0, 0, 0, 0},
			{9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0, 0, 0},
			{35.0 / 384.0, 0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0},
		},
		B:         []float64{35.0 / 384.0, 0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0},
		BStar:     []float64{5179.0 / 57600.0, 0, 7571.0 / 16695.0, 393.0 / 640.0, -92097.0 / 339200.0, 187.0 / 2100.0, 1.0 / 40.0},
		C:         []float64{0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1, 1},
		Order:     5,
		OrderStar: 4,
		Name:      "DOPRI5",
	}
}

// Fehlberg78Tableau returns the high-order Runge-Kutta-Fehlberg RKF78 7(8)
// embedded pair. B is the seventh-order solution and BStar the eighth-order
// companion, so the difference furnishes a seventh-order error estimate.
func Fehlberg78Tableau() *ButcherTableau {
	a := make([][]float64, 13)
	for i := range a {
		a[i] = make([]float64, 13)
	}
	a[1][0] = 2.0 / 27.0
	a[2][0] = 1.0 / 36.0
	a[2][1] = 1.0 / 12.0
	a[3][0] = 1.0 / 24.0
	a[3][2] = 1.0 / 8.0
	a[4][0] = 5.0 / 12.0
	a[4][2] = -25.0 / 16.0
	a[4][3] = 25.0 / 16.0
	a[5][0] = 1.0 / 20.0
	a[5][3] = 1.0 / 4.0
	a[5][4] = 1.0 / 5.0
	a[6][0] = -25.0 / 108.0
	a[6][3] = 125.0 / 108.0
	a[6][4] = -65.0 / 27.0
	a[6][5] = 125.0 / 54.0
	a[7][0] = 31.0 / 300.0
	a[7][4] = 61.0 / 225.0
	a[7][5] = -2.0 / 9.0
	a[7][6] = 13.0 / 900.0
	a[8][0] = 2.0
	a[8][3] = -53.0 / 6.0
	a[8][4] = 704.0 / 45.0
	a[8][5] = -107.0 / 9.0
	a[8][6] = 67.0 / 90.0
	a[8][7] = 3.0
	a[9][0] = -91.0 / 108.0
	a[9][3] = 23.0 / 108.0
	a[9][4] = -976.0 / 135.0
	a[9][5] = 311.0 / 54.0
	a[9][6] = -19.0 / 60.0
	a[9][7] = 17.0 / 6.0
	a[9][8] = -1.0 / 12.0
	a[10][0] = 2383.0 / 4100.0
	a[10][3] = -341.0 / 164.0
	a[10][4] = 4496.0 / 1025.0
	a[10][5] = -301.0 / 82.0
	a[10][6] = 2133.0 / 4100.0
	a[10][7] = 45.0 / 82.0
	a[10][8] = 45.0 / 164.0
	a[10][9] = 18.0 / 41.0
	a[11][0] = 3.0 / 205.0
	a[11][5] = -6.0 / 41.0
	a[11][6] = -3.0 / 205.0
	a[11][7] = -3.0 / 41.0
	a[11][8] = 3.0 / 41.0
	a[11][9] = 6.0 / 41.0
	a[12][0] = -1777.0 / 4100.0
	a[12][3] = -341.0 / 164.0
	a[12][4] = 4496.0 / 1025.0
	a[12][5] = -289.0 / 82.0
	a[12][6] = 2193.0 / 4100.0
	a[12][7] = 51.0 / 82.0
	a[12][8] = 33.0 / 164.0
	a[12][9] = 12.0 / 41.0
	a[12][11] = 1.0
	c := []float64{0, 2.0 / 27.0, 1.0 / 9.0, 1.0 / 6.0, 5.0 / 12.0, 1.0 / 2.0, 5.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0, 1.0 / 3.0, 1.0, 0.0, 1.0}
	// Seventh-order weights.
	b := []float64{41.0 / 840.0, 0, 0, 0, 0, 34.0 / 105.0, 9.0 / 35.0, 9.0 / 35.0, 9.0 / 280.0, 9.0 / 280.0, 41.0 / 840.0, 0, 0}
	// Eighth-order weights.
	bstar := []float64{0, 0, 0, 0, 0, 34.0 / 105.0, 9.0 / 35.0, 9.0 / 35.0, 9.0 / 280.0, 9.0 / 280.0, 0, 41.0 / 840.0, 41.0 / 840.0}
	return &ButcherTableau{
		A:         a,
		B:         b,
		BStar:     bstar,
		C:         c,
		Order:     7,
		OrderStar: 8,
		Name:      "RKF78",
	}
}
