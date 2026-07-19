package queueing

import "math"

// JacksonNetwork models an open Jackson network of M/M/c queues. Node i
// receives external Poisson arrivals at rate External[i], has Servers[i]
// identical servers of rate ServiceRate[i], and routes a completed customer to
// node j with probability Routing[i][j]; the residual probability
// 1 - sum_j Routing[i][j] is the probability of leaving the network. Solving
// the traffic equations yields per-node arrival rates, from which each node is
// analyzed as an independent M/M/c queue (Jackson's theorem).
type JacksonNetwork struct {
	External    []float64   // external arrival rate at each node
	Routing     [][]float64 // routing probabilities Routing[i][j]
	ServiceRate []float64   // per-server service rate at each node
	Servers     []int       // number of servers at each node
}

// NewJacksonNetwork constructs an open [JacksonNetwork]. Servers may be nil, in
// which case every node is single-server. It validates dimensions, sign
// constraints, and that each node's routing probabilities sum to at most one,
// returning an error otherwise.
func NewJacksonNetwork(external []float64, routing [][]float64, serviceRate []float64, servers []int) (JacksonNetwork, error) {
	n := len(external)
	if n == 0 || len(routing) != n || len(serviceRate) != n {
		return JacksonNetwork{}, ErrDimension
	}
	if servers == nil {
		servers = make([]int, n)
		for i := range servers {
			servers[i] = 1
		}
	}
	if len(servers) != n {
		return JacksonNetwork{}, ErrDimension
	}
	for i := 0; i < n; i++ {
		if external[i] < 0 {
			return JacksonNetwork{}, ErrNegative
		}
		if serviceRate[i] <= 0 {
			return JacksonNetwork{}, ErrNonPositiveRate
		}
		if servers[i] <= 0 {
			return JacksonNetwork{}, ErrServers
		}
		if len(routing[i]) != n {
			return JacksonNetwork{}, ErrDimension
		}
		rowSum := 0.0
		for j := 0; j < n; j++ {
			if routing[i][j] < 0 {
				return JacksonNetwork{}, ErrNegative
			}
			rowSum += routing[i][j]
		}
		if rowSum > 1+1e-9 {
			return JacksonNetwork{}, ErrCapacity
		}
	}
	// deep copy
	ext := append([]float64(nil), external...)
	sr := append([]float64(nil), serviceRate...)
	sv := append([]int(nil), servers...)
	rt := make([][]float64, n)
	for i := range routing {
		rt[i] = append([]float64(nil), routing[i]...)
	}
	return JacksonNetwork{External: ext, Routing: rt, ServiceRate: sr, Servers: sv}, nil
}

// Nodes returns the number of nodes in the network.
func (net JacksonNetwork) Nodes() int { return len(net.External) }

// ArrivalRates solves the traffic equations lambda_j = gamma_j + sum_i
// lambda_i r_ij and returns the vector of per-node total arrival rates. It
// returns nil if the linear system is singular.
func (net JacksonNetwork) ArrivalRates() []float64 {
	n := net.Nodes()
	// (I - R^T) lambda = gamma
	a := make([][]float64, n)
	for i := 0; i < n; i++ {
		a[i] = make([]float64, n+1)
		for j := 0; j < n; j++ {
			v := 0.0
			if i == j {
				v = 1
			}
			// R^T[i][j] = Routing[j][i]
			a[i][j] = v - net.Routing[j][i]
		}
		a[i][n] = net.External[i]
	}
	return gaussianSolve(a)
}

// Utilizations returns the per-node utilization rho_i = lambda_i/(c_i mu_i). It
// returns nil if the traffic equations are singular.
func (net JacksonNetwork) Utilizations() []float64 {
	lam := net.ArrivalRates()
	if lam == nil {
		return nil
	}
	rho := make([]float64, len(lam))
	for i := range lam {
		rho[i] = lam[i] / (float64(net.Servers[i]) * net.ServiceRate[i])
	}
	return rho
}

// Stable reports whether every node has utilization strictly below one, the
// condition for the network to admit a product-form steady state.
func (net JacksonNetwork) Stable() bool {
	rho := net.Utilizations()
	if rho == nil {
		return false
	}
	for _, r := range rho {
		if r < 0 || r >= 1 {
			return false
		}
	}
	return true
}

// MeanNumberAt returns the mean number of customers at node i, analyzing it as
// an independent M/M/c queue with the arrival rate from the traffic equations.
// It returns NaN if the node is unstable or the network is singular.
func (net JacksonNetwork) MeanNumberAt(i int) float64 {
	lam := net.ArrivalRates()
	if lam == nil || i < 0 || i >= len(lam) {
		return math.NaN()
	}
	c := net.Servers[i]
	mu := net.ServiceRate[i]
	if lam[i] >= float64(c)*mu {
		return math.NaN()
	}
	if c == 1 {
		q := MM1{Lambda: lam[i], Mu: mu}
		return q.L()
	}
	q := MMc{Lambda: lam[i], Mu: mu, C: c}
	return q.L()
}

// MeanNumbers returns the mean number of customers at every node.
func (net JacksonNetwork) MeanNumbers() []float64 {
	n := net.Nodes()
	out := make([]float64, n)
	for i := 0; i < n; i++ {
		out[i] = net.MeanNumberAt(i)
	}
	return out
}

// TotalMean returns the mean total number of customers in the network,
// summed over all nodes.
func (net JacksonNetwork) TotalMean() float64 {
	sum := 0.0
	for _, l := range net.MeanNumbers() {
		sum += l
	}
	return sum
}

// Throughput returns the network throughput, the total external arrival rate
// (equal in steady state to the total departure rate).
func (net JacksonNetwork) Throughput() float64 {
	sum := 0.0
	for _, g := range net.External {
		sum += g
	}
	return sum
}

// SojournTime returns the mean time a customer spends in the network,
// TotalMean/Throughput, by Little's law applied to the whole network.
func (net JacksonNetwork) SojournTime() float64 {
	thr := net.Throughput()
	if thr <= 0 {
		return math.NaN()
	}
	return net.TotalMean() / thr
}

// gaussianSolve solves the augmented n x (n+1) system by Gaussian elimination
// with partial pivoting and returns the solution vector, or nil if the matrix
// is singular.
func gaussianSolve(a [][]float64) []float64 {
	n := len(a)
	for col := 0; col < n; col++ {
		piv := col
		best := math.Abs(a[col][col])
		for r := col + 1; r < n; r++ {
			if math.Abs(a[r][col]) > best {
				best = math.Abs(a[r][col])
				piv = r
			}
		}
		if best < 1e-14 {
			return nil
		}
		a[col], a[piv] = a[piv], a[col]
		for r := 0; r < n; r++ {
			if r == col {
				continue
			}
			f := a[r][col] / a[col][col]
			for k := col; k <= n; k++ {
				a[r][k] -= f * a[col][k]
			}
		}
	}
	x := make([]float64, n)
	for i := 0; i < n; i++ {
		x[i] = a[i][n] / a[i][i]
	}
	return x
}
