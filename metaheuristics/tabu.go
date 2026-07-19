package metaheuristics

import "math"

// TabuList is a bounded first-in-first-out set of recently visited moves,
// identified by string keys. A move present in the list is "tabu" (forbidden)
// until it ages out.
type TabuList struct {
	tenure int
	order  []string
	set    map[string]int
}

// NewTabuList creates an empty tabu list that remembers up to tenure moves.
func NewTabuList(tenure int) *TabuList {
	if tenure < 0 {
		tenure = 0
	}
	return &TabuList{tenure: tenure, set: make(map[string]int)}
}

// Contains reports whether the given move key is currently tabu.
func (t *TabuList) Contains(key string) bool {
	_, ok := t.set[key]
	return ok
}

// Add records a move key, evicting the oldest entry when the tenure is
// exceeded.
func (t *TabuList) Add(key string) {
	if t.tenure == 0 {
		return
	}
	if _, ok := t.set[key]; ok {
		return
	}
	t.order = append(t.order, key)
	t.set[key] = len(t.order)
	if len(t.order) > t.tenure {
		old := t.order[0]
		t.order = t.order[1:]
		delete(t.set, old)
	}
}

// Len returns the number of moves currently in the list.
func (t *TabuList) Len() int { return len(t.set) }

// Tenure returns the maximum number of moves the list retains.
func (t *TabuList) Tenure() int { return t.tenure }

// Clear empties the tabu list.
func (t *TabuList) Clear() {
	t.order = t.order[:0]
	t.set = make(map[string]int)
}

// TabuMove describes a candidate neighbour in a combinatorial tabu search: the
// resulting state, its objective value, and a string key identifying the move
// (used for tabu bookkeeping).
type TabuMove struct {
	State []int
	Value float64
	Key   string
}

// TabuSearchConfig configures [TabuSearchCombinatorial].
type TabuSearchConfig struct {
	// Tenure is the tabu list length.
	Tenure int
	// Iterations is the number of search steps.
	Iterations int
	// RecordHistory enables per-iteration best-value recording.
	RecordHistory bool
}

// TabuSearchCombinatorial minimizes an objective over a discrete state space
// (states represented as int slices). At each step it evaluates the neighbours
// produced by neighbors and moves to the best non-tabu neighbour, unless a tabu
// neighbour beats the incumbent best (aspiration criterion). It returns the best
// state found and its value.
func TabuSearchCombinatorial(
	initial []int,
	neighbors func(state []int) []TabuMove,
	cfg TabuSearchConfig,
) (best []int, bestValue float64, history []float64, err error) {
	if cfg.Iterations <= 0 {
		return nil, 0, nil, ErrInvalidConfig
	}
	tabu := NewTabuList(cfg.Tenure)
	cur := append([]int(nil), initial...)
	best = append([]int(nil), cur...)
	// Evaluate incumbent by probing a single-element neighbourhood is not
	// possible generically; require neighbors to include current value via the
	// first move's provenance. Instead compute best from the first neighbour
	// scan by tracking as we go.
	bestValue = math.Inf(1)
	for it := 0; it < cfg.Iterations; it++ {
		moves := neighbors(cur)
		if len(moves) == 0 {
			break
		}
		chosen := -1
		var chosenVal float64
		for i, mv := range moves {
			isTabu := tabu.Contains(mv.Key)
			aspirate := mv.Value < bestValue
			if isTabu && !aspirate {
				continue
			}
			if chosen < 0 || mv.Value < chosenVal {
				chosen = i
				chosenVal = mv.Value
			}
		}
		if chosen < 0 {
			// all neighbours tabu; pick overall best to keep moving
			for i, mv := range moves {
				if i == 0 || mv.Value < chosenVal {
					chosen = i
					chosenVal = mv.Value
				}
			}
		}
		mv := moves[chosen]
		cur = append([]int(nil), mv.State...)
		tabu.Add(mv.Key)
		if mv.Value < bestValue {
			bestValue = mv.Value
			best = append([]int(nil), mv.State...)
		}
		if cfg.RecordHistory {
			history = append(history, bestValue)
		}
	}
	return best, bestValue, history, nil
}

// ContinuousTabuConfig configures [ContinuousTabuSearch].
type ContinuousTabuConfig struct {
	// Bounds is the search box.
	Bounds Bounds
	// Iterations is the number of search steps.
	Iterations int
	// Neighbors is the number of candidate moves sampled per step.
	Neighbors int
	// StepSize is the relative neighbourhood radius as a fraction of the box
	// width.
	StepSize float64
	// Tenure is the tabu list length; a move is tabu if its rounded grid cell
	// was recently visited.
	Tenure int
	// GridResolution controls how finely positions are quantized into tabu
	// cells; larger values make more positions distinct.
	GridResolution int
	// RecordHistory enables per-iteration best-value recording.
	RecordHistory bool
}

// DefaultContinuousTabuConfig returns a reasonable configuration for the box.
func DefaultContinuousTabuConfig(b Bounds) ContinuousTabuConfig {
	return ContinuousTabuConfig{
		Bounds:         b,
		Iterations:     1000,
		Neighbors:      20,
		StepSize:       0.1,
		Tenure:         30,
		GridResolution: 50,
	}
}

func (c ContinuousTabuConfig) validate() error {
	if !c.Bounds.Valid() {
		return ErrEmptyBounds
	}
	if c.Iterations <= 0 || c.Neighbors <= 0 || c.StepSize <= 0 {
		return ErrInvalidConfig
	}
	return nil
}

// ContinuousTabuSearch minimizes f over the box by, at each step, sampling
// several Gaussian neighbours and moving to the best whose quantized grid cell
// is not tabu (with aspiration when it improves on the incumbent). It is
// deterministic given rng.
func ContinuousTabuSearch(f ObjectiveFunc, cfg ContinuousTabuConfig, start []float64, rng *RNG) (Result, error) {
	if err := cfg.validate(); err != nil {
		return Result{}, err
	}
	if start == nil {
		start = cfg.Bounds.Center()
	}
	if len(start) != cfg.Bounds.Dim() {
		return Result{}, ErrDimMismatch
	}
	width := cfg.Bounds.Width()
	res := gridResGuard(cfg.GridResolution)
	tabu := NewTabuList(cfg.Tenure)
	cur := cfg.Bounds.Clip(start)
	curF := f(cur)
	best := VecCopy(cur)
	bestF := curF
	evals := 1
	out := Result{}
	tabu.Add(gridKey(cur, cfg.Bounds, res))
	for it := 0; it < cfg.Iterations; it++ {
		var chosen []float64
		chosenF := math.Inf(1)
		chosenKey := ""
		for k := 0; k < cfg.Neighbors; k++ {
			cand := make([]float64, len(cur))
			for i := range cand {
				cand[i] = cur[i] + rng.NormFloat64()*cfg.StepSize*width[i]
			}
			cfg.Bounds.ClipInPlace(cand)
			cf := f(cand)
			evals++
			key := gridKey(cand, cfg.Bounds, res)
			if tabu.Contains(key) && cf >= bestF {
				continue
			}
			if cf < chosenF {
				chosen = cand
				chosenF = cf
				chosenKey = key
			}
		}
		if chosen == nil {
			continue
		}
		cur = chosen
		curF = chosenF
		tabu.Add(chosenKey)
		if curF < bestF {
			best = VecCopy(cur)
			bestF = curF
		}
		if cfg.RecordHistory {
			out.History = append(out.History, bestF)
		}
	}
	out.X = best
	out.F = bestF
	out.Iterations = cfg.Iterations
	out.Evaluations = evals
	return out, nil
}

func gridResGuard(r int) int {
	if r < 1 {
		return 1
	}
	return r
}

func gridKey(x []float64, b Bounds, res int) string {
	buf := make([]byte, 0, len(x)*3)
	for i := range x {
		w := b.Upper[i] - b.Lower[i]
		var cell int
		if w > 0 {
			cell = int(math.Floor((x[i] - b.Lower[i]) / w * float64(res)))
		}
		buf = appendInt(buf, cell)
		buf = append(buf, ',')
	}
	return string(buf)
}

func appendInt(b []byte, v int) []byte {
	if v == 0 {
		return append(b, '0')
	}
	if v < 0 {
		b = append(b, '-')
		v = -v
	}
	var tmp [20]byte
	i := len(tmp)
	for v > 0 {
		i--
		tmp[i] = byte('0' + v%10)
		v /= 10
	}
	return append(b, tmp[i:]...)
}
