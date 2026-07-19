package networkflow

import "math"

// CostMatrix is a rectangular matrix of assignment costs with rows (workers)
// and columns (jobs). It is the input type for the Hungarian algorithm.
type CostMatrix struct {
	rows int
	cols int
	c    [][]float64
}

// NewCostMatrix returns a rows-by-cols cost matrix initialized to zero. It
// panics if either dimension is negative.
func NewCostMatrix(rows, cols int) *CostMatrix {
	if rows < 0 || cols < 0 {
		panic("networkflow: negative dimension")
	}
	c := make([][]float64, rows)
	for i := range c {
		c[i] = make([]float64, cols)
	}
	return &CostMatrix{rows: rows, cols: cols, c: c}
}

// CostMatrixFrom wraps a rectangular slice as a [CostMatrix] (copying it). It
// returns [ErrDimensionMismatch] if the rows are ragged.
func CostMatrixFrom(m [][]float64) (*CostMatrix, error) {
	rows := len(m)
	cols := 0
	if rows > 0 {
		cols = len(m[0])
	}
	for _, r := range m {
		if len(r) != cols {
			return nil, ErrDimensionMismatch
		}
	}
	cm := NewCostMatrix(rows, cols)
	for i := range m {
		copy(cm.c[i], m[i])
	}
	return cm, nil
}

// Rows returns the number of rows (workers).
func (m *CostMatrix) Rows() int { return m.rows }

// Cols returns the number of columns (jobs).
func (m *CostMatrix) Cols() int { return m.cols }

// Set sets the cost of assigning row i to column j.
func (m *CostMatrix) Set(i, j int, cost float64) { m.c[i][j] = cost }

// At returns the cost of assigning row i to column j.
func (m *CostMatrix) At(i, j int) float64 { return m.c[i][j] }

// Clone returns a deep copy of the cost matrix.
func (m *CostMatrix) Clone() *CostMatrix {
	c := NewCostMatrix(m.rows, m.cols)
	for i := range m.c {
		copy(c.c[i], m.c[i])
	}
	return c
}

// Transpose returns a new cost matrix with rows and columns exchanged.
func (m *CostMatrix) Transpose() *CostMatrix {
	t := NewCostMatrix(m.cols, m.rows)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			t.c[j][i] = m.c[i][j]
		}
	}
	return t
}

// Negate returns a new cost matrix with every entry negated; solving a
// minimum-cost assignment on it maximizes the original weights.
func (m *CostMatrix) Negate() *CostMatrix {
	c := NewCostMatrix(m.rows, m.cols)
	for i := range m.c {
		for j := range m.c[i] {
			c.c[i][j] = -m.c[i][j]
		}
	}
	return c
}

// AssignmentResult describes an optimal assignment.
type AssignmentResult struct {
	// Assignment maps each row to the column it is assigned to, or -1 if the
	// row is unassigned (possible only when there are more rows than columns).
	Assignment []int
	// Cost is the total cost of the assignment.
	Cost float64
}

// Pairs returns the assigned (row, column) pairs.
func (a *AssignmentResult) Pairs() [][2]int {
	var out [][2]int
	for i, j := range a.Assignment {
		if j >= 0 {
			out = append(out, [2]int{i, j})
		}
	}
	return out
}

// hungarian solves the rectangular minimum-cost assignment for a cost matrix
// with rows <= cols. It returns rowAssign (column chosen by each row) and the
// total cost. It uses the O(n^3) Kuhn-Munkres shortest-augmenting-path form.
func hungarian(cost [][]float64) ([]int, float64) {
	n := len(cost)
	if n == 0 {
		return nil, 0
	}
	m := len(cost[0])
	const inf = math.MaxFloat64

	u := make([]float64, n+1)
	v := make([]float64, m+1)
	p := make([]int, m+1) // p[j] = row assigned to column j (1-indexed rows)
	way := make([]int, m+1)

	for i := 1; i <= n; i++ {
		p[0] = i
		j0 := 0
		minv := make([]float64, m+1)
		used := make([]bool, m+1)
		for j := 0; j <= m; j++ {
			minv[j] = inf
		}
		for {
			used[j0] = true
			i0 := p[j0]
			delta := inf
			j1 := -1
			for j := 1; j <= m; j++ {
				if used[j] {
					continue
				}
				cur := cost[i0-1][j-1] - u[i0] - v[j]
				if cur < minv[j] {
					minv[j] = cur
					way[j] = j0
				}
				if minv[j] < delta {
					delta = minv[j]
					j1 = j
				}
			}
			for j := 0; j <= m; j++ {
				if used[j] {
					u[p[j]] += delta
					v[j] -= delta
				} else {
					minv[j] -= delta
				}
			}
			j0 = j1
			if p[j0] == 0 {
				break
			}
		}
		for j0 != 0 {
			j1 := way[j0]
			p[j0] = p[j1]
			j0 = j1
		}
	}

	rowAssign := make([]int, n)
	for i := range rowAssign {
		rowAssign[i] = -1
	}
	for j := 1; j <= m; j++ {
		if p[j] >= 1 && p[j] <= n {
			rowAssign[p[j]-1] = j - 1
		}
	}
	var total float64
	for i := 0; i < n; i++ {
		if rowAssign[i] >= 0 {
			total += cost[i][rowAssign[i]]
		}
	}
	return rowAssign, total
}

// HungarianMinCost solves the minimum-cost assignment problem for the given
// cost matrix using the Hungarian (Kuhn-Munkres) algorithm in O(n^3) time. Each
// row is assigned to a distinct column so as to minimize total cost. When there
// are more rows than columns the matrix is transposed internally so that every
// column (the smaller side) is covered. The input is left unchanged.
func HungarianMinCost(m *CostMatrix) *AssignmentResult {
	if m.rows == 0 || m.cols == 0 {
		return &AssignmentResult{Assignment: make([]int, m.rows), Cost: 0}
	}
	if m.rows <= m.cols {
		assign, cost := hungarian(m.c)
		return &AssignmentResult{Assignment: assign, Cost: cost}
	}
	// More rows than columns: solve on the transpose, then invert the map.
	t := m.Transpose()
	colAssign, cost := hungarian(t.c) // colAssign[j] = row for column j
	rowAssign := make([]int, m.rows)
	for i := range rowAssign {
		rowAssign[i] = -1
	}
	for j, i := range colAssign {
		if i >= 0 {
			rowAssign[i] = j
		}
	}
	return &AssignmentResult{Assignment: rowAssign, Cost: cost}
}

// MaxWeightAssignment solves the maximum-weight assignment problem: it assigns
// rows to distinct columns to maximize total weight, by negating the matrix and
// running [HungarianMinCost]. The returned Cost is the (positive) total weight.
// The input is left unchanged.
func MaxWeightAssignment(m *CostMatrix) *AssignmentResult {
	res := HungarianMinCost(m.Negate())
	res.Cost = -res.Cost
	return res
}

// KuhnMunkres is an alias for [HungarianMinCost].
func KuhnMunkres(m *CostMatrix) *AssignmentResult { return HungarianMinCost(m) }

// AssignmentCost returns the total cost of the given row-to-column assignment
// under matrix m. Rows assigned to -1 contribute nothing.
func AssignmentCost(m *CostMatrix, assignment []int) float64 {
	var total float64
	for i, j := range assignment {
		if j >= 0 {
			total += m.c[i][j]
		}
	}
	return total
}

// IsValidAssignment reports whether assignment maps rows to distinct, in-range
// columns (allowing -1 for unassigned rows).
func IsValidAssignment(m *CostMatrix, assignment []int) bool {
	if len(assignment) != m.rows {
		return false
	}
	seen := make(map[int]bool)
	for _, j := range assignment {
		if j == -1 {
			continue
		}
		if j < 0 || j >= m.cols || seen[j] {
			return false
		}
		seen[j] = true
	}
	return true
}
