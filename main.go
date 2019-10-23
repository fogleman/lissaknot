package main

import (
	"fmt"
	"math"
	"math/rand"
	"os"

	. "github.com/fogleman/fauxgl"
)

const (
	maxFrequency = 8
	phaseDivisor = 360

	tubeRadius       = 0.125
	tubeSteps        = 1024
	tubeSectionSteps = 64
)

func gcd(a, b int) int {
	for b != 0 {
		a, b = b, a%b
	}
	return a
}

func validFrequencies(fx, fy, fz int) bool {
	return gcd(fx, fy) == 1 && gcd(fx, fz) == 1 && gcd(fy, fz) == 1
}

func randomValidFrequencies(maxFrequency int) (fx, fy, fz int) {
	for {
		fx = rand.Intn(maxFrequency) + 1
		fy = rand.Intn(maxFrequency) + 1
		fz = rand.Intn(maxFrequency) + 1
		if validFrequencies(fx, fy, fz) {
			return
		}
	}
}

func phaseScore(fx, fy, fz int, px, py, pz float64) float64 {
	ffx := float64(fx)
	ffy := float64(fy)
	ffz := float64(fz)
	ppx := px * 2 * math.Pi
	ppy := py * 2 * math.Pi
	ppz := pz * 2 * math.Pi
	m1 := math.Mod(math.Abs(ffx*ppy-ffy*ppx), math.Pi) / math.Pi
	m2 := math.Mod(math.Abs(ffy*ppz-ffz*ppy), math.Pi) / math.Pi
	m3 := math.Mod(math.Abs(ffz*ppx-ffx*ppz), math.Pi) / math.Pi
	q1 := math.Pow((m1-0.5)*2, 2)
	q2 := math.Pow((m2-0.5)*2, 2)
	q3 := math.Pow((m3-0.5)*2, 2)
	return q1 + q2 + q3
}

func bestPhaseShifts(fx, fy, fz, divisor int) []Vector {
	var result []Vector
	var best = 1e9
	const eps = 1e-9
	for y := 0; y < divisor; y++ {
		for x := 0; x < divisor; x++ {
			px := float64(x) / float64(divisor)
			py := float64(y) / float64(divisor)
			score := phaseScore(fx, fy, fz, px, py, 0)
			delta := score - best
			if delta > eps {
				continue
			}
			if delta < -eps {
				result = result[:0]
				best = score
			}
			result = append(result, Vector{px, py, 0})
		}
	}
	return result
}

type Knot struct {
	FrequencyX int
	FrequencyY int
	FrequencyZ int
	PhaseX     float64
	PhaseY     float64
	PhaseZ     float64
}

func NewKnot(fx, fy, fz int, px, py, pz float64) *Knot {
	return &Knot{fx, fy, fz, px, py, pz}
}

func NewRandomKnot(maxFrequency, phaseDivisor int) *Knot {
	fx, fy, fz := randomValidFrequencies(maxFrequency)
	if fy < fx {
		fx, fy = fy, fx
	}
	ps := bestPhaseShifts(fx, fy, fz, phaseDivisor)
	p := ps[rand.Intn(len(ps))]
	px := p.X
	py := p.Y

	// for {
	// 	px = rand.Float64()
	// 	py = rand.Float64()
	// 	if phaseScore(fx, fy, fz, px, py, 0) < 0.5 {
	// 		break
	// 	}
	// }

	return &Knot{fx, fy, fz, px, py, 0}
}

func (k *Knot) Score() float64 {
	return phaseScore(
		k.FrequencyX, k.FrequencyY, k.FrequencyZ,
		k.PhaseX, k.PhaseY, k.PhaseZ)
}

func (k *Knot) At(t float64) Vector {
	x := math.Cos(float64(k.FrequencyX)*t + k.PhaseX*2*math.Pi)
	y := math.Cos(float64(k.FrequencyY)*t + k.PhaseY*2*math.Pi)
	z := math.Cos(float64(k.FrequencyZ)*t + k.PhaseZ*2*math.Pi)
	return Vector{x, y, z}
}

func (k *Knot) DirectionAt(t float64) Vector {
	const eps = 1e-9
	a := k.At(t - eps)
	b := k.At(t + eps)
	return b.Sub(a).Normalize()
}

func (k *Knot) CrossSectionAt(t, r float64, n int, up Vector, result []Vector) Vector {
	p := k.At(t)
	w := k.DirectionAt(t)
	u := up.Cross(w).Normalize()
	v := w.Cross(u)
	for i := range result {
		a := 2 * math.Pi * float64(i) / float64(n)
		q := p
		q = q.Add(u.MulScalar(math.Cos(a) * r))
		q = q.Add(v.MulScalar(math.Sin(a) * r))
		result[i] = q
	}
	return v
}

func (k *Knot) Mesh(r float64, n, m int) *Mesh {
	triangles := make([]*Triangle, 0, (n+1)*m*2)
	c0 := make([]Vector, m)
	c1 := make([]Vector, m)
	up := k.DirectionAt(0).Perpendicular()
	up = k.CrossSectionAt(0, r, m, up, c0)
	for i := 0; i < n; i++ {
		t := float64(i+1) / float64(n) * 2 * math.Pi
		up = k.CrossSectionAt(t, r, m, up, c1)
		for j0 := 0; j0 < m; j0++ {
			j1 := (j0 + 1) % m
			v00 := c0[j0]
			v10 := c1[j0]
			v01 := c0[j1]
			v11 := c1[j1]
			triangles = append(triangles, NewTriangleForPoints(v11, v10, v00))
			triangles = append(triangles, NewTriangleForPoints(v01, v11, v00))
		}
		copy(c0, c1)
	}
	return NewTriangleMesh(triangles)
}

func fileExists(path string) bool {
	_, err := os.Stat(path)
	return !os.IsNotExist(err)
}

func main() {
	args := os.Args[1:]
	if len(args) == 5 {
		a := ParseFloats(args)
		fx := int(a[0])
		fy := int(a[1])
		fz := int(a[2])
		px := a[3] / (2 * math.Pi)
		py := a[4] / (2 * math.Pi)
		k := Knot{fx, fy, fz, px, py, 0}
		mesh := k.Mesh(tubeRadius, tubeSteps, tubeSectionSteps)
		mesh.SaveSTL("out.stl")
		fmt.Println(k.Score())
		return
	}
	for {
		k := NewRandomKnot(maxFrequency, phaseDivisor)
		fx := k.FrequencyX
		fy := k.FrequencyY
		fz := k.FrequencyZ
		px := int(math.Round(k.PhaseX * phaseDivisor))
		py := int(math.Round(k.PhaseY * phaseDivisor))
		name := fmt.Sprintf("%d.%d.%d.%d.%d", fx, fy, fz, px, py)
		path := name + ".stl"
		if fileExists(path) {
			fmt.Println("SKIP")
			continue
		}
		mesh := k.Mesh(tubeRadius, tubeSteps, tubeSectionSteps)
		mesh.SaveSTL(name + ".stl")
		fmt.Println(name, k.Score())
	}
}
