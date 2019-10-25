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
	// knots cannot be self-intersecting
	// so all frequency pairs must be co-prime
	return gcd(fx, fy) == 1 && gcd(fx, fz) == 1 && gcd(fy, fz) == 1
}

// randomValidFrequencies returns random frequencies for x, y, and z
// the frequencies are pairwise co-prime and are in [1, maxFrequency]
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

// phaseScore returns a measure of how close to self-intersecting the knot
// with the specified parameters will be (lower is better)
func phaseScore(fx, fy, fz int, px, py, pz float64) float64 {
	// the following quantities may not be integer multiples of pi
	//   fx * py - fy * px
	//   fy * pz - fz * py
	//   fz * px - fx * pz
	// this "score" computes how far away from being a multiple these
	// quantities are, as a measure of how close to self-intersecting
	// the knot will be
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

// bestPhaseShifts returns a list of (px, py) tuples representing
// phase shifts with the best score (multiple results in case of ties)
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

func randomPhaseShifts(fx, fy, fz, divisor int) (px, py float64) {
	ps := bestPhaseShifts(fx, fy, fz, phaseDivisor)
	p := ps[rand.Intn(len(ps))]
	return p.X, p.Y
}

type Knot struct {
	FrequencyX  int
	FrequencyY  int
	FrequencyZ  int
	PhaseShiftX float64
	PhaseShiftY float64
	PhaseShiftZ float64
}

func NewKnot(fx, fy, fz int, px, py, pz float64) *Knot {
	return &Knot{fx, fy, fz, px, py, pz}
}

func NewRandomKnot(maxFrequency, phaseDivisor int) *Knot {
	// pick some random frequencies
	fx, fy, fz := randomValidFrequencies(maxFrequency)

	// normalize X and Y order
	// (Z is special as we force its phase shift to zero below)
	if fy < fx {
		fx, fy = fy, fx
	}

	// find best phase shifts for these frequencies
	px, py := randomPhaseShifts(fx, fy, fz, phaseDivisor)

	// sometimes it may be desireable to use other phase shifts
	// in that case the code below may be used instead
	// for {
	// 	// px = rand.Float64()
	// 	// py = rand.Float64()
	// 	px = float64(rand.Intn(phaseDivisor)) / float64(phaseDivisor)
	// 	py = float64(rand.Intn(phaseDivisor)) / float64(phaseDivisor)
	// 	if phaseScore(fx, fy, fz, px, py, 0) < 1 {
	// 		break
	// 	}
	// }

	return &Knot{fx, fy, fz, px, py, 0}
}

func (k *Knot) Score() float64 {
	return phaseScore(
		k.FrequencyX, k.FrequencyY, k.FrequencyZ,
		k.PhaseShiftX, k.PhaseShiftY, k.PhaseShiftZ)
}

// Position returns the X, Y, Z position of the curve at time t
func (k *Knot) Position(t float64) Vector {
	x := math.Cos(float64(k.FrequencyX)*t + k.PhaseShiftX*2*math.Pi)
	y := math.Cos(float64(k.FrequencyY)*t + k.PhaseShiftY*2*math.Pi)
	z := math.Cos(float64(k.FrequencyZ)*t + k.PhaseShiftZ*2*math.Pi)
	return Vector{x, y, z}
}

// Derivative returns the first derivative at time t
func (k *Knot) Derivative(t float64) Vector {
	x := -math.Sin(float64(k.FrequencyX)*t + k.PhaseShiftX*2*math.Pi)
	y := -math.Sin(float64(k.FrequencyY)*t + k.PhaseShiftY*2*math.Pi)
	z := -math.Sin(float64(k.FrequencyZ)*t + k.PhaseShiftZ*2*math.Pi)
	return Vector{x, y, z}.Normalize()
}

// CrossSectionAt computes the cross-sectional profile of the tube at time t
// the points are stored in the `result` buffer
func (k *Knot) CrossSectionAt(t float64, profile, result []Vector) {
	p := k.Position(t)
	d := k.Derivative(t)
	up := p.Normalize()
	u := up.Cross(d).Normalize()
	v := d.Cross(u)
	for i, q := range profile {
		result[i] = p.Add(u.MulScalar(q.X)).Add(v.MulScalar(q.Y))
	}
}

// Mesh computes a 3D mesh for this Knot
// n = number of slices from 0 to 2 pi
func (k *Knot) Mesh(n int, profile []Vector) *Mesh {
	m := len(profile)
	triangles := make([]*Triangle, 0, n*m*2)
	c0 := make([]Vector, m)
	c1 := make([]Vector, m)
	k.CrossSectionAt(0, profile, c0)
	for i := 0; i < n; i++ {
		t := float64(i+1) / float64(n) * 2 * math.Pi
		k.CrossSectionAt(t, profile, c1)
		for j0 := 0; j0 < m; j0++ {
			j1 := (j0 + 1) % m
			v00 := c0[j0]
			v10 := c1[j0]
			v01 := c0[j1]
			v11 := c1[j1]
			triangles = append(triangles, NewTriangleForPoints(v11, v10, v00))
			triangles = append(triangles, NewTriangleForPoints(v01, v11, v00))
		}
		c1, c0 = c0, c1
	}
	return NewTriangleMesh(triangles)
}

func (k *Knot) Name() string {
	fx := k.FrequencyX
	fy := k.FrequencyY
	fz := k.FrequencyZ
	px := int(math.Round(k.PhaseShiftX * phaseDivisor))
	py := int(math.Round(k.PhaseShiftY * phaseDivisor))
	return fmt.Sprintf("%d.%d.%d.%d.%d", fx, fy, fz, px, py)
}

func ellipticalProfile(n int, a0, rx, ry float64) []Vector {
	result := make([]Vector, n)
	for i := 0; i < n; i++ {
		a := a0 + float64(i)/float64(n)*2*math.Pi
		result[i].X = math.Cos(a) * rx
		result[i].Y = math.Sin(a) * ry
	}
	return result
}

func fileExists(path string) bool {
	_, err := os.Stat(path)
	return !os.IsNotExist(err)
}

func main() {
	profile := ellipticalProfile(tubeSectionSteps, 0, tubeRadius, tubeRadius)

	args := os.Args[1:]
	if len(args) == 0 {
		for {
			k := NewRandomKnot(maxFrequency, phaseDivisor)
			name := k.Name()
			path := name + ".stl"
			if fileExists(path) {
				fmt.Println("SKIP")
				continue
			}
			mesh := k.Mesh(tubeSteps, profile)
			mesh.SaveSTL(path)
			fmt.Println(name, k.Score())
		}
	}

	a := ParseFloats(args)
	fx := int(a[0])
	fy := int(a[1])
	fz := int(a[2])

	var px, py float64
	if len(args) > 3 {
		px = a[3] / phaseDivisor
		py = a[4] / phaseDivisor
	} else {
		px, py = randomPhaseShifts(fx, fy, fz, phaseDivisor)
	}

	k := Knot{fx, fy, fz, px, py, 0}
	name := k.Name()
	mesh := k.Mesh(tubeSteps, profile)
	path := name + ".stl"
	mesh.SaveSTL(path)
	fmt.Println(name, k.Score())
}
