# lissaknot

Simple Go program to create 3D models (STLs) of [Lissajous knots](https://en.wikipedia.org/wiki/Lissajous_knot).

### Installation

    $ go get -u github.com/fogleman/lissaknot

### Usage

Running the program with no arguments generates random knots and saves them to disk. See the [source](https://github.com/fogleman/lissaknot/blob/master/main.go) for more details.

    $ go run main.go

You can also specify knot parameters to generate a single knot STL:

    $ go run main.go fx fy fz [px py]

### Examples

![2.3.7](https://www.michaelfogleman.com/static/lissaknots/2.3.7.gif)
![3.4.5](https://www.michaelfogleman.com/static/lissaknots/3.4.5.gif)
![3.5.7](https://www.michaelfogleman.com/static/lissaknots/3.5.7.gif)
![3.8.7](https://www.michaelfogleman.com/static/lissaknots/3.8.7.gif)
![4.5.7](https://www.michaelfogleman.com/static/lissaknots/4.5.7.gif)
![5.6.7](https://www.michaelfogleman.com/static/lissaknots/5.6.7.gif)
