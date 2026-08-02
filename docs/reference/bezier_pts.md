# Generate Bezier curve points

Generates points on a Bezier curve based on start point, end point, and
control points (for smooth ribbons)

## Usage

``` r
bezier_pts(p0, p3, c1, c2, n = 100)
```

## Arguments

- p0:

  Numeric vector (length 2), start point coordinates (x, y)

- p3:

  Numeric vector (length 2), end point coordinates (x, y)

- c1:

  Numeric vector (length 2), first control point coordinates (x, y)

- c2:

  Numeric vector (length 2), second control point coordinates (x, y)

- n:

  Integer, number of curve points (controls smoothness), default 100

## Value

data.frame containing columns x, y (coordinates of points on the Bezier
curve)
