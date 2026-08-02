# Generate curved sequence paths

Generates smooth sequence paths (supporting straight lines, arcs, and
custom curvatures) based on start angle, end angle, radius, and
curvature.

## Usage

``` r
generate_curvature_path(
  start_angle,
  end_angle,
  radius,
  curvature,
  n_points = 100
)
```

## Arguments

- start_angle:

  Numeric, start angle (in radians)

- end_angle:

  Numeric, end angle (in radians)

- radius:

  Numeric, path radius

- curvature:

  Numeric, curvature (0 = straight line, 1 = standard arc, \>1 = more
  curved)

- n_points:

  Integer, number of points in the path (controls smoothness), default
  100

## Value

data.frame containing columns x, y (coordinates of points on the path)
