# Numerical Interpolation Library

## Overview

This Python class provides implementations of various numerical interpolation and finite difference methods. The library can be used to calculate interpolating polynomials, approximate values at intermediate points, and compute derivatives using Newton's forward and backward difference formulas.

## Features

- **Lagrange's Interpolation**
  - Constructs an interpolating polynomial for a given set of points using Lagrange's method.

- **Newton's Divided Difference Interpolation**
  - Computes the polynomial using Newton's divided differences for unequally spaced intervals.

- **Newton's Forward and Backward Difference Interpolation**
  - Handles equally spaced intervals to compute interpolating polynomials.

- **Derivative Calculation**
  - Calculates the nth derivative at a given point using Newton's forward or backward difference formulas.

## Methods

### `lagranges_interpolation()`
- Computes and returns the Lagrange interpolating polynomial for a given set of points (`x` and `fx`).

### `newton_divided_difference()`
- Computes and returns the Newton divided difference polynomial for unequally spaced points.

### `nfb(method="forward")`
- Helper function that generates the polynomial for Newton's forward or backward difference interpolation based on the provided method.
  - **Parameters**:
    - `method`: `"forward"` (default) or `"backward"` to indicate the interpolation method.

### `newton_forward_interpolation(input)`
- Evaluates the Newton forward difference interpolation polynomial at a given point.
  - **Parameters**:
    - `input`: The value at which to evaluate the polynomial.

### `newton_backward_interpolation(input)`
- Evaluates the Newton backward difference interpolation polynomial at a given point.
  - **Parameters**:
    - `input`: The value at which to evaluate the polynomial.

### `newton_forward_diff(nth, input)`
- Computes the nth derivative using Newton's forward difference method.
  - **Parameters**:
    - `nth`: The order of the derivative.
    - `input`: The point at which to evaluate the derivative.

### `newton_backward_diff(nth, input)`
- Computes the nth derivative using Newton's backward difference method.
  - **Parameters**:
    - `nth`: The order of the derivative.
    - `input`: The point at which to evaluate the derivative.

