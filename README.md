# lau-complex-analysis

> Complex analysis for agents: holomorphic functions, contour integration, residue calculus, Julia set decision boundaries, and conformal mapping.

## What This Does

This crate provides a comprehensive complex analysis library covering holomorphic function checking (Cauchy-Riemann equations), contour integration (Cauchy's theorem and integral formula), residue calculus (Laurent series, argument principle, Rouché's theorem), Riemann surfaces (branch cuts, monodromy), Möbius transformations and conformal mapping, harmonic functions (Poisson kernel, Dirichlet problem), entire functions (Liouville's theorem, Weierstrass factorization), Julia and Mandelbrot sets as agent decision boundaries, and complex potential theory (logarithmic potential, equilibrium measure).

Use this when you need rigorous complex-analytic tools — computing contour integrals, classifying singularities, mapping domains conformally, or using fractal sets as decision boundaries for agent systems.

## The Key Idea

Complex analysis studies functions of a complex variable that are **holomorphic** (complex differentiable). Unlike real differentiation, complex differentiability is extraordinarily restrictive — it implies infinite differentiability, convergent power series, and the Cauchy integral formula. This rigidity makes complex analysis a powerful computational tool: contour integrals evaluate real integrals, residues sum series, and conformal maps transform difficult domains into easy ones.

## Install

```bash
cargo add lau-complex-analysis
```

## Quick Start

```rust
use lau_complex_analysis::*;
use num_complex::Complex64;

fn c(re: f64, im: f64) -> Complex64 { Complex64::new(re, im) }

fn main() {
    // Check holomorphicity via Cauchy-Riemann equations
    let checker = HolomorphicCheck::new(1e-6);
    let f = |z: Complex64| z * z;  // f(z) = z²
    let result = checker.check_at(&f, 1.0, 1.0);
    println!("f(z)=z² is holomorphic: {}", result.is_holomorphic); // true

    // Contour integral: ∮ 1/z dz = 2πi around origin
    let circle = CircleContour::new(c(0.0, 0.0), 2.0);
    let f = |z: Complex64| 1.0 / z;
    let integral = contour::contour_integral(&f, &circle, 10000);
    println!("∮ 1/z dz = {:.4} (expected 2πi ≈ 6.2832)", integral);

    // Residue at a simple pole
    let f = |z: Complex64| 1.0 / z;
    let res = residue::residue(&f, c(0.0, 0.0), 1, 1e-6);
    println!("Res(1/z, 0) = {}", res); // ≈ 1

    // Mandelbrot set membership
    let (in_set, escape) = julia::mandelbrot_membership(c(0.3, 0.5), 100, 10.0);
    println!("c=0.3+0.5i in Mandelbrot: {}", in_set);

    // Möbius transformation
    let mobius = MobiusTransformation::from_three_points(
        c(0.0, 0.0), c(1.0, 0.0), c(0.0, 1.0),
        c(0.0, 0.0), c(1.0, 0.0), c(0.0, f64::INFINITY),
    );
    println!("T(0) = {}", mobius.apply(c(0.0, 0.0)));
}
```

## API Reference

### Holomorphic Functions

#### `HolomorphicCheck`
Check holomorphicity via Cauchy-Riemann equations using numerical differentiation.

```rust
let checker = HolomorphicCheck::new(1e-6);
let result = checker.check_at(|z| z.powi(3), 1.0, 2.0);
// result.u_x, result.u_y, result.v_x, result.v_y — partial derivatives
// result.cr1_error, result.cr2_error — CR equation violations
// result.is_holomorphic

checker.check_region(&f, &[(0.0, 0.0), (1.0, 1.0)]);
checker.is_entire(&f, (-5.0, 5.0), 20);
HolomorphicCheck::derivative(&f, z);  // f'(z) numerically
```

#### `PowerSeries`
Power series representation of holomorphic functions.

```rust
let series = PowerSeries::from_coefficients(vec![c(1.0, 0.0), c(0.0, 1.0), c(1.0, 0.0)]);
series.evaluate(z);       // evaluate at z
series.radius_of_convergence();
series.truncate(n);       // keep first n terms
```

### Contour Integration

#### Contour Types

```rust
// Circle: γ(t) = center + r·e^{2πit}
let circle = CircleContour::new(center, radius);

// Line segment
let line = LineSegment::new(start, end);

// Polygon
let polygon = PolygonContour::new(vertices);
```

#### Integration Functions

```rust
// Numerical contour integral via midpoint rule
contour_integral(&f, &contour, n_points);

// Cauchy integral formula: f(a) = (1/2πi) ∮ f(z)/(z-a) dz
cauchy_integral_formula(&f, a, &contour, n_points);

// Nth derivative via Cauchy: f^(n)(a) = n!/(2πi) ∮ f(z)/(z-a)^(n+1) dz
cauchy_nth_derivative(&f, a, n_order, &contour, n_points);

// Verify Cauchy's theorem: ∮ f(z) dz = 0 for holomorphic f
cauchy_theorem_check(&f, &contour, n, tolerance);

// Winding number: n(γ, a) = (1/2πi) ∮ 1/(z-a) dz
winding_number(&contour, a, n_points);
```

### Residue Calculus

#### `residue`
Compute the residue of f at z₀.

```rust
// Simple pole (order 1): Res = lim_{z→z0} (z-z0)f(z)
residue(&f, z0, 1, h);

// Higher order pole: Res = (1/(n-1)!) d^{n-1}/dz^{n-1} [(z-z0)^n f(z)]|_{z0}
residue(&f, z0, order, h);
```

#### `classify_singularity`

```rust
let sing = classify_singularity(&f, z0, h);
// SingularityType::Removable | Pole { order } | Essential | NotIsolated
```

#### `laurent_coefficients`

```rust
let (negative, positive) = laurent_coefficients(&f, z0, inner_r, outer_r, n_neg, n_pos);
// f(z) = Σ a_k (z-z0)^k, k from -n_neg to n_pos
```

#### Argument Principle & Rouché's Theorem

```rust
// N - P = (1/2πi) ∮ f'(z)/f(z) dz
let zeros_minus_poles = argument_principle(&f, &f_prime, &contour, n_points);

// If |f-g| < |f| on γ, then f and g have same number of zeros
let same_zeros = rouche_check(&f, &g, &contour, n_points);
```

### Riemann Surfaces

```rust
// Pre-built surfaces
RiemannSurface::sqrt_surface();       // √z: 2 sheets
RiemannSurface::log_surface();        // log(z): ∞ sheets
RiemannSurface::nth_root_surface(n);  // z^{1/n}: n sheets

// Evaluate on specific sheets
RiemannSurface::sqrt_on_sheet(z, 0);   // principal √z
RiemannSurface::sqrt_on_sheet(z, 1);   // -√z
RiemannSurface::log_on_sheet(z, k);    // Log(z) + 2πik

// Branch points and cuts
surface.branch_points;
surface.branch_cuts;
surface.n_sheets;

// Monodromy: track function value along a path
Monodromy::track(&surface, &path, initial_sheet);
```

### Conformal Mapping

#### `MobiusTransformation`
T(z) = (az + b) / (cz + d)

```rust
let t = MobiusTransformation::new(a, b, c, d);
t.apply(z);
t.inverse();
t.compose(&other);
t.fixed_points();
t.determinant();

// Map three points to three points (determines unique Möbius)
MobiusTransformation::from_three_points(z1, z2, z3, w1, w2, w3);
MobiusTransformation::identity();
```

### Harmonic Functions

```rust
// Mean value property: u(x0,y0) = (1/2π) ∫ u(x0+r cos θ, y0+r sin θ) dθ
let (center, mean) = mean_value_property(&u, x0, y0, r, 200);

// Check if harmonic: ∇²u = 0
is_harmonic(&u, x0, y0, r, tolerance);
laplacian(&u, x, y, h);

// Poisson kernel on unit disk
poisson_kernel(r, theta);

// Solve Dirichlet problem on unit disk
dirichlet_disk(&boundary_values, r, theta, n_points);

// Solve Dirichlet problem on rectangle (finite differences)
dirichlet_rectangle(nx, ny, dx, dy, &boundary, max_iter, tolerance);
```

### Entire Functions

```rust
// Order: ρ = lim sup (log log M(r)) / (log r)
let order = compute_order(&f, &radii, n_theta);

// Type: σ = lim sup (log M(r)) / r^ρ
let sigma = compute_type(&f, order, &radii, n_theta);

// Liouville's theorem: entire + bounded ⟹ constant
let (is_constant, max_mod) = check_liouville(&f, &radii, n_theta);

// Weierstrass product for canonical factors
WeierstrassProduct::factor(z, p);
```

### Julia Sets & Mandelbrot

```rust
// Iterate z → z² + c
let result = iterate_quadratic(z0, c, max_iter, escape_radius);
// result.orbit, result.escaped, result.escape_iteration, result.converged

// Set membership
mandelbrot_membership(c, max_iter, escape_r);  // (in_set, escape_iter)
julia_membership(z0, c, max_iter, escape_r);

// Compute grids
let julia_grid = filled_julia_set(c, (-2.0, 2.0), (-2.0, 2.0), 100, 200, 10.0);
let mandelbrot_grid = mandelbrot_set((-2.5, 1.0), (-1.5, 1.5), 100, 200, 10.0);
```

### Potential Theory

```rust
// Logarithmic potential: U^μ(z) = ∫ log|z-t| dμ(t)
logarithmic_potential(z, &atoms);
potential_gradient(z, &atoms);

// Discrete measure
let mu = DiscreteMeasure::uniform(points);
mu.energy();          // I(μ) = -∫∫ log|z-t| dμ(z)dμ(t)
mu.capacity();        // transfinite diameter
mu.potential(z);
mu.robin_constant();

// Equilibrium on unit circle
let eq = unit_circle_equilibrium(100);
```

### Agent Integration

#### `AgentComplex`
Unified API for complex-analytic agent reasoning.

```rust
let mut agent = AgentComplex::with_julia_boundary(c, initial_position);
let classification = agent.classify_position(100);
// classification.in_set, classification.escape_iteration, classification.boundary_type

let orbit = agent.compute_orbit(100);

// With conformal map
agent.conformal_map = Some(mobius);

// With potential theory
agent.measure = Some(measure);
```

#### `BoundaryClassification`

```rust
pub struct BoundaryClassification {
    pub in_set: bool,
    pub escape_iteration: Option<usize>,
    pub boundary_type: String,
}
```

## How It Works

**Holomorphicity** is checked numerically via the Cauchy-Riemann equations: ∂u/∂x = ∂v/∂y and ∂u/∂y = -∂v/∂x. Partial derivatives use finite differences with step h ≈ 10⁻⁷.

**Contour integration** uses the midpoint rule: ∮ f(z) dz ≈ Σ f(γ(tₖ))·γ'(tₖ)·Δt. For Cauchy's theorem, the integral of a holomorphic function around a closed contour is verified to be ≈ 0. The winding number n(γ, a) counts how many times γ winds around a.

**Residues** at simple poles use Richardson extrapolation on the limit lim_{z→z₀} (z-z₀)f(z). Higher-order poles use numerical nth derivatives via Cauchy's integral formula on small circles. Laurent coefficients are computed by numerical integration around an annulus.

**Möbius transformations** are composed via 2×2 matrix multiplication of (a,b;c,d). Three-point determination uses the cross-ratio. Fixed points solve the quadratic cz² + (d-a)z - b = 0.

**Harmonic functions** satisfy ∇²u = 0. The Dirichlet problem on a disk uses the Poisson integral; on a rectangle, Gauss-Seidel iteration solves the finite-difference Laplace equation.

**Julia/Mandelbrot sets** iterate z → z² + c, recording escape time. Points that don't escape within max_iterations are classified as in the set. The filled Julia set grid provides escape-time data for visualization.

**Potential theory** computes logarithmic potentials U^μ(z) = Σ wᵢ·log|z-tᵢ| and energy I(μ) = -ΣΣ wᵢwⱼ·log|zᵢ-tⱼ|. Capacity is exp(-I(μ)).

## The Math

### Cauchy-Riemann Equations

f = u + iv is holomorphic ⟺ u_x = v_y and u_y = -v_x.

### Cauchy Integral Formula

$$f(a) = \frac{1}{2\pi i} \oint_\gamma \frac{f(z)}{z - a} \, dz$$

### Cauchy's Theorem

If f is holomorphic on and inside a closed contour γ:

$$\oint_\gamma f(z) \, dz = 0$$

### Residue Theorem

$$\oint_\gamma f(z) \, dz = 2\pi i \sum_k \text{Res}(f, z_k)$$

where the sum is over poles enclosed by γ.

### Argument Principle

$$N - P = \frac{1}{2\pi i} \oint_\gamma \frac{f'(z)}{f(z)} \, dz$$

### Laurent Series

$$f(z) = \sum_{k=-\infty}^{\infty} a_k (z - z_0)^k, \quad a_k = \frac{1}{2\pi i} \oint \frac{f(z)}{(z-z_0)^{k+1}} \, dz$$

### Liouville's Theorem

Every bounded entire function is constant.

### Riemann Mapping Theorem

Every simply connected domain (≠ ℂ) is conformally equivalent to the unit disk.

### Möbius Transformation

$$T(z) = \frac{az + b}{cz + d}, \quad ad - bc \neq 0$$

Forms the automorphism group of the Riemann sphere.

### Poisson Kernel

$$P(r, \theta) = \frac{1 - r^2}{1 - 2r\cos\theta + r^2}$$

Solves the Dirichlet problem: $u(re^{i\theta}) = \frac{1}{2\pi} \int_0^{2\pi} P(r, \theta - \phi) \, f(e^{i\phi}) \, d\phi$.

### Weierstrass Factorization

$$f(z) = z^m e^{g(z)} \prod_{n=1}^{\infty} E_{p_n}\left(\frac{z}{a_n}\right)$$

where $E_p(z) = (1-z)\exp\left(z + \frac{z^2}{2} + \cdots + \frac{z^p}{p}\right)$.

## License

MIT
