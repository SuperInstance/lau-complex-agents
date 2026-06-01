# lau-complex-agents

**Complex analysis and holomorphic dynamics for agents** — contour integration, residue theory, Riemann surfaces, Julia sets as decision boundaries, and potential-theoretic navigation.

113 tests · MIT license · `nalgebra` + `num-complex`

---

## What This Does

This crate treats the complex plane as an **agent decision landscape**. An `AgentComplex` positions itself in ℂ, classifies its location relative to fractal decision boundaries (Julia/Mandelbrot sets), navigates using conformal maps and potential gradients, and reasons about its environment through the full machinery of complex analysis:

- **Holomorphic function theory** — Cauchy-Riemann checks, power series, analytic continuation
- **Contour integration** — Cauchy's integral formula, winding numbers, Cauchy's theorem verification
- **Residue calculus** — pole classification, Laurent series, argument principle, Rouché's theorem
- **Conformal mapping** — Möbius transformations, Schwarz-Christoffel, Riemann mapping approximation
- **Harmonic analysis** — mean value property, Dirichlet problem (disk + rectangle), Poisson kernel, Green's function
- **Entire functions** — Liouville's theorem, order/type computation, Weierstrass factorization, Hadamard factorization
- **Holomorphic dynamics** — Julia/Mandelbrot sets, Fatou components, fixed point classification, distance estimation
- **Potential theory** — logarithmic potential, equilibrium measures, Fekete points, Dirichlet energy
- **Riemann surfaces** — branch points/cuts, monodromy, covering spaces, Euler characteristic

---

## Key Idea

The central metaphor: **a complex-analytic function IS a decision boundary**. An agent's position in the complex plane is classified by whether it lies in the filled Julia set (in-set → bounded orbit → "safe"), the basin of infinity (escape → "detected"), or on the boundary itself. Conformal maps transform decision landscapes into canonical domains; potential theory provides gradient-based navigation.

This turns classical complex analysis into a computational toolkit for spatial reasoning.

---

## Install

```toml
[dependencies]
lau-complex-agents = "0.1.0"
```

Or clone directly:

```bash
git clone https://github.com/SuperInstance/lau-complex-agents.git
cd lau-complex-agents
cargo build
```

Dependencies: `nalgebra = "0.33"`, `num-complex = "0.4"` (with serde), `serde = "1"`.

---

## Quick Start

```rust
use lau_complex_agents::{AgentComplex, ReasoningReport, BoundaryClassification};
use num_complex::Complex64;

fn c(re: f64, im: f64) -> Complex64 {
    Complex64::new(re, im)
}

fn main() {
    // Create an agent on a Julia set boundary
    let mut agent = AgentComplex::with_julia_boundary(
        c(-0.7, 0.27015),  // Julia parameter
        c(0.0, 0.0),       // starting position
    );

    // Classify: is this point in the filled Julia set?
    let classification = agent.classify_position(200);
    println!("In set: {}, boundary: {}", classification.in_set, classification.boundary_type);

    // Full complex-analytic reasoning
    agent.setup_equilibrium_measure(vec![c(0.0, 0.0), c(1.0, 0.0)]);
    let report = agent.reason(&|z| z.exp(), 200);
    println!("Holomorphic: {}", report.holomorphic.is_holomorphic);
    println!("Singularity: {}", report.singularity_type);
}
```

---

## API Reference

### Top-Level Exports

| Type | Module | Description |
|------|--------|-------------|
| `AgentComplex` | `agent_complex` | Main agent struct — positions in ℂ, navigates decision boundaries |
| `ReasoningReport` | `agent_complex` | Full analysis: holomorphicity, singularity, residue, potential |
| `BoundaryClassification` | `agent_complex` | In-set / escape-iteration classification result |

### `AgentComplex`

The central type. Wraps all complex analysis into an agent API.

```rust
pub struct AgentComplex {
    pub boundary: DecisionBoundary,
    pub state: AgentState,
    pub mode: ReasoningMode,
    pub conformal_map: Option<MobiusTransformation>,
    pub measure: Option<DiscreteMeasure>,
}
```

**Constructors:**
- `with_julia_boundary(c: Complex64, initial_position: Complex64)` — Julia set decision boundary
- `with_mandelbrot_boundary(initial_c: Complex64)` — Mandelbrot set decision boundary

**Methods:**
- `classify_position(max_iterations) → BoundaryClassification` — classify current position
- `compute_orbit(max_iterations) → Vec<Complex64>` — iterate quadratic map from position
- `analyze_holomorphic(f) → HolomorphicAnalysis` — check CR equations, compute derivative
- `contour_integrate(f, radius, n_points) → Complex64` — circle contour integral around position
- `compute_residue(f, pole_order) → Complex64` — residue at current position
- `classify_singularity(f) → SingularityType` — removable/pole/essential classification
- `apply_conformal_map() → Complex64` — apply stored Möbius transformation
- `setup_mobius_map(a, b, c, d)` — set Möbius transformation coefficients
- `setup_equilibrium_measure(points)` — set potential-theoretic measure
- `compute_potential() → f64` — evaluate potential at current position
- `gradient_navigate(step_size) → Complex64` — move along potential gradient
- `check_harmonic_region(u, radius) → bool` — test mean value property
- `analyze_dynamics(c) → DynamicsAnalysis` — fixed points + multipliers of quadratic map
- `riemann_structure(n_branches) → RiemannSurface` — n-th root Riemann surface
- `reason(f, max_iterations) → ReasoningReport` — full reasoning pipeline

### Supporting Types

**`AgentState`** — position, velocity, potential, iteration count, basin label, boundary flag

**`DecisionBoundary`** — enum: `JuliaSet { c, resolution }`, `MandelbrotSet { resolution }`, `Custom { name }`

**`ReasoningMode`** — enum: `Holomorphic`, `Contour`, `Residue`, `Potential`, `Dynamics`, `Conformal`

**`BoundaryClassification`** — `in_set: bool`, `escape_iteration: Option<usize>`, `boundary_type: String`

**`HolomorphicAnalysis`** — `is_holomorphic: bool`, `derivative: Complex64`, `value: Complex64`

**`DynamicsAnalysis`** — `fixed_points`, `multipliers`, `fixed_point_types` (vectors)

**`ReasoningReport`** — combines classification, holomorphic analysis, singularity type, residue, potential

### Module: `holomorphic`

| Function | Description |
|----------|-------------|
| `cauchy_riemann_check(u, v, x, y, h) → bool` | Numerically verify Cauchy-Riemann equations |
| `complex_derivative(f, z, h) → Complex64` | Central-difference complex derivative |
| `is_holomorphic(f, z, h) → bool` | Holomorphicity test via CR + differentiability |
| `analytic_continuation(series, new_center, terms) → PowerSeries` | Re-center a power series by Taylor expansion |

**`PowerSeries`** — stores `center: Complex64` and `coefficients: Vec<Complex64>`:
- `evaluate(z)`, `radius_of_convergence()`, `differentiate()`, `integrate()`

### Module: `contour`

| Type | Description |
|------|-------------|
| `trait Contour` | `evaluate(t: f64) → Complex64`, `derivative(t: f64) → Complex64` |
| `CircleContour` | Circle centered at a point with given radius |
| `LineSegment` | Straight line from start to end |
| `PolygonContour` | Closed polygon through vertices |

| Function | Description |
|----------|-------------|
| `contour_integral(f, contour, n) → Complex64` | Trapezoidal-rule contour integral |
| `cauchy_integral_formula(f, a, contour, n) → Complex64` | CIF: f(a) via contour integral |
| `cauchy_nth_derivative(f, a, n_order, contour, n_points) → Complex64` | CIF for nth derivative |
| `cauchy_theorem_check(f, contour, n, tol) → bool` | Verify ∮f = 0 for holomorphic f |
| `winding_number(contour, a, n) → i64` | Winding number of contour around point |

### Module: `residue`

| Type/Function | Description |
|---------------|-------------|
| `SingularityType` | Enum: `Removable`, `Pole { order }`, `Essential`, `NotIsolated` |
| `residue(f, z0, order, h) → Complex64` | Compute residue at pole (simple or higher-order) |
| `classify_singularity(f, z0, h) → SingularityType` | Classify singularity type |
| `laurent_coefficients(f, z0, r_inner, r_outer, n_neg, n_pos) → (Vec, Vec)` | Compute Laurent series coefficients |
| `argument_principle(f, f', contour, n) → i64` | Number of zeros minus poles inside contour |
| `rouche_check(f, g, contour, n) → bool` | Verify Rouché's theorem conditions |

### Module: `conformal_map`

**`MobiusTransformation`** — stores `a, b, c, d: Complex64`:
- `identity()`, `from_three_points(z1,z2,z3,w1,w2,w3)`, `disk_to_halfplane()`
- `apply(z)`, `inverse()`, `compose(other)`, `determinant()`
- `fixed_points() → Vec<Complex64>`, `preserves_unit_circle() → bool`

**`SchwarzChristoffel`** — prevertices, interior angles, integration constant:
- `to_rectangle(a)`, `integrand(z)`, `evaluate(z, n_steps)`

| Function | Description |
|----------|-------------|
| `cross_ratio(z1, z2, z3, z4) → Complex64` | Cross-ratio (Möbius invariant) |
| `approximate_riemann_mapping(boundary, center, n_terms) → Vec<Complex64>` | Laurent coefficient fit |

### Module: `harmonic`

| Function | Description |
|----------|-------------|
| `mean_value_property(u, x0, y0, r, n) → (f64, f64)` | Returns (center value, circle average) |
| `is_harmonic(u, x0, y0, r, tol) → bool` | Mean value property check |
| `laplacian(u, x, y, h) → f64` | Numerical Laplacian ∇²u |
| `poisson_kernel(r, theta) → f64` | Poisson kernel for unit disk |
| `dirichlet_disk(boundary_fn, r, theta, n) → f64` | Solve Dirichlet problem on disk |
| `dirichlet_rectangle(nx, ny, dx, dy, boundary, max_iter, tol) → DMatrix` | Gauss-Seidel Dirichlet solver |
| `harmonic_conjugate(u, x0, y0, x, y, h) → f64` | Find harmonic conjugate v via path integration |
| `greens_function_disk(z, ζ) → f64` | Green's function for unit disk |

### Module: `entire`

| Function/Type | Description |
|---------------|-------------|
| `compute_order(f, radii, n_theta) → f64` | Order ρ via lim sup log log M(r) / log r |
| `check_liouville(f, radii, n_theta) → (bool, f64)` | Boundedness → const check |
| `compute_type(f, order, radii, n_theta) → f64` | Type σ via lim sup log M(r) / r^ρ |
| `weierstrass_factor(z, p) → Complex64` | Canonical factor E_p(z) |
| `WeierstrassProduct` | Evaluates product with `evaluate(z)`, computes `genus` |
| `hadamard_factorization(zeros, origin_order, order) → WeierstrassProduct` | Construct Hadamard product |
| `fundamental_theorem_of_algebra(coeffs) → Vec<Complex64>` | Durand-Kerner root finding |

### Module: `julia`

| Function/Type | Description |
|---------------|-------------|
| `IterationResult` | Orbit, escaped flag, escape iteration, converged flag |
| `iterate_quadratic(z0, c, max_iter, escape_r) → IterationResult` | Iterate z → z² + c |
| `mandelbrot_membership(c, max_iter, escape_r) → (bool, Option<usize>)` | Is c in Mandelbrot set? |
| `julia_membership(z0, c, max_iter, escape_r) → (bool, Option<usize>)` | Is z0 in Julia set? |
| `filled_julia_set(c, x_range, y_range, res, max_iter, escape_r) → Vec<Vec<usize>>` | Grid computation |
| `mandelbrot_set(x_range, y_range, res, max_iter, escape_r) → Vec<Vec<usize>>` | Grid computation |
| `FatouComponent` | Enum: AttractingFixedPoint, PeriodicCycle, BasinOfInfinity, SiegelDisk, Unknown |
| `quadratic_fixed_points(c) → (Complex64, Complex64)` | Fixed points of z² + c |
| `find_periodic_points(c, period, n_initial) → Vec<Complex64>` | Find periodic cycles |
| `multiplier(f, f', fixed_point) → Complex64` | Multiplier λ = f'(z*) |
| `classify_fixed_point(lambda) → FixedPointType` | SuperAttracting/Attracting/Repelling/Indifferent/Parabolic |
| `mandelbrot_distance(c, max_iter, escape_r) → f64` | Distance estimation to boundary |

### Module: `potential`

| Function/Type | Description |
|---------------|-------------|
| `logarithmic_potential(z, atoms) → f64` | U^μ(z) = Σ w_i log|z - t_i| |
| `potential_gradient(z, atoms) → Complex64` | Gradient of logarithmic potential |
| `DiscreteMeasure` | Weighted atoms with `total_mass()`, `energy()`, `potential(z)`, `robin_constant()`, `capacity()` |
| `unit_circle_equilibrium(n) → DiscreteMeasure` | Uniform measure on unit circle |
| `chebyshev_nodes(n) → Vec<Complex64>` | Chebyshev nodes on [-1, 1] |
| `approximate_fekete_points(n, iterations, lr) → Vec<Complex64>` | Energy-minimization Fekete points |
| `exterior_green_function(z) → f64` | Green's function for complement of unit disk |
| `dirichlet_energy(u, x_range, y_range, nx, ny, h) → f64` | ∫|∇u|² dA |
| `superharmonic_check(u, v, test_points) → bool` | Verify u ≥ v |

### Module: `riemann_surface`

| Type | Description |
|------|-------------|
| `BranchPoint` | location + ramification index |
| `BranchCut` | start/end points, `is_on_cut(z, tol)` |
| `Sheet` | index + label |
| `RiemannSurface` | branch_points, branch_cuts, n_sheets; factory methods: `sqrt_surface()`, `log_surface()`, `nth_root_surface(n)` |
| `Monodromy` | Permutation tracking sheet transitions; `identity(n)`, `transposition(n,i,j)`, `compose()`, `apply()` |
| `CoveringSpace` | base points + fibers |

| Function | Description |
|----------|-------------|
| `euler_characteristic(n_sheets, n_branch_points) → i64` | χ = 2n - b |
| `genus(euler_char) → i64` | g = (2 - χ) / 2 |

---

## How It Works

### Architecture

The crate is organized as independent mathematical modules feeding into the unified `AgentComplex` API:

```
holomorphic ──┐
contour ──────┤
residue ──────┤
conformal_map─┤──→ agent_complex (AgentComplex)
harmonic ─────┤
entire ───────┤
julia ────────┤
potential ────┤
riemann_surface┘
```

Each module is self-contained and can be used independently. `AgentComplex` composes them into a coherent agent that:

1. **Lives in the complex plane** at position `z ∈ ℂ`
2. **Has a decision boundary** (Julia set, Mandelbrot set, or custom)
3. **Classifies its position** by iterating the quadratic map z → z² + c
4. **Navigates** using potential gradients and conformal maps
5. **Reasons** by applying the full toolkit of complex analysis

### Numerical Methods

- **Contour integration**: midpoint trapezoidal rule with configurable resolution
- **Cauchy integral formula**: direct application of contour integration
- **Derivatives**: central finite differences (complex step for Cauchy-Riemann)
- **Residue computation**: Richardson extrapolation for simple poles, Cauchy integral for higher-order poles
- **Laurent series**: numerical integration on annulus
- **Dirichlet problem**: Gauss-Seidel iteration (rectangles), Poisson integral (disk)
- **Root finding**: Durand-Kerner method for polynomial roots
- **Julia/Mandelbrot**: escape-time algorithm with distance estimation
- **Fekete points**: gradient descent on logarithmic energy
- **Singularity classification**: multi-directional limit analysis

---

## The Math

### Cauchy-Riemann Equations

A function f(z) = u(x,y) + iv(x,y) is **holomorphic** when:

```
∂u/∂x = ∂v/∂y    and    ∂u/∂y = -∂v/∂x
```

This crate numerically verifies these at any point using central differences.

### Cauchy's Integral Formula

For holomorphic f inside a contour γ:

```
f(a) = (1/2πi) ∮_γ f(z)/(z-a) dz
```

The nth derivative:

```
f^(n)(a) = n!/(2πi) ∮ f(z)/(z-a)^(n+1) dz
```

### Residue Theorem

```
∮_γ f(z) dz = 2πi Σ Res(f, z_k)
```

where the sum is over poles z_k inside γ. Residues are computed via limits for simple poles and Cauchy integrals for higher-order poles.

### Liouville's Theorem

Every bounded entire function is constant. The crate checks this by sampling M(r) = max_{|z|=r} |f(z)| over multiple radii.

### Order and Type of Entire Functions

```
ρ = lim sup_{r→∞} log log M(r) / log r       (order)
σ = lim sup_{r→∞} log M(r) / r^ρ              (type)
```

e^z has order 1 and type 1; e^{z²} has order 2.

### Möbius Transformations

```
T(z) = (az + b)/(cz + d),   ad - bc ≠ 0
```

These are the conformal automorphisms of the Riemann sphere. They preserve cross-ratios:

```
[z₁, z₂, z₃, z₄] = ((z₁-z₃)(z₂-z₄)) / ((z₁-z₄)(z₂-z₃))
```

### Julia Sets and the Mandelbrot Set

For f_c(z) = z² + c, the **filled Julia set** K_c is the set of z whose orbit remains bounded. The **Mandelbrot set** M is the set of c for which the orbit of 0 remains bounded.

Fixed points: z = z² + c → z² - z + c = 0 → z = (1 ± √(1-4c))/2

Multiplier: λ = f'(z*) = 2z*. Classification:
- |λ| < 1: attracting (|λ| = 0: superattracting)
- |λ| > 1: repelling
- |λ| = 1: indifferent (λ = 1: parabolic)

### Potential Theory

The **logarithmic potential** of a measure μ:

```
U^μ(z) = ∫ log|z-t| dμ(t)
```

The **equilibrium measure** minimizes the energy I(μ) = -∫∫ log|z-t| dμ(z) dμ(t). The **logarithmic capacity** is cap(K) = e^{-I(μ*)} where μ* is the equilibrium measure.

### Riemann Surfaces

For w = z^{1/n}, the surface has n sheets connected at the branch point z = 0. **Monodromy** tracks which sheet you end up on after traversing a loop:

```
sheet → permutation[loop_index](sheet)
```

Euler characteristic: χ = 2n - b (n sheets, b branch points), genus g = (2 - χ)/2.

### Green's Function (Unit Disk)

```
G(z, ζ) = (1/2π)(log|1 - z̄ζ| - log|z - ζ|)
```

Vanishes on the boundary |z| = 1, fundamental solution to Laplace's equation.

---

## License

MIT
