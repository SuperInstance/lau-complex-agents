# lau-complex-agents

**Complex analysis and holomorphic dynamics for autonomous agents.**

A Rust library that provides the full machinery of complex analysis — contour integration, residue calculus, Riemann surfaces, conformal mapping, harmonic functions, entire functions, Julia/Mandelbrot dynamics, and potential theory — unified under a single `AgentComplex` API designed for agent decision-making over complex decision landscapes.

[![113 tests passing](https://img.shields.io/badge/tests-113%20passing-brightgreen)]()

---

## What This Does

This crate treats the complex plane as an agent's decision space. Fractal boundaries (Julia sets, Mandelbrot set) become decision boundaries. Conformal maps transform agent domains. Potential theory provides gradient-based navigation. Holomorphic analysis classifies the local geometry of decision functions.

In practice: you create an `AgentComplex`, give it a decision boundary (Julia set, Mandelbrot, or custom), position it in the complex plane, and then use 2,500+ years of complex analysis to reason about where it is and where it should go.

## Key Idea

**Complex analysis is the natural language for agent reasoning on continuous domains.** A holomorphic function preserves angles (conformality), and this local structure gives agents rich geometric information. Contour integrals measure global properties. Residues capture singular behavior. Julia sets partition the plane into basins of attraction — natural decision boundaries.

The library implements this as composable building blocks, each independently usable, then unifies them under `AgentComplex` for a single high-level API.

## Install

```toml
[dependencies]
lau-complex-agents = "0.1.0"
```

Requires Rust 2021 edition. Dependencies: `nalgebra`, `num-complex` (with serde), `serde`.

## Quick Start

```rust
use lau_complex_agents::{AgentComplex, ReasoningReport};
use num_complex::Complex64;

fn c(re: f64, im: f64) -> Complex64 { Complex64::new(re, im) }

fn main() {
    // Create an agent with a Julia set decision boundary
    let mut agent = AgentComplex::with_julia_boundary(
        c(-0.7, 0.27015),  // Julia parameter c
        c(0.0, 0.0),        // starting position
    );

    // Classify current position: inside or outside the Julia set?
    let classification = agent.classify_position(200);
    println!("In set: {}, escape: {:?}", classification.in_set, classification.escape_iteration);

    // Full reasoning report: holomorphic analysis, singularities, residues
    let report = agent.reason(&|z| z.exp(), 100);
    println!("Holomorphic at position: {}", report.holomorphic.is_holomorphic);
    println!("Singularity type: {}", report.singularity_type);
    println!("Residue: {:?}", report.residue);

    // Set up potential field and navigate
    agent.setup_equilibrium_measure(vec![c(0.0, 0.0), c(1.0, 0.0), c(0.5, 0.5)]);
    let new_pos = agent.gradient_navigate(0.01);
    println!("Navigated to: {:?}", new_pos);
}
```

## API Reference

### Module: `holomorphic`

Holomorphic functions, Cauchy-Riemann equations, power series, analytic continuation.

| Type / Function | Description |
|---|---|
| `cauchy_riemann_check(u, v, x, y, h) → bool` | Verify ∂u/∂x = ∂v/∂y and ∂u/∂y = −∂v/∂x at (x, y) via central differences |
| `complex_derivative(f, z, h) → Complex64` | Compute f'(z) numerically |
| `is_holomorphic(f, z, h) → bool` | Check Cauchy-Riemann equations at z |
| `PowerSeries` | Power series Σ aₙ(z − c)ⁿ with `center` and `coefficients` fields |
| `PowerSeries::new(center, coefficients)` | Constructor |
| `PowerSeries::evaluate(z) → Complex64` | Evaluate the series at z |
| `PowerSeries::radius_of_convergence() → Option<f64>` | Estimate via root test: 1/lim sup |aₙ|^{1/n} |
| `PowerSeries::differentiate() → PowerSeries` | Term-by-term differentiation: n·aₙ → (n+1)·a_{n+1} |
| `PowerSeries::integrate() → PowerSeries` | Term-by-term integration (constant = 0): aₙ → aₙ/(n+1) |
| `analytic_continuation(series, new_center, terms) → PowerSeries` | Re-center a power series at a new point via Taylor expansion |

### Module: `contour`

Contour integration, Cauchy integral formula, winding numbers.

| Type / Function | Description |
|---|---|
| `trait Contour` | A parameterized curve γ: [0,1] → ℂ with `evaluate(t)` and `derivative(t)` |
| `CircleContour` | Circle contour (center + radius, counterclockwise) |
| `LineSegment` | Straight line from start to end |
| `PolygonContour` | Closed polygon through vertices |
| `contour_integral(f, contour, n) → Complex64` | Numerical integration ∫_γ f(z) dz via midpoint rule with n points |
| `cauchy_integral_formula(f, a, contour, n) → Complex64` | f(a) = (1/2πi) ∮ f(z)/(z−a) dz |
| `cauchy_nth_derivative(f, a, n_order, contour, n_points) → Complex64` | f⁽ⁿ⁾(a) via Cauchy's formula |
| `cauchy_theorem_check(f, contour, n, tol) → bool` | Verify ∮ f(z) dz ≈ 0 (holomorphic f, simply connected domain) |
| `winding_number(contour, a, n) → i64` | n(γ, a) = (1/2πi) ∮ 1/(z−a) dz |

### Module: `residue`

Residue calculus, Laurent series, singularity classification, argument principle.

| Type / Function | Description |
|---|---|
| `enum SingularityType` | `Removable`, `Pole { order }`, `Essential`, `NotIsolated` |
| `residue(f, z₀, order, h) → Complex64` | Compute Res(f, z₀) for a pole of given order |
| `classify_singularity(f, z₀, h) → SingularityType` | Determine singularity type by analyzing behavior near z₀ |
| `laurent_coefficients(f, z₀, r_inner, r_outer, n_neg, n_pos) → (Vec<Complex64>, Vec<Complex64>)` | Laurent series coefficients via contour integration on an annulus |
| `argument_principle(f, f', contour, n) → i64` | Number of zeros minus poles: (1/2πi) ∮ f'(z)/f(z) dz |
| `rouche_check(f, g, contour, n) → bool` | Rouché's theorem: if |f−g| < |f| on γ, then f and g have same zero count |

### Module: `riemann_surface`

Riemann surfaces, branch cuts, covering spaces, monodromy.

| Type / Function | Description |
|---|---|
| `BranchPoint` | Location + ramification index |
| `BranchCut` | Line segment between two branch points; `is_on_cut(z, tol)` |
| `RiemannSurface` | Multivalued function structure: branch points, cuts, number of sheets |
| `RiemannSurface::sqrt_surface()` | Surface for √z (2 sheets) |
| `RiemannSurface::log_surface()` | Surface for log(z) (∞ sheets) |
| `RiemannSurface::nth_root_surface(n)` | Surface for z^{1/n} (n sheets) |
| `RiemannSurface::sqrt_on_sheet(z, sheet) → Complex64` | Evaluate √z on sheet 0 or 1 |
| `RiemannSurface::log_on_sheet(z, sheet) → Complex64` | log(z) + 2πi·sheet |
| `RiemannSurface::nth_root_on_sheet(z, n, sheet) → Complex64` | z^{1/n} on specified sheet |
| `Monodromy` | Sheet permutation when looping around a branch point |
| `Monodromy::identity(n)` | No-op permutation |
| `Monodromy::transposition(n, i, j)` | Swap sheets i and j |
| `Monodromy::compose(other)` | Compose two permutations |
| `Monodromy::apply(sheet) → usize` | Apply permutation to sheet index |
| `Monodromy::is_identity() → bool` | Check if identity |
| `euler_characteristic(n_sheets, n_branch_points) → i64` | χ = 2n − b |
| `genus(euler_char) → i64` | g = (2 − χ)/2 |

### Module: `conformal_map`

Möbius transformations, Schwarz-Christoffel, Riemann mapping.

| Type / Function | Description |
|---|---|
| `MobiusTransformation` | T(z) = (az + b)/(cz + d) with fields a, b, c, d |
| `MobiusTransformation::new(a, b, c, d)` | Constructor |
| `MobiusTransformation::identity()` | T(z) = z |
| `MobiusTransformation::from_three_points(z₁..z₃, w₁..w₃)` | Unique Möbius map sending zᵢ → wᵢ |
| `MobiusTransformation::apply(z) → Complex64` | Evaluate T(z) |
| `MobiusTransformation::determinant() → Complex64` | ad − bc |
| `MobiusTransformation::inverse() → MobiusTransformation` | T⁻¹ |
| `MobiusTransformation::compose(other) → MobiusTransformation` | T₂ ∘ T₁ |
| `MobiusTransformation::fixed_points() → Vec<Complex64>` | Solve T(z) = z |
| `MobiusTransformation::disk_to_halfplane()` | Map unit disk to upper half-plane |
| `MobiusTransformation::preserves_unit_circle() → bool` | Check if unit circle is invariant |
| `cross_ratio(z₁, z₂, z₃, z₄) → Complex64` | Möbius-invariant cross-ratio |
| `SchwarzChristoffel` | S-C mapping with prevertices, interior angles, integration constant |
| `SchwarzChristoffel::to_rectangle(a)` | Map upper half-plane to rectangle |
| `SchwarzChristoffel::integrand(z) → Complex64` | S-C integrand at z |
| `SchwarzChristoffel::evaluate(z, n_steps) → Complex64` | Numerical S-C integral from 0 to z |
| `approximate_riemann_mapping(boundary, center, n_terms) → Vec<Complex64>` | Approximate Riemann mapping coefficients |

### Module: `harmonic`

Harmonic functions, Dirichlet problem, Poisson kernel, Green's function.

| Type / Function | Description |
|---|---|
| `mean_value_property(u, x₀, y₀, r, n) → (f64, f64)` | Returns (u(center), mean on circle) |
| `is_harmonic(u, x₀, y₀, r, tol) → bool` | Check mean value property |
| `laplacian(u, x, y, h) → f64` | ∇²u via finite differences |
| `poisson_kernel(r, θ) → f64` | P(r, θ) = (1−r²)/(1−2r cos θ+r²) |
| `dirichlet_disk(boundary, r, θ, n) → f64` | Solve Dirichlet problem on unit disk via Poisson integral |
| `dirichlet_rectangle(nx, ny, dx, dy, boundary, max_iter, tol) → DMatrix<f64>` | Gauss-Seidel solver on rectangle |
| `harmonic_conjugate(u, x₀, y₀, x, y, h) → f64` | Find v such that u + iv is analytic |
| `greens_function_disk(z, ζ) → f64` | Green's function for unit disk |

### Module: `entire`

Entire functions, order and type, Liouville's theorem, Weierstrass/Hadamard factorization.

| Type / Function | Description |
|---|---|
| `compute_order(f, radii, n_θ) → f64` | ρ = lim sup (log log M(r))/(log r) |
| `check_liouville(f, radii, n_θ) → (bool, f64)` | Check if entire f is constant (bounded + low variation) |
| `compute_type(f, order, radii, n_θ) → f64` | σ = lim sup (log M(r))/r^ρ |
| `weierstrass_factor(z, p) → Complex64` | Canonical factor Eₚ(z) = (1−z)exp(z + z²/2 + ⋯ + zᵖ/p) |
| `WeierstrassProduct` | f(z) = zᵐ e^{g(z)} ∏ Eₚ(z/aₙ) |
| `WeierstrassProduct::new(origin_order, zeros, genus)` | Constructor |
| `WeierstrassProduct::evaluate(z) → Complex64` | Evaluate the product |
| `WeierstrassProduct::compute_genus(zeros) → usize` | Smallest p such that Σ 1/|aₙ|^{p+1} converges |
| `hadamard_factorization(zeros, origin_order, order) → WeierstrassProduct` | Construct Hadamard product from zeros |
| `fundamental_theorem_of_algebra(coefficients) → Vec<Complex64>` | Find all roots via Durand-Kerner |

### Module: `julia`

Julia sets, Mandelbrot set, Fatou components, iteration dynamics.

| Type / Function | Description |
|---|---|
| `IterationResult` | Orbit, escape status, escape iteration, convergence |
| `iterate_quadratic(z₀, c, max_iter, escape_r) → IterationResult` | Iterate f(z) = z² + c |
| `mandelbrot_membership(c, max_iter, escape_r) → (bool, Option<usize>)` | Is c in the Mandelbrot set? |
| `julia_membership(z₀, c, max_iter, escape_r) → (bool, Option<usize>)` | Is z₀ in the Julia set for parameter c? |
| `filled_julia_set(c, x_range, y_range, res, max_iter, escape_r) → Vec<Vec<usize>>` | Grid of escape iteration counts |
| `mandelbrot_set(x_range, y_range, res, max_iter, escape_r) → Vec<Vec<usize>>` | Grid of Mandelbrot escape counts |
| `enum FatouComponent` | AttractingFixedPoint, PeriodicCycle, BasinOfInfinity, SiegelDisk, Unknown |
| `quadratic_fixed_points(c) → (Complex64, Complex64)` | Solve z = z² + c |
| `find_periodic_points(c, period, n_initial) → Vec<Complex64>` | Find periodic cycles by sampling |
| `multiplier(f, f', fixed_point) → Complex64` | λ = f'(z*) |
| `enum FixedPointType` | SuperAttracting, Attracting, Repelling, Indifferent, Parabolic |
| `classify_fixed_point(λ) → FixedPointType` | Classify by multiplier magnitude |
| `mandelbrot_distance(c, max_iter, escape_r) → f64` | Distance estimation to Mandelbrot boundary |

### Module: `potential`

Complex potential theory, logarithmic potential, equilibrium measures, Fekete points.

| Type / Function | Description |
|---|---|
| `logarithmic_potential(z, atoms) → f64` | U^μ(z) = Σ wᵢ log|z − tᵢ| |
| `potential_gradient(z, atoms) → Complex64` | ∇U^μ(z) |
| `DiscreteMeasure` | Discrete probability measure with atoms (point, weight) |
| `DiscreteMeasure::new(atoms)` | Constructor |
| `DiscreteMeasure::uniform(points)` | Equal-weight measure |
| `DiscreteMeasure::total_mass() → f64` | Σ wᵢ |
| `DiscreteMeasure::energy() → f64` | I(μ) = −∫∫ log|z−t| dμ(z) dμ(t) |
| `DiscreteMeasure::potential(z) → f64` | U^μ(z) |
| `DiscreteMeasure::robin_constant() → f64` | −I(μ) |
| `DiscreteMeasure::capacity() → f64` | exp(−I(μ)) (transfinite diameter) |
| `unit_circle_equilibrium(n) → DiscreteMeasure` | Uniform measure on unit circle |
| `chebyshev_nodes(n) → Vec<Complex64>` | Chebyshev nodes on [−1, 1] |
| `approximate_fekete_points(n, iterations, lr) → Vec<Complex64>` | Energy-minimizing points via gradient descent |
| `exterior_green_function(z) → f64` | max(log|z|, 0) |
| `dirichlet_energy(u, x_range, y_range, nx, ny, h) → f64` | ∫|∇u|² dA |
| `superharmonic_check(u, v, points) → bool` | Check u ≥ v on test set |

### Module: `agent_complex` (unified API)

| Type / Function | Description |
|---|---|
| `enum DecisionBoundary` | JuliaSet, MandelbrotSet, Custom |
| `enum ReasoningMode` | Holomorphic, Contour, Residue, Potential, Dynamics, Conformal |
| `AgentState` | position, velocity, potential, iteration_count, on_boundary, basin_label |
| `AgentComplex` | Main type — complex-analytic reasoning agent |
| `AgentComplex::with_julia_boundary(c, position)` | Create agent with Julia set boundary |
| `AgentComplex::with_mandelbrot_boundary(initial_c)` | Create agent with Mandelbrot boundary |
| `AgentComplex::classify_position(max_iter) → BoundaryClassification` | Inside/outside decision boundary + escape time |
| `AgentComplex::compute_orbit(max_iter) → Vec<Complex64>` | Full orbit from current position |
| `AgentComplex::analyze_holomorphic(f) → HolomorphicAnalysis` | Check holomorphicity + derivative at position |
| `AgentComplex::contour_integrate(f, radius, n) → Complex64` | Integral around position |
| `AgentComplex::compute_residue(f, order) → Complex64` | Residue at position |
| `AgentComplex::classify_singularity(f) → SingularityType` | Singularity type at position |
| `AgentComplex::apply_conformal_map() → Complex64` | Transform position via Möbius map |
| `AgentComplex::setup_mobius_map(a, b, c, d)` | Configure conformal map |
| `AgentComplex::compute_potential() → f64` | Logarithmic potential at position |
| `AgentComplex::setup_equilibrium_measure(points)` | Set up discrete measure |
| `AgentComplex::gradient_navigate(step) → Complex64` | Move along potential gradient |
| `AgentComplex::check_harmonic_region(u, radius) → bool` | Is position in a harmonic region? |
| `AgentComplex::analyze_dynamics(c) → DynamicsAnalysis` | Fixed points + multipliers for z²+c |
| `AgentComplex::riemann_structure(n) → RiemannSurface` | Riemann surface with n sheets |
| `AgentComplex::reason(f, max_iter) → ReasoningReport` | Full reasoning: classification + holomorphicity + singularity + residue + potential |
| `BoundaryClassification` | in_set, escape_iteration, boundary_type |
| `HolomorphicAnalysis` | is_holomorphic, derivative, value |
| `DynamicsAnalysis` | fixed_points, multipliers, fixed_point_types |
| `ReasoningReport` | Full report: position + classification + holomorphic + singularity + residue + potential |

## How It Works

The library is layered from foundational to applied:

1. **`holomorphic`** — The analytic primitives. Power series, derivatives, Cauchy-Riemann checks. Everything rests on the fact that holomorphic functions are rigid: knowing them on a tiny disk determines them everywhere.

2. **`contour`** — Integration along curves. The Cauchy integral formula converts global information (integral around a loop) into local information (value at a point). This is the engine that powers residue calculus.

3. **`residue`** — Singularities carry information. The residue is the coefficient of 1/(z−z₀) in the Laurent expansion, and the residue theorem converts a sum of residues into a contour integral (and vice versa).

4. **`riemann_surface`** — Multivalued functions (√z, log z, z^{1/n}) live on covering spaces. Branch cuts choose one value. Monodromy tracks which sheet you're on as you loop around branch points.

5. **`conformal_map`** — Möbius transformations are the automorphisms of the Riemann sphere. They map circles to circles, preserve angles, and provide the simplest nontrivial conformal equivalences. Schwarz-Christoffel extends this to polygonal domains.

6. **`harmonic`** — The real and imaginary parts of holomorphic functions are harmonic (Laplacian = 0). The Dirichlet problem (find u with given boundary values and ∇²u = 0) is solved by the Poisson kernel on the disk and by Gauss-Seidel iteration on rectangles.

7. **`entire`** — Functions holomorphic on all of ℂ. Their growth is classified by order ρ and type σ. Weierstrass factorization reconstructs them from their zeros. Liouville's theorem says bounded entire functions are constant.

8. **`julia`** — Iterating f(z) = z² + c partitions the plane: points that stay bounded form the filled Julia set, points that escape form basins of infinity. The Mandelbrot set is the set of c values where the Julia set is connected.

9. **`potential`** — Logarithmic potentials measure how much "charge" is nearby. Equilibrium measures minimize energy. Fekete points are the optimal discrete approximation. Capacity quantifies how "large" a set is from the potential theory perspective.

10. **`agent_complex`** — Unifies all modules into an agent that lives on the complex plane with a decision boundary, uses holomorphic analysis to understand its local geometry, contour integration for global properties, potential theory for navigation, and dynamical systems for boundary classification.

## The Math

### Cauchy-Riemann Equations
A function f(z) = u(x,y) + iv(x,y) is holomorphic if ∂u/∂x = ∂v/∂y and ∂u/∂y = −∂v/∂x. These equations guarantee that f has a complex derivative f'(z) that is independent of direction.

### Cauchy Integral Formula
For holomorphic f inside a contour γ: **f(a) = (1/2πi) ∮_γ f(z)/(z−a) dz**. This means the value at any interior point is determined by values on the boundary. The formula generalizes to derivatives: **f⁽ⁿ⁾(a) = n!/(2πi) ∮ f(z)/(z−a)^{n+1} dz**.

### Residue Theorem
If f has isolated singularities z₁, …, zₖ inside γ: **∮_γ f(z) dz = 2πi Σ Res(f, zₖ)**. This converts contour integrals into algebraic computations.

### Winding Number
n(γ, a) = (1/2πi) ∮_γ 1/(z−a) dz counts how many times γ wraps around a. It's the topological degree of the map γ → S¹.

### Laurent Series
Around an isolated singularity: **f(z) = Σ_{n=−∞}^{∞} aₙ(z−z₀)ⁿ**. The coefficient a₋₁ is the residue. Coefficients are computed via: **aₙ = (1/2πi) ∮ f(z)/(z−z₀)^{n+1} dz** on an annulus.

### Power Series
A holomorphic function has a Taylor expansion **f(z) = Σ aₙ(z−c)ⁿ** with radius of convergence R = 1/lim sup |aₙ|^{1/n}. Differentiation and integration work term-by-term.

### Riemann Surfaces
Multivalued functions like √z require "sheets" — copies of ℂ glued along branch cuts. The monodromy group describes how sheets permute when you loop around branch points. The Euler characteristic χ = 2n − b gives the topology (genus g = (2−χ)/2).

### Möbius Transformations
T(z) = (az+b)/(cz+d) are the conformal automorphisms of the Riemann sphere. They form a group under composition, are determined by 3 points, and preserve cross-ratios. The determinant ad−bc must be nonzero.

### Schwarz-Christoffel
Maps the upper half-plane to a polygon with interior angles παₖ via the integral **f(z) = C ∫ ∏(z−xₖ)^{αₖ−1} dz**. This is the canonical way to construct conformal maps to polygonal regions.

### Poisson Kernel
Solves the Dirichlet problem on the unit disk: **u(r, θ) = (1/2π) ∫₀^{2π} P(r, θ−φ) f(φ) dφ** where P(r, θ) = (1−r²)/(1−2r cos θ+r²). The kernel is positive and integrates to 1.

### Entire Functions
Order: ρ = lim sup (log log M(r))/(log r). Type: σ = lim sup (log M(r))/r^ρ. Weierstrass factorization: f(z) = zᵐ e^{g(z)} ∏ Eₚ(z/aₙ) where Eₚ are canonical factors ensuring convergence.

### Julia/Mandelbrot Dynamics
For f(z) = z² + c: the Julia set J_c is the boundary between bounded and escaping orbits. The Mandelbrot set M = {c : J_c is connected}. Fixed points satisfy z² − z + c = 0. The multiplier λ = 2z* classifies: |λ| < 1 attracting, |λ| > 1 repelling, |λ| = 1 indifferent.

### Potential Theory
The logarithmic potential U^μ(z) = ∫ log|z−t| dμ(t) is subharmonic. The energy I(μ) = −∫∫ log|z−t| dμ(z)dμ(t) is minimized by the equilibrium measure. The Robin constant V(K) = −inf_μ I(μ), and capacity cap(K) = e^{−V(K)}.

## License

MIT
