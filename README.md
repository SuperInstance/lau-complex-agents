# lau-complex-agents

> Complex analysis and holomorphic dynamics for agents

## What This Does

Complex analysis and holomorphic dynamics for agents. Part of the PLATO/LAU ecosystem — a mathematically rigorous framework for building educational agents that learn, teach, and evolve.

## The Key Idea

This crate implements the core abstractions needed for its domain, with a focus on correctness, composability, and conservation guarantees. Every public type is serializable (serde), every algorithm is tested, and every invariant is verified.

## Install

```bash
cargo add lau-complex-agents
```

## Quick Start

See the API Reference below for complete usage. Key entry points:

```rust
use lau_complex_agents::*;
// See types and methods below for complete usage
```

## API Reference

```rust
pub enum DecisionBoundary 
pub enum ReasoningMode 
pub struct AgentState 
pub struct AgentComplex 
    pub fn with_julia_boundary(c: Complex64, initial_position: Complex64) -> Self 
    pub fn with_mandelbrot_boundary(initial_c: Complex64) -> Self 
    pub fn classify_position(&mut self, max_iterations: usize) -> BoundaryClassification 
    pub fn compute_orbit(&self, max_iterations: usize) -> Vec<Complex64> 
    pub fn analyze_holomorphic(
    pub fn contour_integrate(
    pub fn compute_residue(
    pub fn classify_singularity(
    pub fn apply_conformal_map(&mut self) -> Complex64 
    pub fn setup_mobius_map(&mut self, a: Complex64, b: Complex64, c: Complex64, d: Complex64) 
    pub fn compute_potential(&mut self) -> f64 
    pub fn setup_equilibrium_measure(&mut self, points: Vec<Complex64>) 
    pub fn gradient_navigate(&mut self, step_size: f64) -> Complex64 
    pub fn check_harmonic_region(
    pub fn analyze_dynamics(&self, c: Complex64) -> DynamicsAnalysis 
    pub fn riemann_structure(&self, n_branches: usize) -> RiemannSurface 
    pub fn reason(
pub struct BoundaryClassification 
pub struct HolomorphicAnalysis 
pub struct DynamicsAnalysis 
pub struct ReasoningReport 
pub struct MobiusTransformation 
    pub fn new(a: Complex64, b: Complex64, c: Complex64, d: Complex64) -> Self 
    pub fn identity() -> Self 
    pub fn from_three_points(z1: Complex64, z2: Complex64, z3: Complex64, w1: Complex64, w2: Complex64, w3: Complex64) -> Self 
    pub fn apply(&self, z: Complex64) -> Complex64 
    pub fn determinant(&self) -> Complex64 
    pub fn inverse(&self) -> Self 
    pub fn compose(&self, other: &MobiusTransformation) -> Self 
    pub fn fixed_points(&self) -> Vec<Complex64> 
    pub fn disk_to_halfplane() -> Self 
    pub fn preserves_unit_circle(&self) -> bool 
pub fn cross_ratio(z1: Complex64, z2: Complex64, z3: Complex64, z4: Complex64) -> Complex64 
pub struct SchwarzChristoffel 
    pub fn new(prevertices: Vec<f64>, alpha: Vec<f64>, c: Complex64) -> Self 
    pub fn to_rectangle(a: f64) -> Self 
    pub fn integrand(&self, z: Complex64) -> Complex64 
    pub fn evaluate(&self, z: Complex64, n_steps: usize) -> Complex64 
pub fn approximate_riemann_mapping(
pub fn compute_order(
pub fn check_liouville(
pub fn compute_type(
pub fn weierstrass_factor(z: Complex64, p: usize) -> Complex64 
pub struct WeierstrassProduct 
    pub fn new(origin_order: usize, zeros: Vec<Complex64>, genus: usize) -> Self 
    pub fn evaluate(&self, z: Complex64) -> Complex64 
    pub fn compute_genus(zeros: &[Complex64]) -> usize 
pub fn hadamard_factorization(
pub fn fundamental_theorem_of_algebra(coefficients: &[f64]) -> Vec<Complex64> 
pub fn logarithmic_potential(z: Complex64, atoms: &[(Complex64, f64)]) -> f64 
pub fn potential_gradient(z: Complex64, atoms: &[(Complex64, f64)]) -> Complex64 
pub struct DiscreteMeasure 
    pub fn new(atoms: Vec<(Complex64, f64)>) -> Self 
    pub fn uniform(points: Vec<Complex64>) -> Self 
    pub fn total_mass(&self) -> f64 
    pub fn energy(&self) -> f64 
```

## How It Works

Read the source in `src/` for full implementation details. All algorithms are documented with inline comments explaining the mathematical foundations.

## The Math

This crate implements formal mathematical constructs. See the source documentation for theorem statements and proofs of correctness.

## Testing

**113 tests** covering construction, serialization, correctness properties, edge cases, and composability with other lau-* crates.

## License

MIT
