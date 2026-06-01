//! Unified API — AgentComplex with complex-analytic reasoning

use num_complex::Complex64;
use serde::{Deserialize, Serialize};

use crate::holomorphic::{self, PowerSeries};
use crate::contour::{self, CircleContour, Contour};
use crate::residue::{self, SingularityType};
use crate::riemann_surface::{self, RiemannSurface, Monodromy};
use crate::conformal_map::{self, MobiusTransformation};
use crate::harmonic;
use crate::entire::{self, WeierstrassProduct};
use crate::julia;
use crate::potential::{self, DiscreteMeasure};

/// Decision boundary classification for agents.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum DecisionBoundary {
    JuliaSet { c: Complex64, resolution: usize },
    MandelbrotSet { resolution: usize },
    Custom { name: String },
}

/// Agent reasoning mode.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ReasoningMode {
    /// Analyze via holomorphic function theory.
    Holomorphic,
    /// Contour integration based reasoning.
    Contour,
    /// Residue calculus.
    Residue,
    /// Potential theory.
    Potential,
    /// Dynamical systems (Julia/Mandelbrot).
    Dynamics,
    /// Conformal mapping.
    Conformal,
}

/// Complex-analytic agent state.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AgentState {
    pub position: Complex64,
    pub velocity: Complex64,
    pub potential: f64,
    pub iteration_count: usize,
    pub on_decision_boundary: bool,
    pub basin_label: Option<usize>,
}

/// The main AgentComplex type — provides complex-analytic reasoning for agents.
pub struct AgentComplex {
    /// Decision boundary for the agent.
    pub boundary: DecisionBoundary,
    /// Current state.
    pub state: AgentState,
    /// Reasoning mode.
    pub mode: ReasoningMode,
    /// Conformal map for domain transformation.
    pub conformal_map: Option<MobiusTransformation>,
    /// Equilibrium measure for potential theory.
    pub measure: Option<DiscreteMeasure>,
}

impl AgentComplex {
    /// Create a new agent with Julia set decision boundary.
    pub fn with_julia_boundary(c: Complex64, initial_position: Complex64) -> Self {
        Self {
            boundary: DecisionBoundary::JuliaSet { c, resolution: 100 },
            state: AgentState {
                position: initial_position,
                velocity: Complex64::new(0.0, 0.0),
                potential: 0.0,
                iteration_count: 0,
                on_decision_boundary: false,
                basin_label: None,
            },
            mode: ReasoningMode::Dynamics,
            conformal_map: None,
            measure: None,
        }
    }

    /// Create a new agent with Mandelbrot decision boundary.
    pub fn with_mandelbrot_boundary(initial_c: Complex64) -> Self {
        Self {
            boundary: DecisionBoundary::MandelbrotSet { resolution: 100 },
            state: AgentState {
                position: initial_c,
                velocity: Complex64::new(0.0, 0.0),
                potential: 0.0,
                iteration_count: 0,
                on_decision_boundary: false,
                basin_label: None,
            },
            mode: ReasoningMode::Dynamics,
            conformal_map: None,
            measure: None,
        }
    }

    /// Classify the current position relative to the decision boundary.
    pub fn classify_position(&mut self, max_iterations: usize) -> BoundaryClassification {
        match &self.boundary {
            DecisionBoundary::JuliaSet { c, .. } => {
                let (in_set, escape_iter) = julia::julia_membership(
                    self.state.position, *c, max_iterations, 10.0,
                );
                self.state.on_decision_boundary = in_set;
                self.state.basin_label = if in_set { Some(0) } else { escape_iter };
                BoundaryClassification {
                    in_set,
                    escape_iteration: escape_iter,
                    boundary_type: "Julia".to_string(),
                }
            }
            DecisionBoundary::MandelbrotSet { .. } => {
                let (in_set, escape_iter) = julia::mandelbrot_membership(
                    self.state.position, max_iterations, 10.0,
                );
                self.state.on_decision_boundary = in_set;
                BoundaryClassification {
                    in_set,
                    escape_iteration: escape_iter,
                    boundary_type: "Mandelbrot".to_string(),
                }
            }
            DecisionBoundary::Custom { name } => {
                BoundaryClassification {
                    in_set: false,
                    escape_iteration: None,
                    boundary_type: name.clone(),
                }
            }
        }
    }

    /// Compute the orbit from current position.
    pub fn compute_orbit(&self, max_iterations: usize) -> Vec<Complex64> {
        match &self.boundary {
            DecisionBoundary::JuliaSet { c, .. } => {
                let result = julia::iterate_quadratic(
                    self.state.position, *c, max_iterations, 1e10,
                );
                result.orbit
            }
            DecisionBoundary::MandelbrotSet { .. } => {
                let result = julia::iterate_quadratic(
                    Complex64::new(0.0, 0.0), self.state.position, max_iterations, 1e10,
                );
                result.orbit
            }
            DecisionBoundary::Custom { .. } => vec![self.state.position],
        }
    }

    /// Analyze the decision landscape using holomorphic function theory.
    pub fn analyze_holomorphic(
        &self,
        f: &dyn Fn(Complex64) -> Complex64,
    ) -> HolomorphicAnalysis {
        let z = self.state.position;
        HolomorphicAnalysis {
            is_holomorphic: holomorphic::is_holomorphic(f, z, 1e-7),
            derivative: holomorphic::complex_derivative(f, z, 1e-8),
            value: f(z),
        }
    }

    /// Evaluate a contour integral around the current position.
    pub fn contour_integrate(
        &self,
        f: &dyn Fn(Complex64) -> Complex64,
        radius: f64,
        n_points: usize,
    ) -> Complex64 {
        let circle = CircleContour::new(self.state.position, radius);
        contour::contour_integral(f, &circle, n_points)
    }

    /// Compute residue at current position.
    pub fn compute_residue(
        &self,
        f: &dyn Fn(Complex64) -> Complex64,
        pole_order: usize,
    ) -> Complex64 {
        residue::residue(f, self.state.position, pole_order, 1e-6)
    }

    /// Classify singularity at current position.
    pub fn classify_singularity(
        &self,
        f: &dyn Fn(Complex64) -> Complex64,
    ) -> SingularityType {
        residue::classify_singularity(f, self.state.position, 1e-6)
    }

    /// Apply conformal map to move to canonical domain.
    pub fn apply_conformal_map(&mut self) -> Complex64 {
        if let Some(ref map) = self.conformal_map {
            let new_pos = map.apply(self.state.position);
            self.state.position = new_pos;
            new_pos
        } else {
            self.state.position
        }
    }

    /// Set up a Möbius transformation mapping current region.
    pub fn setup_mobius_map(&mut self, a: Complex64, b: Complex64, c: Complex64, d: Complex64) {
        self.conformal_map = Some(MobiusTransformation::new(a, b, c, d));
    }

    /// Compute potential at current position.
    pub fn compute_potential(&mut self) -> f64 {
        if let Some(ref measure) = self.measure {
            let pot = measure.potential(self.state.position);
            self.state.potential = pot;
            pot
        } else {
            0.0
        }
    }

    /// Set up equilibrium measure from a set of points.
    pub fn setup_equilibrium_measure(&mut self, points: Vec<Complex64>) {
        self.measure = Some(DiscreteMeasure::uniform(points));
    }

    /// Navigate using gradient of the potential field.
    pub fn gradient_navigate(&mut self, step_size: f64) -> Complex64 {
        if let Some(ref measure) = self.measure {
            let grad = potential::potential_gradient(self.state.position, &measure.atoms);
            let new_pos = self.state.position + step_size * grad;
            self.state.velocity = step_size * grad;
            self.state.position = new_pos;
            new_pos
        } else {
            self.state.position
        }
    }

    /// Check if current position is in a harmonic region.
    pub fn check_harmonic_region(
        &self,
        u: &dyn Fn(f64, f64) -> f64,
        radius: f64,
    ) -> bool {
        harmonic::is_harmonic(u, self.state.position.re, self.state.position.im, radius, 1e-3)
    }

    /// Analyze decision dynamics using fixed point analysis.
    pub fn analyze_dynamics(&self, c: Complex64) -> DynamicsAnalysis {
        let (fp1, fp2) = julia::quadratic_fixed_points(c);
        let f_prime = |z: Complex64| 2.0 * z;
        let lambda1 = f_prime(fp1);
        let lambda2 = f_prime(fp2);

        DynamicsAnalysis {
            fixed_points: vec![fp1, fp2],
            multipliers: vec![lambda1, lambda2],
            fixed_point_types: vec![
                julia::classify_fixed_point(lambda1),
                julia::classify_fixed_point(lambda2),
            ],
        }
    }

    /// Compute the Riemann surface structure for the decision function.
    pub fn riemann_structure(&self, n_branches: usize) -> RiemannSurface {
        RiemannSurface::nth_root_surface(n_branches)
    }

    /// Full complex-analytic reasoning report.
    pub fn reason(
        &mut self,
        f: &dyn Fn(Complex64) -> Complex64,
        max_iterations: usize,
    ) -> ReasoningReport {
        let classification = self.classify_position(max_iterations);
        let holomorphic = self.analyze_holomorphic(f);
        let singularity = self.classify_singularity(f);
        let residue_val = self.compute_residue(f, 1);

        ReasoningReport {
            position: self.state.position,
            classification,
            holomorphic,
            singularity_type: format!("{:?}", singularity),
            residue: residue_val,
            potential: self.state.potential,
        }
    }
}

/// Boundary classification result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BoundaryClassification {
    pub in_set: bool,
    pub escape_iteration: Option<usize>,
    pub boundary_type: String,
}

/// Holomorphic analysis result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HolomorphicAnalysis {
    pub is_holomorphic: bool,
    pub derivative: Complex64,
    pub value: Complex64,
}

/// Dynamics analysis result.
#[derive(Debug, Clone)]
pub struct DynamicsAnalysis {
    pub fixed_points: Vec<Complex64>,
    pub multipliers: Vec<Complex64>,
    pub fixed_point_types: Vec<julia::FixedPointType>,
}

/// Full reasoning report.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReasoningReport {
    pub position: Complex64,
    pub classification: BoundaryClassification,
    pub holomorphic: HolomorphicAnalysis,
    pub singularity_type: String,
    pub residue: Complex64,
    pub potential: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn c(re: f64, im: f64) -> Complex64 {
        Complex64::new(re, im)
    }

    #[test]
    fn test_agent_julia_creation() {
        let agent = AgentComplex::with_julia_boundary(c(-0.7, 0.27015), c(0.0, 0.0));
        assert_eq!(agent.state.position, c(0.0, 0.0));
    }

    #[test]
    fn test_agent_classify_inside_julia() {
        let mut agent = AgentComplex::with_julia_boundary(c(0.0, 0.0), c(0.0, 0.0));
        let result = agent.classify_position(100);
        assert!(result.in_set); // origin is in Julia set for c=0
    }

    #[test]
    fn test_agent_classify_outside_julia() {
        let mut agent = AgentComplex::with_julia_boundary(c(0.0, 0.0), c(3.0, 0.0));
        let result = agent.classify_position(100);
        assert!(!result.in_set);
    }

    #[test]
    fn test_agent_mandelbrot_inside() {
        let mut agent = AgentComplex::with_mandelbrot_boundary(c(0.0, 0.0));
        let result = agent.classify_position(100);
        assert!(result.in_set);
    }

    #[test]
    fn test_agent_orbit() {
        let agent = AgentComplex::with_julia_boundary(c(0.0, 0.0), c(1.0, 0.0));
        let orbit = agent.compute_orbit(5);
        // 1 -> 1 -> 1 -> ... (fixed point)
        assert_eq!(orbit[0], c(1.0, 0.0));
    }

    #[test]
    fn test_agent_holomorphic_analysis() {
        let agent = AgentComplex::with_julia_boundary(c(0.0, 0.0), c(1.0, 0.0));
        let analysis = agent.analyze_holomorphic(&|z| z.exp());
        assert!(analysis.is_holomorphic);
    }

    #[test]
    fn test_agent_contour_integrate() {
        let agent = AgentComplex::with_julia_boundary(c(0.0, 0.0), c(0.0, 0.0));
        let integral = agent.contour_integrate(&|z| z.exp(), 1.0, 1000);
        assert!(integral.norm() < 0.01); // holomorphic => ~0
    }

    #[test]
    fn test_agent_conformal_map() {
        let mut agent = AgentComplex::with_julia_boundary(c(0.0, 0.0), c(1.0, 0.0));
        agent.setup_mobius_map(c(2.0, 0.0), c(0.0, 0.0), c(0.0, 0.0), c(1.0, 0.0));
        let new_pos = agent.apply_conformal_map();
        assert!((new_pos - c(2.0, 0.0)).norm() < 1e-10);
    }

    #[test]
    fn test_agent_potential() {
        let mut agent = AgentComplex::with_julia_boundary(c(0.0, 0.0), c(0.5, 0.0));
        agent.setup_equilibrium_measure(vec![c(0.0, 0.0), c(1.0, 0.0)]);
        let pot = agent.compute_potential();
        assert!(pot.is_finite());
    }

    #[test]
    fn test_agent_dynamics_analysis() {
        let agent = AgentComplex::with_julia_boundary(c(0.0, 0.0), c(0.0, 0.0));
        let analysis = agent.analyze_dynamics(c(0.0, 0.0));
        assert_eq!(analysis.fixed_points.len(), 2);
        assert_eq!(analysis.fixed_point_types.len(), 2);
    }

    #[test]
    fn test_agent_riemann_structure() {
        let agent = AgentComplex::with_julia_boundary(c(0.0, 0.0), c(0.0, 0.0));
        let surface = agent.riemann_structure(3);
        assert_eq!(surface.n_sheets, 3);
    }

    #[test]
    fn test_agent_full_reasoning() {
        let mut agent = AgentComplex::with_julia_boundary(c(-0.5, 0.0), c(0.0, 0.0));
        agent.setup_equilibrium_measure(vec![c(0.0, 0.0), c(1.0, 0.0)]);
        let report = agent.reason(&|z| z.exp(), 100);
        assert!(report.holomorphic.is_holomorphic);
    }

    #[test]
    fn test_agent_residue() {
        let mut agent = AgentComplex::with_julia_boundary(c(0.0, 0.0), c(0.0, 0.0));
        // 1/z at origin
        let res = agent.compute_residue(&|z| 1.0 / z, 1);
        assert!((res - c(1.0, 0.0)).norm() < 0.5);
    }

    #[test]
    fn test_agent_harmonic_check() {
        let agent = AgentComplex::with_julia_boundary(c(0.0, 0.0), c(1.0, 0.0));
        let u = |x: f64, y: f64| x * x - y * y;
        assert!(agent.check_harmonic_region(&u, 0.1));
    }
}
