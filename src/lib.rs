//! # lau-complex-agents
//!
//! Complex analysis and holomorphic dynamics for agents —
//! contour integration, residue theory, Riemann surfaces, and Julia sets as agent decision boundaries.

pub mod holomorphic;
pub mod contour;
pub mod residue;
pub mod riemann_surface;
pub mod conformal_map;
pub mod harmonic;
pub mod entire;
pub mod julia;
pub mod potential;
pub mod agent_complex;

pub use agent_complex::{AgentComplex, ReasoningReport, BoundaryClassification};
