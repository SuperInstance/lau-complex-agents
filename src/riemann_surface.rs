//! Riemann surfaces, branch cuts, covering spaces, monodromy

use num_complex::Complex64;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

/// Branch point of a multivalued function.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BranchPoint {
    pub location: Complex64,
    pub ramification_index: usize,
}

impl BranchPoint {
    pub fn new(location: Complex64, ramification_index: usize) -> Self {
        Self { location, ramification_index }
    }
}

/// A branch cut connecting two branch points.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BranchCut {
    pub start: Complex64,
    pub end: Complex64,
}

impl BranchCut {
    pub fn new(start: Complex64, end: Complex64) -> Self {
        Self { start, end }
    }

    /// Check if a point is on the branch cut (within tolerance).
    pub fn is_on_cut(&self, z: Complex64, tol: f64) -> bool {
        let d1 = (z - self.start).norm();
        let d2 = (z - self.end).norm();
        let length = (self.end - self.start).norm();
        d1 + d2 < length + tol
    }
}

/// Sheet of a Riemann surface.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Sheet {
    pub index: usize,
    pub label: String,
}

/// Riemann surface for a multivalued function.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RiemannSurface {
    pub branch_points: Vec<BranchPoint>,
    pub branch_cuts: Vec<BranchCut>,
    pub n_sheets: usize,
}

impl RiemannSurface {
    pub fn new(branch_points: Vec<BranchPoint>, branch_cuts: Vec<BranchCut>, n_sheets: usize) -> Self {
        Self { branch_points, branch_cuts, n_sheets }
    }

    /// The Riemann surface for sqrt(z) — 2 sheets, branch point at 0 and ∞.
    pub fn sqrt_surface() -> Self {
        Self {
            branch_points: vec![BranchPoint::new(Complex64::new(0.0, 0.0), 2)],
            branch_cuts: vec![BranchCut::new(Complex64::new(0.0, 0.0), Complex64::new(f64::NEG_INFINITY, 0.0))],
            n_sheets: 2,
        }
    }

    /// The Riemann surface for log(z) — ∞ sheets, branch point at 0.
    pub fn log_surface() -> Self {
        Self {
            branch_points: vec![BranchPoint::new(Complex64::new(0.0, 0.0), 0)], // 0 = infinite ramification
            branch_cuts: vec![BranchCut::new(Complex64::new(0.0, 0.0), Complex64::new(f64::NEG_INFINITY, 0.0))],
            n_sheets: usize::MAX,
        }
    }

    /// The Riemann surface for z^(1/n) — n sheets.
    pub fn nth_root_surface(n: usize) -> Self {
        Self {
            branch_points: vec![BranchPoint::new(Complex64::new(0.0, 0.0), n)],
            branch_cuts: vec![BranchCut::new(Complex64::new(0.0, 0.0), Complex64::new(f64::NEG_INFINITY, 0.0))],
            n_sheets: n,
        }
    }

    /// Compute principal value of sqrt(z) on specified sheet.
    pub fn sqrt_on_sheet(z: Complex64, sheet: usize) -> Complex64 {
        let principal = z.sqrt();
        if sheet == 0 { principal } else { -principal }
    }

    /// Compute log(z) on specified sheet.
    pub fn log_on_sheet(z: Complex64, sheet: usize) -> Complex64 {
        let principal = z.ln();
        principal + Complex64::new(0.0, 2.0 * PI * sheet as f64)
    }

    /// Compute z^(1/n) on specified sheet.
    pub fn nth_root_on_sheet(z: Complex64, n: usize, sheet: usize) -> Complex64 {
        let r = z.norm();
        let theta = z.arg();
        let new_r = r.powf(1.0 / n as f64);
        let new_theta = (theta + 2.0 * PI * sheet as f64) / n as f64;
        Complex64::from_polar(new_r, new_theta)
    }
}

/// Monodromy: track sheet transitions when looping around branch points.
#[derive(Debug, Clone)]
pub struct Monodromy {
    /// Permutation of sheets when looping around a branch point.
    /// permutation[i] = j means sheet i transitions to sheet j.
    pub permutation: Vec<usize>,
}

impl Monodromy {
    pub fn new(permutation: Vec<usize>) -> Self {
        Self { permutation }
    }

    /// Identity monodromy (no sheet transition).
    pub fn identity(n: usize) -> Self {
        Self { permutation: (0..n).collect() }
    }

    /// Simple transposition swapping sheets i and j.
    pub fn transposition(n: usize, i: usize, j: usize) -> Self {
        let mut perm: Vec<usize> = (0..n).collect();
        perm.swap(i, j);
        Self { permutation: perm }
    }

    /// Compose two monodromy permutations.
    pub fn compose(&self, other: &Monodromy) -> Monodromy {
        let n = self.permutation.len();
        let result: Vec<usize> = (0..n).map(|i| other.permutation[self.permutation[i]]).collect();
        Monodromy::new(result)
    }

    /// Apply monodromy to a sheet index.
    pub fn apply(&self, sheet: usize) -> usize {
        self.permutation[sheet]
    }

    /// Check if this is the identity permutation.
    pub fn is_identity(&self) -> bool {
        self.permutation.iter().enumerate().all(|(i, &p)| i == p)
    }
}

/// Covering space representation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CoveringSpace {
    pub base_points: Vec<Complex64>,
    /// For each base point, the fiber (points above it).
    pub fibers: Vec<Vec<Complex64>>,
    pub projection: Option<String>, // Label for the projection function
}

/// Euler characteristic of a Riemann surface.
/// For an n-sheeted covering of the sphere with b branch points:
/// χ = 2n - b
pub fn euler_characteristic(n_sheets: usize, n_branch_points: usize) -> i64 {
    2 * n_sheets as i64 - n_branch_points as i64
}

/// Genus from Euler characteristic: g = (2 - χ) / 2
pub fn genus(euler_char: i64) -> i64 {
    (2 - euler_char) / 2
}

#[cfg(test)]
mod tests {
    use super::*;

    fn c(re: f64, im: f64) -> Complex64 {
        Complex64::new(re, im)
    }

    #[test]
    fn test_sqrt_surface() {
        let surface = RiemannSurface::sqrt_surface();
        assert_eq!(surface.n_sheets, 2);
        assert_eq!(surface.branch_points.len(), 1);
    }

    #[test]
    fn test_sqrt_on_sheet_0() {
        let z = c(4.0, 0.0);
        let result = RiemannSurface::sqrt_on_sheet(z, 0);
        assert!((result - c(2.0, 0.0)).norm() < 1e-10);
    }

    #[test]
    fn test_sqrt_on_sheet_1() {
        let z = c(4.0, 0.0);
        let result = RiemannSurface::sqrt_on_sheet(z, 1);
        assert!((result - c(-2.0, 0.0)).norm() < 1e-10);
    }

    #[test]
    fn test_log_on_sheet_0() {
        let z = c(1.0, 0.0);
        let result = RiemannSurface::log_on_sheet(z, 0);
        assert!(result.re.abs() < 1e-10);
        assert!(result.im.abs() < 1e-10);
    }

    #[test]
    fn test_log_on_sheet_1() {
        let z = c(1.0, 0.0);
        let result = RiemannSurface::log_on_sheet(z, 1);
        assert!((result.im - 2.0 * PI).abs() < 1e-10);
    }

    #[test]
    fn test_nth_root_on_sheet() {
        let z = c(1.0, 0.0);
        // z^(1/3) on sheet 0 = 1
        let r0 = RiemannSurface::nth_root_on_sheet(z, 3, 0);
        assert!((r0 - c(1.0, 0.0)).norm() < 1e-10);
        // z^(1/3) on sheet 1 = e^(2πi/3)
        let r1 = RiemannSurface::nth_root_on_sheet(z, 3, 1);
        let expected = Complex64::from_polar(1.0, 2.0 * PI / 3.0);
        assert!((r1 - expected).norm() < 1e-10);
    }

    #[test]
    fn test_monodromy_identity() {
        let m = Monodromy::identity(3);
        assert!(m.is_identity());
        assert_eq!(m.apply(0), 0);
        assert_eq!(m.apply(1), 1);
    }

    #[test]
    fn test_monodromy_transposition() {
        let m = Monodromy::transposition(3, 0, 1);
        assert_eq!(m.apply(0), 1);
        assert_eq!(m.apply(1), 0);
        assert_eq!(m.apply(2), 2);
        assert!(!m.is_identity());
    }

    #[test]
    fn test_monodromy_compose() {
        let m1 = Monodromy::transposition(3, 0, 1);
        let m2 = Monodromy::transposition(3, 1, 2);
        let composed = m1.compose(&m2);
        assert_eq!(composed.apply(0), 2); // 0 -> 1 (by m1) -> 2 (by m2)
        assert_eq!(composed.apply(1), 0); // 1 -> 0 (by m1) -> 0 (by m2)
    }

    #[test]
    fn test_euler_characteristic() {
        // Sphere: n=1, b=0 -> χ=2
        assert_eq!(euler_characteristic(1, 0), 2);
        // Torus-like: n=2, b=4 -> χ=0
        assert_eq!(euler_characteristic(2, 4), 0);
    }

    #[test]
    fn test_genus() {
        assert_eq!(genus(2), 0); // sphere
        assert_eq!(genus(0), 1); // torus
    }

    #[test]
    fn test_branch_cut_on_cut() {
        let cut = BranchCut::new(c(-1.0, 0.0), c(1.0, 0.0));
        assert!(cut.is_on_cut(c(0.0, 0.0), 1e-6));
        assert!(!cut.is_on_cut(c(0.0, 1.0), 1e-6));
    }

    #[test]
    fn test_nth_root_surface() {
        let surface = RiemannSurface::nth_root_surface(4);
        assert_eq!(surface.n_sheets, 4);
        assert_eq!(surface.branch_points[0].ramification_index, 4);
    }
}
