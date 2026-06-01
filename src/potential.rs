//! Complex potential theory, logarithmic potential, equilibrium measure

use num_complex::Complex64;
use serde::{Deserialize, Serialize};

/// Logarithmic potential: U^μ(z) = ∫ log|z - t| dμ(t)
/// For discrete measure with atoms at points with weights.
pub fn logarithmic_potential(z: Complex64, atoms: &[(Complex64, f64)]) -> f64 {
    atoms
        .iter()
        .map(|(t, w)| w * (z - *t).norm().ln().max(-50.0))
        .sum()
}

/// Gradient of the logarithmic potential.
pub fn potential_gradient(z: Complex64, atoms: &[(Complex64, f64)]) -> Complex64 {
    atoms
        .iter()
        .map(|(t, w)| {
            let diff = z - *t;
            let r2 = diff.norm_sqr();
            if r2 < 1e-15 {
                Complex64::new(0.0, 0.0)
            } else {
                w / diff * diff.conj() / r2
            }
        })
        .fold(Complex64::new(0.0, 0.0), |acc, x| acc + x)
}

/// Discrete probability measure on the complex plane.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DiscreteMeasure {
    pub atoms: Vec<(Complex64, f64)>,
}

impl DiscreteMeasure {
    pub fn new(atoms: Vec<(Complex64, f64)>) -> Self {
        Self { atoms }
    }

    /// Create a uniform measure on n points.
    pub fn uniform(points: Vec<Complex64>) -> Self {
        let n = points.len() as f64;
        Self {
            atoms: points.into_iter().map(|p| (p, 1.0 / n)).collect(),
        }
    }

    /// Total mass.
    pub fn total_mass(&self) -> f64 {
        self.atoms.iter().map(|(_, w)| w).sum()
    }

    /// Energy of the measure: I(μ) = -∫∫ log|z - t| dμ(z) dμ(t).
    pub fn energy(&self) -> f64 {
        let n = self.atoms.len();
        let mut total = 0.0;
        for i in 0..n {
            for j in 0..n {
                if i != j {
                    let (zi, wi) = self.atoms[i];
                    let (zj, wj) = self.atoms[j];
                    total -= wi * wj * (zi - zj).norm().ln();
                }
            }
        }
        total
    }

    /// Potential at z.
    pub fn potential(&self, z: Complex64) -> f64 {
        logarithmic_potential(z, &self.atoms)
    }

    /// Robin constant for a compact set (approximation).
    pub fn robin_constant(&self) -> f64 {
        -self.energy()
    }

    /// Logarithmic capacity (transfinite diameter).
    pub fn capacity(&self) -> f64 {
        let energy = self.energy();
        (-energy).exp()
    }
}

/// Compute the equilibrium measure on the unit circle.
/// For the unit circle, it's the uniform distribution.
pub fn unit_circle_equilibrium(n_points: usize) -> DiscreteMeasure {
    let atoms: Vec<(Complex64, f64)> = (0..n_points)
        .map(|k| {
            let theta = 2.0 * std::f64::consts::PI * k as f64 / n_points as f64;
            (Complex64::from_polar(1.0, theta), 1.0 / n_points as f64)
        })
        .collect();
    DiscreteMeasure::new(atoms)
}

/// Chebyshev nodes on [-1, 1] mapped to complex plane.
pub fn chebyshev_nodes(n: usize) -> Vec<Complex64> {
    (0..n)
        .map(|k| {
            let theta = std::f64::consts::PI * (2 * k + 1) as f64 / (2 * n) as f64;
            Complex64::new(theta.cos(), 0.0)
        })
        .collect()
}

/// Fekete points (approximation via energy minimization).
pub fn approximate_fekete_points(n: usize, iterations: usize, learning_rate: f64) -> Vec<Complex64> {
    let mut points: Vec<Complex64> = (0..n)
        .map(|k| Complex64::from_polar(1.0, 2.0 * std::f64::consts::PI * k as f64 / n as f64))
        .collect();

    for _ in 0..iterations {
        for i in 0..n {
            let mut gradient = Complex64::new(0.0, 0.0);
            for j in 0..n {
                if i != j {
                    let diff = points[i] - points[j];
                    let r2 = diff.norm_sqr();
                    if r2 > 1e-15 {
                        // d/dz_i [-log|z_i - z_j|] = -1/(2*(z_i - z_j))
                        gradient -= 1.0 / (2.0 * diff);
                    }
                }
            }
            points[i] += learning_rate * gradient;
        }
    }

    points
}

/// Green's function for the complement of the unit disk.
/// g(z, ∞) = log|z|
pub fn exterior_green_function(z: Complex64) -> f64 {
    z.norm().ln().max(0.0)
}

/// Reduced modulus of a domain.
pub fn reduced_modulus(potential_at_boundary: f64, capacity: f64) -> f64 {
    potential_at_boundary + capacity.ln()
}

/// Superharmonic ordering: check if u ≥ v on a set.
pub fn superharmonic_check(
    u: &dyn Fn(Complex64) -> f64,
    v: &dyn Fn(Complex64) -> f64,
    test_points: &[Complex64],
) -> bool {
    test_points.iter().all(|&z| u(z) >= v(z) - 1e-10)
}

/// Dirichlet energy: ∫|∇u|² dA (approximated on a grid).
pub fn dirichlet_energy(
    u: &dyn Fn(f64, f64) -> f64,
    x_range: (f64, f64),
    y_range: (f64, f64),
    nx: usize,
    ny: usize,
    h: f64,
) -> f64 {
    let dx = (x_range.1 - x_range.0) / nx as f64;
    let dy = (y_range.1 - y_range.0) / ny as f64;
    let mut energy = 0.0;

    for i in 0..nx {
        for j in 0..ny {
            let x = x_range.0 + i as f64 * dx;
            let y = y_range.0 + j as f64 * dy;
            let ux = (u(x + h, y) - u(x - h, y)) / (2.0 * h);
            let uy = (u(x, y + h) - u(x, y - h)) / (2.0 * h);
            energy += (ux * ux + uy * uy) * dx * dy;
        }
    }
    energy
}

#[cfg(test)]
mod tests {
    use super::*;

    fn c(re: f64, im: f64) -> Complex64 {
        Complex64::new(re, im)
    }

    #[test]
    fn test_logarithmic_potential() {
        let atoms = vec![(c(0.0, 0.0), 1.0)];
        let u = logarithmic_potential(c(1.0, 0.0), &atoms);
        // log|1-0| = ln(1) = 0
        assert!(u.abs() < 1e-10);
    }

    #[test]
    fn test_logarithmic_potential_two_atoms() {
        let atoms = vec![(c(1.0, 0.0), 0.5), (c(-1.0, 0.0), 0.5)];
        let u = logarithmic_potential(c(0.0, 0.0), &atoms);
        // 0.5 * ln|0-1| + 0.5 * ln|0-(-1)| = 0.5 * 0 + 0.5 * 0 = 0
        assert!(u.abs() < 1e-10);
    }

    #[test]
    fn test_discrete_measure_total_mass() {
        let m = DiscreteMeasure::uniform(vec![c(0.0, 0.0), c(1.0, 0.0), c(2.0, 0.0)]);
        assert!((m.total_mass() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_discrete_measure_energy() {
        let m = DiscreteMeasure::uniform(vec![c(0.0, 0.0), c(1.0, 0.0)]);
        // I(μ) = -2 * (1/2) * (1/2) * ln(1) = 0
        let e = m.energy();
        assert!(e.abs() < 1e-10);
    }

    #[test]
    fn test_unit_circle_equilibrium() {
        let m = unit_circle_equilibrium(100);
        assert!((m.total_mass() - 1.0).abs() < 1e-10);
        // All points should be on the unit circle
        for (p, _) in &m.atoms {
            assert!((p.norm() - 1.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_chebyshev_nodes() {
        let nodes = chebyshev_nodes(5);
        assert_eq!(nodes.len(), 5);
        // All should be in [-1, 1]
        for n in &nodes {
            assert!(n.re >= -1.0 && n.re <= 1.0);
        }
    }

    #[test]
    fn test_fekete_points() {
        let points = approximate_fekete_points(4, 100, 0.01);
        assert_eq!(points.len(), 4);
        // Points should roughly stay on unit circle for capacity of disk
        for p in &points {
            assert!(p.norm() < 3.0);
        }
    }

    #[test]
    fn test_exterior_green_function() {
        let g = exterior_green_function(c(2.0, 0.0));
        assert!((g - 2.0_f64.ln()).abs() < 1e-10);
    }

    #[test]
    fn test_exterior_green_function_inside() {
        let g = exterior_green_function(c(0.5, 0.0));
        assert_eq!(g, 0.0);
    }

    #[test]
    fn test_dirichlet_energy() {
        // u(x,y) = x, energy should be area
        let u = |x: f64, _y: f64| x;
        let energy = dirichlet_energy(&u, (-1.0, 1.0), (-1.0, 1.0), 50, 50, 1e-4);
        // |∇u|² = 1, area = 4, so energy ≈ 4
        assert!((energy - 4.0).abs() < 0.5);
    }

    #[test]
    fn test_superharmonic_check() {
        let u = |z: Complex64| z.norm().ln();
        let v = |z: Complex64| z.norm().ln() - 1.0;
        let points: Vec<Complex64> = (0..10)
            .map(|k| Complex64::from_polar(2.0, std::f64::consts::PI * k as f64 / 10.0))
            .collect();
        assert!(superharmonic_check(&u, &v, &points));
    }

    #[test]
    fn test_capacity() {
        let m = DiscreteMeasure::uniform(vec![c(0.0, 0.0), c(1.0, 0.0)]);
        let cap = m.capacity();
        assert!(cap > 0.0);
    }
}
