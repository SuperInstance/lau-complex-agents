//! Julia sets, Mandelbrot set, Fatou components, iteration dynamics

use num_complex::Complex64;
use serde::{Deserialize, Serialize};

/// Result of iterating a complex function.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IterationResult {
    pub orbit: Vec<Complex64>,
    pub escaped: bool,
    pub escape_iteration: Option<usize>,
    pub converged: bool,
}

/// Iterate f(z) = z² + c starting from z0.
pub fn iterate_quadratic(z0: Complex64, c: Complex64, max_iterations: usize, escape_radius: f64) -> IterationResult {
    let mut orbit = Vec::with_capacity(max_iterations + 1);
    let mut z = z0;
    orbit.push(z);

    let mut escaped = false;
    let mut escape_iteration = None;
    let mut converged = false;

    for i in 0..max_iterations {
        z = z * z + c;
        orbit.push(z);

        if z.norm() > escape_radius {
            escaped = true;
            escape_iteration = Some(i + 1);
            break;
        }
        if orbit.len() > 2 && (z - orbit[orbit.len() - 2]).norm() < 1e-12 {
            converged = true;
            break;
        }
    }

    IterationResult {
        orbit,
        escaped,
        escape_iteration,
        converged,
    }
}

/// Check if c is in the Mandelbrot set (z0 = 0).
pub fn mandelbrot_membership(c: Complex64, max_iterations: usize, escape_radius: f64) -> (bool, Option<usize>) {
    let result = iterate_quadratic(Complex64::new(0.0, 0.0), c, max_iterations, escape_radius);
    (!result.escaped, result.escape_iteration)
}

/// Check if z0 is in the Julia set for parameter c.
pub fn julia_membership(z0: Complex64, c: Complex64, max_iterations: usize, escape_radius: f64) -> (bool, Option<usize>) {
    let result = iterate_quadratic(z0, c, max_iterations, escape_radius);
    (!result.escaped, result.escape_iteration)
}

/// Compute the filled Julia set on a grid.
/// Returns a 2D array where each cell is the iteration count (0 = in set).
pub fn filled_julia_set(
    c: Complex64,
    x_range: (f64, f64),
    y_range: (f64, f64),
    resolution: usize,
    max_iterations: usize,
    escape_radius: f64,
) -> Vec<Vec<usize>> {
    let mut grid = Vec::with_capacity(resolution);
    let dx = (x_range.1 - x_range.0) / resolution as f64;
    let dy = (y_range.1 - y_range.0) / resolution as f64;

    for j in 0..resolution {
        let mut row = Vec::with_capacity(resolution);
        for i in 0..resolution {
            let z0 = Complex64::new(
                x_range.0 + i as f64 * dx,
                y_range.0 + j as f64 * dy,
            );
            let result = iterate_quadratic(z0, c, max_iterations, escape_radius);
            row.push(result.escape_iteration.unwrap_or(0));
        }
        grid.push(row);
    }
    grid
}

/// Compute the Mandelbrot set on a grid.
pub fn mandelbrot_set(
    x_range: (f64, f64),
    y_range: (f64, f64),
    resolution: usize,
    max_iterations: usize,
    escape_radius: f64,
) -> Vec<Vec<usize>> {
    let mut grid = Vec::with_capacity(resolution);
    let dx = (x_range.1 - x_range.0) / resolution as f64;
    let dy = (y_range.1 - y_range.0) / resolution as f64;

    for j in 0..resolution {
        let mut row = Vec::with_capacity(resolution);
        for i in 0..resolution {
            let c = Complex64::new(
                x_range.0 + i as f64 * dx,
                y_range.0 + j as f64 * dy,
            );
            let result = iterate_quadratic(Complex64::new(0.0, 0.0), c, max_iterations, escape_radius);
            row.push(result.escape_iteration.unwrap_or(0));
        }
        grid.push(row);
    }
    grid
}

/// Classify a Fatou component.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum FatouComponent {
    /// Basin of attraction of an attracting fixed point.
    AttractingFixedPoint { fixed_point: Complex64 },
    /// Basin of attraction of a periodic cycle.
    PeriodicCycle { period: usize, points: Vec<Complex64> },
    /// Basin of attraction of infinity.
    BasinOfInfinity,
    /// Siegel disk.
    SiegelDisk,
    /// Unknown.
    Unknown,
}

/// Find fixed points of f(z) = z² + c: z = z² + c => z² - z + c = 0.
pub fn quadratic_fixed_points(c: Complex64) -> (Complex64, Complex64) {
    // z = (1 ± sqrt(1 - 4c)) / 2
    let disc = (Complex64::new(1.0, 0.0) - 4.0 * c).sqrt();
    let z1 = (Complex64::new(1.0, 0.0) + disc) / 2.0;
    let z2 = (Complex64::new(1.0, 0.0) - disc) / 2.0;
    (z1, z2)
}

/// Find periodic points of period n for f(z) = z² + c.
/// We iterate f^n(z) = z.
pub fn find_periodic_points(c: Complex64, period: usize, n_initial: usize) -> Vec<Complex64> {
    let mut periodic = Vec::new();

    // Sample points and check if they converge to periodic cycles
    for i in 0..n_initial {
        let angle = 2.0 * std::f64::consts::PI * i as f64 / n_initial as f64;
        let r = 0.5;
        let mut z = Complex64::from_polar(r, angle);

        // Iterate to find attracting cycles
        for _ in 0..1000 {
            for _ in 0..period {
                z = z * z + c;
            }
        }

        // Check if z is approximately periodic
        let z_after = {
            let mut zz = z;
            for _ in 0..period {
                zz = zz * zz + c;
            }
            zz
        };

        if (z_after - z).norm() < 1e-6 {
            // Check if this is a new cycle
            let is_new = !periodic.iter().any(|p: &Complex64| (p - z).norm() < 1e-3);
            if is_new {
                // Collect all points in the cycle
                let mut cycle = vec![z];
                let mut zz = z;
                for _ in 1..period {
                    zz = zz * zz + c;
                    cycle.push(zz);
                }
                periodic.extend(cycle);
            }
        }
    }
    periodic
}

/// Compute the multiplier of a fixed point: λ = f'(z*).
pub fn multiplier(f: &dyn Fn(Complex64) -> Complex64, f_prime: &dyn Fn(Complex64) -> Complex64, fixed_point: Complex64) -> Complex64 {
    f_prime(fixed_point)
}

/// Classify a fixed point by its multiplier.
#[derive(Debug, Clone, PartialEq)]
pub enum FixedPointType {
    SuperAttracting,
    Attracting,
    Repelling,
    Indifferent,
    Parabolic,
}

pub fn classify_fixed_point(lambda: Complex64) -> FixedPointType {
    let r = lambda.norm();
    if r < 0.99 {
        if r < 0.01 {
            FixedPointType::SuperAttracting
        } else {
            FixedPointType::Attracting
        }
    } else if r > 1.01 {
        FixedPointType::Repelling
    } else if (lambda.im).abs() < 0.01 {
        FixedPointType::Parabolic
    } else {
        FixedPointType::Indifferent
    }
}

/// Compute the Mandelbrot set boundary approximation via distance estimation.
pub fn mandelbrot_distance(c: Complex64, max_iterations: usize, escape_radius: f64) -> f64 {
    let mut z = Complex64::new(0.0, 0.0);
    let mut dz = Complex64::new(1.0, 0.0);

    for _ in 0..max_iterations {
        dz = 2.0 * z * dz + Complex64::new(1.0, 0.0);
        z = z * z + c;
        if z.norm() > escape_radius {
            let dist = z.norm() * (z.norm().ln() / dz.norm().ln()) / dz.norm();
            return dist;
        }
    }
    0.0 // Inside the set
}

#[cfg(test)]
mod tests {
    use super::*;

    fn c(re: f64, im: f64) -> Complex64 {
        Complex64::new(re, im)
    }

    #[test]
    fn test_mandelbrot_center() {
        let (in_set, _) = mandelbrot_membership(c(0.0, 0.0), 100, 10.0);
        assert!(in_set);
    }

    #[test]
    fn test_mandelbrot_outside() {
        let (in_set, esc) = mandelbrot_membership(c(2.0, 0.0), 100, 10.0);
        assert!(!in_set);
        assert!(esc.is_some());
    }

    #[test]
    fn test_mandelbrot_main_cardioid() {
        // c = -0.5 is inside the main cardioid
        let (in_set, _) = mandelbrot_membership(c(-0.5, 0.0), 200, 10.0);
        assert!(in_set);
    }

    #[test]
    fn test_julia_escaped() {
        let (in_set, _) = julia_membership(c(2.0, 0.0), c(0.0, 0.0), 100, 10.0);
        assert!(!in_set);
    }

    #[test]
    fn test_julia_bounded() {
        // For c = 0, Julia set is the unit circle
        let (in_set, _) = julia_membership(c(0.5, 0.0), c(0.0, 0.0), 100, 10.0);
        assert!(in_set);
    }

    #[test]
    fn test_iteration_orbit() {
        let result = iterate_quadratic(c(0.0, 0.0), c(1.0, 0.0), 10, 100.0);
        // 0 -> 1 -> 2 -> 5 -> 26 -> ...
        assert_eq!(result.orbit[0], c(0.0, 0.0));
        assert_eq!(result.orbit[1], c(1.0, 0.0));
        assert_eq!(result.orbit[2], c(2.0, 0.0));
    }

    #[test]
    fn test_quadratic_fixed_points() {
        let (z1, z2) = quadratic_fixed_points(c(0.0, 0.0));
        // z^2 = z => z = 0 or z = 1
        assert!((z1 - c(1.0, 0.0)).norm() < 1e-10 || (z2 - c(1.0, 0.0)).norm() < 1e-10);
        assert!((z1 - c(0.0, 0.0)).norm() < 1e-10 || (z2 - c(0.0, 0.0)).norm() < 1e-10);
    }

    #[test]
    fn test_classify_superattracting() {
        assert_eq!(classify_fixed_point(c(0.0, 0.0)), FixedPointType::SuperAttracting);
    }

    #[test]
    fn test_classify_attracting() {
        assert_eq!(classify_fixed_point(c(0.5, 0.0)), FixedPointType::Attracting);
    }

    #[test]
    fn test_classify_repelling() {
        assert_eq!(classify_fixed_point(c(2.0, 0.0)), FixedPointType::Repelling);
    }

    #[test]
    fn test_filled_julia_set_size() {
        let grid = filled_julia_set(
            c(-1.0, 0.0),
            (-2.0, 2.0),
            (-2.0, 2.0),
            50,
            100,
            10.0,
        );
        assert_eq!(grid.len(), 50);
        assert_eq!(grid[0].len(), 50);
    }

    #[test]
    fn test_mandelbrot_set_size() {
        let grid = mandelbrot_set(
            (-2.5, 1.0),
            (-1.5, 1.5),
            50,
            100,
            10.0,
        );
        assert_eq!(grid.len(), 50);
    }

    #[test]
    fn test_mandelbrot_distance_inside() {
        let d = mandelbrot_distance(c(0.0, 0.0), 100, 10.0);
        assert_eq!(d, 0.0);
    }

    #[test]
    fn test_mandelbrot_distance_outside() {
        let d = mandelbrot_distance(c(3.0, 0.0), 100, 10.0);
        assert!(d > 0.0);
    }
}
