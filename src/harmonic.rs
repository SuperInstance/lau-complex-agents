//! Harmonic functions, Dirichlet problem, Poisson kernel, mean value property

use nalgebra::DMatrix;
use num_complex::Complex64;

/// Check the mean value property for a function u at point (x0, y0).
/// u(x0, y0) should equal (1/2π) ∫ u(x0 + r cos θ, y0 + r sin θ) dθ
pub fn mean_value_property(
    u: &dyn Fn(f64, f64) -> f64,
    x0: f64,
    y0: f64,
    r: f64,
    n_points: usize,
) -> (f64, f64) {
    let u_center = u(x0, y0);
    let pi = std::f64::consts::PI;
    let mut integral = 0.0;
    for k in 0..n_points {
        let theta = 2.0 * pi * k as f64 / n_points as f64;
        integral += u(x0 + r * theta.cos(), y0 + r * theta.sin());
    }
    let mean = integral / n_points as f64;
    (u_center, mean)
}

/// Check if a function satisfies the mean value property (i.e., is harmonic).
pub fn is_harmonic(
    u: &dyn Fn(f64, f64) -> f64,
    x0: f64,
    y0: f64,
    r: f64,
    tolerance: f64,
) -> bool {
    let (center, mean) = mean_value_property(u, x0, y0, r, 200);
    (center - mean).abs() < tolerance
}

/// Compute the Laplacian ∇²u numerically.
pub fn laplacian(u: &dyn Fn(f64, f64) -> f64, x: f64, y: f64, h: f64) -> f64 {
    let d2udx2 = (u(x + h, y) - 2.0 * u(x, y) + u(x - h, y)) / (h * h);
    let d2udy2 = (u(x, y + h) - 2.0 * u(x, y) + u(x, y - h)) / (h * h);
    d2udx2 + d2udy2
}

/// Poisson kernel for the unit disk.
/// P(r, θ) = (1 - r²) / (1 - 2r cos θ + r²)
pub fn poisson_kernel(r: f64, theta: f64) -> f64 {
    let denom = 1.0 - 2.0 * r * theta.cos() + r * r;
    (1.0 - r * r) / denom
}

/// Solve the Dirichlet problem on the unit disk using Poisson integral.
/// Given boundary values on the unit circle, compute u at interior point.
pub fn dirichlet_disk(
    boundary_values: &dyn Fn(f64) -> f64, // f(θ) gives boundary value at angle θ
    r: f64,
    theta: f64,
    n_points: usize,
) -> f64 {
    let pi = std::f64::consts::PI;
    let mut sum = 0.0;
    let dphi = 2.0 * pi / n_points as f64;
    for k in 0..n_points {
        let phi = k as f64 * dphi;
        let p = poisson_kernel(r, theta - phi);
        sum += p * boundary_values(phi) * dphi;
    }
    sum / (2.0 * pi)
}

/// Solve the Dirichlet problem on a rectangle using finite differences.
pub fn dirichlet_rectangle(
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
    boundary: &dyn Fn(usize, usize) -> f64,
    max_iterations: usize,
    tolerance: f64,
) -> DMatrix<f64> {
    let mut grid = DMatrix::zeros(nx, ny);

    // Set boundary conditions
    for i in 0..nx {
        grid[(i, 0)] = boundary(i, 0);
        grid[(i, ny - 1)] = boundary(i, ny - 1);
    }
    for j in 0..ny {
        grid[(0, j)] = boundary(0, j);
        grid[(nx - 1, j)] = boundary(nx - 1, j);
    }

    let rx = 1.0 / (dx * dx);
    let ry = 1.0 / (dy * dy);
    let factor = 1.0 / (2.0 * (rx + ry));

    // Gauss-Seidel iteration
    for _ in 0..max_iterations {
        let mut max_change = 0.0;
        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                let old = grid[(i, j)];
                let new_val = factor * (
                    rx * (grid[(i - 1, j)] + grid[(i + 1, j)]) +
                    ry * (grid[(i, j - 1)] + grid[(i, j + 1)])
                );
                grid[(i, j)] = new_val;
                let change = (new_val - old).abs();
                if change > max_change {
                    max_change = change;
                }
            }
        }
        if max_change < tolerance {
            break;
        }
    }

    grid
}

/// Find the harmonic conjugate v of u (such that u + iv is analytic).
/// Uses path integration along horizontal then vertical.
pub fn harmonic_conjugate(
    u: &dyn Fn(f64, f64) -> f64,
    x0: f64,
    y0: f64,
    x: f64,
    y: f64,
    h: f64,
) -> f64 {
    // v_x = -u_y and v_y = u_x (Cauchy-Riemann)
    let mut v = 0.0;
    // Integrate along x first (y = y0 fixed)
    let n_steps_x = ((x - x0) / h).abs().ceil() as usize;
    let dx = (x - x0) / n_steps_x.max(1) as f64;
    let mut cx = x0;
    for _ in 0..n_steps_x {
        // v_x = -u_y
        let duy = (u(cx, y0 + h) - u(cx, y0 - h)) / (2.0 * h);
        v += -duy * dx;
        cx += dx;
    }
    // Then along y (x fixed)
    let n_steps_y = ((y - y0) / h).abs().ceil() as usize;
    let dy_step = (y - y0) / n_steps_y.max(1) as f64;
    let mut cy = y0;
    for _ in 0..n_steps_y {
        // v_y = u_x
        let dux = (u(x, cy + h) - u(x, cy - h)) / (2.0 * h);
        v += dux * dy_step;
        cy += dy_step;
    }
    v
}

/// Green's function for the unit disk.
/// G(z, ζ) = (1/(2π)) * log|1 - z̄ζ| - (1/(2π)) * log|z - ζ|
pub fn greens_function_disk(z: Complex64, zeta: Complex64) -> f64 {
    let term1 = (1.0 - z.conj() * zeta).norm().ln();
    let term2 = (z - zeta).norm().ln();
    (term1 - term2) / (2.0 * std::f64::consts::PI)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn c(re: f64, im: f64) -> Complex64 {
        Complex64::new(re, im)
    }

    #[test]
    fn test_mean_value_harmonic() {
        // u(x,y) = x^2 - y^2 is harmonic
        let u = |x: f64, y: f64| x * x - y * y;
        let (center, mean) = mean_value_property(&u, 1.0, 1.0, 0.1, 1000);
        assert!((center - mean).abs() < 1e-4);
    }

    #[test]
    fn test_is_harmonic() {
        let u = |x: f64, y: f64| x * x - y * y;
        assert!(is_harmonic(&u, 0.0, 0.0, 0.1, 1e-4));
    }

    #[test]
    fn test_laplacian_harmonic() {
        let u = |x: f64, y: f64| x * x - y * y;
        let lap = laplacian(&u, 1.0, 1.0, 1e-5);
        assert!(lap.abs() < 1e-4);
    }

    #[test]
    fn test_laplacian_non_harmonic() {
        let u = |x: f64, y: f64| x * x * x; // ∇²u = 6x ≠ 0
        let lap = laplacian(&u, 1.0, 0.0, 1e-5);
        assert!((lap - 6.0).abs() < 0.1);
    }

    #[test]
    fn test_poisson_kernel_at_center() {
        // P(0, θ) = 1 for all θ
        assert!((poisson_kernel(0.0, 0.0) - 1.0).abs() < 1e-10);
        assert!((poisson_kernel(0.0, PI) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_poisson_kernel_positive() {
        for r in &[0.0, 0.3, 0.5, 0.9] {
            for theta in &[0.0, PI / 4.0, PI / 2.0, PI] {
                assert!(poisson_kernel(*r, *theta) > 0.0);
            }
        }
    }

    #[test]
    fn test_dirichlet_disk_constant() {
        let boundary = |_theta: f64| 5.0;
        let result = dirichlet_disk(&boundary, 0.5, PI / 4.0, 1000);
        assert!((result - 5.0).abs() < 0.1);
    }

    #[test]
    fn test_dirichlet_disk_cos() {
        // If f(θ) = cos θ, then u(r, θ) = r cos θ
        let boundary = |theta: f64| theta.cos();
        let r = 0.5;
        let theta = PI / 4.0;
        let result = dirichlet_disk(&boundary, r, theta, 1000);
        let expected = r * theta.cos();
        assert!((result - expected).abs() < 0.05);
    }

    #[test]
    fn test_dirichlet_rectangle() {
        let boundary = |_i: usize, _j: usize| 0.0;
        let grid = dirichlet_rectangle(10, 10, 0.1, 0.1, &boundary, 1000, 1e-6);
        // All zeros
        assert!(grid.iter().all(|&v| v.abs() < 1e-6));
    }

    #[test]
    fn test_greens_function_disk() {
        // G(z, ζ) = 0 when z is on the boundary |z| = 1
        let z = Complex64::from_polar(1.0, 0.0);
        let zeta = c(0.5, 0.0);
        let g = greens_function_disk(z, zeta);
        assert!(g.abs() < 0.1);
    }

    #[test]
    fn test_greens_function_disk_interior() {
        let z = c(0.0, 0.0);
        let zeta = c(0.5, 0.0);
        let g = greens_function_disk(z, zeta);
        // G should be 0 on boundary, positive for exterior, negative for interior
        // For z inside disk (z≠ζ): G = (ln|1 - z̄ζ| - ln|z-ζ|) / 2π
        // At z=0: G = (ln|1| - ln|0.5|) / 2π = (0 + ln2) / 2π > 0
        assert!(g > 0.0);
    }
}
