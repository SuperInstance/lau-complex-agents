//! Holomorphic functions, Cauchy-Riemann equations, conformal maps, analytic continuation

use num_complex::Complex64;
use serde::{Deserialize, Serialize};

/// Check if a complex function satisfies the Cauchy-Riemann equations at a point.
/// Given u(x,y) and v(x,y), checks ∂u/∂x = ∂v/∂y and ∂u/∂y = -∂v/∂x
pub fn cauchy_riemann_check(
    u: &dyn Fn(f64, f64) -> f64,
    v: &dyn Fn(f64, f64) -> f64,
    x: f64,
    y: f64,
    h: f64,
) -> bool {
    let dudx = (u(x + h, y) - u(x - h, y)) / (2.0 * h);
    let dudy = (u(x, y + h) - u(x, y - h)) / (2.0 * h);
    let dvdx = (v(x + h, y) - v(x - h, y)) / (2.0 * h);
    let dvdy = (v(x, y + h) - v(x, y - h)) / (2.0 * h);
    (dudx - dvdy).abs() < 1e-6 && (dudy + dvdx).abs() < 1e-6
}

/// Compute the complex derivative f'(z) using central differences.
pub fn complex_derivative(f: &dyn Fn(Complex64) -> Complex64, z: Complex64, h: f64) -> Complex64 {
    (f(z + h) - f(z - h)) / (2.0 * h)
}

/// Check if a function is holomorphic at a point (Cauchy-Riemann + differentiable).
pub fn is_holomorphic(
    f: &dyn Fn(Complex64) -> Complex64,
    z: Complex64,
    h: f64,
) -> bool {
    let u = |x: f64, y: f64| f(Complex64::new(x, y)).re;
    let v = |x: f64, y: f64| f(Complex64::new(x, y)).im;
    cauchy_riemann_check(&u, &v, z.re, z.im, h)
}

/// Power series representation of a holomorphic function.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PowerSeries {
    /// Center of the series
    pub center: Complex64,
    /// Coefficients a_n such that f(z) = Σ a_n * (z - center)^n
    pub coefficients: Vec<Complex64>,
}

impl PowerSeries {
    pub fn new(center: Complex64, coefficients: Vec<Complex64>) -> Self {
        Self { center, coefficients }
    }

    /// Evaluate the power series at z.
    pub fn evaluate(&self, z: Complex64) -> Complex64 {
        let w = z - self.center;
        let mut result = Complex64::new(0.0, 0.0);
        let mut w_power = Complex64::new(1.0, 0.0);
        for coeff in &self.coefficients {
            result += coeff * w_power;
            w_power *= w;
        }
        result
    }

    /// Radius of convergence estimated via root test.
    pub fn radius_of_convergence(&self) -> Option<f64> {
        if self.coefficients.is_empty() {
            return None;
        }
        let n = self.coefficients.len() as f64;
        let max_coeff = self
            .coefficients
            .iter()
            .map(|c| c.norm())
            .fold(0.0_f64, f64::max);
        if max_coeff == 0.0 {
            return Some(f64::INFINITY);
        }
        Some(1.0 / max_coeff.powf(1.0 / n))
    }

    /// Differentiate the power series term-by-term.
    pub fn differentiate(&self) -> PowerSeries {
        if self.coefficients.is_empty() {
            return PowerSeries::new(self.center, vec![]);
        }
        let new_coeffs: Vec<Complex64> = self
            .coefficients
            .iter()
            .enumerate()
            .skip(1)
            .map(|(n, c)| c * n as f64)
            .collect();
        PowerSeries::new(self.center, new_coeffs)
    }

    /// Integrate the power series term-by-term (constant of integration = 0).
    pub fn integrate(&self) -> PowerSeries {
        let mut new_coeffs = vec![Complex64::new(0.0, 0.0)];
        for (n, c) in self.coefficients.iter().enumerate() {
            new_coeffs.push(c / (n as f64 + 1.0));
        }
        PowerSeries::new(self.center, new_coeffs)
    }
}

/// Analytic continuation by re-centering a power series.
pub fn analytic_continuation(series: &PowerSeries, new_center: Complex64, terms: usize) -> PowerSeries {
    // Taylor expand at new_center using the original series
    let mut coeffs = Vec::with_capacity(terms);
    let mut current_series = series.clone();
    let mut factorial = 1.0;
    for k in 0..terms {
        if k > 0 {
            factorial *= k as f64;
        }
        coeffs.push(current_series.evaluate(new_center) / factorial);
        current_series = current_series.differentiate();
    }
    PowerSeries::new(new_center, coeffs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex64;
    use std::f64::consts::PI;

    fn c(re: f64, im: f64) -> Complex64 {
        Complex64::new(re, im)
    }

    #[test]
    fn test_cauchy_riemann_exp() {
        let f = |z: Complex64| z.exp();
        assert!(is_holomorphic(&f, c(1.0, 1.0), 1e-7));
    }

    #[test]
    fn test_cauchy_riemann_conjugate() {
        // f(z) = conj(z) is NOT holomorphic
        let f = |z: Complex64| z.conj();
        assert!(!is_holomorphic(&f, c(1.0, 1.0), 1e-7));
    }

    #[test]
    fn test_complex_derivative_exp() {
        let f = |z: Complex64| z.exp();
        let z = c(1.0, 0.5);
        let deriv = complex_derivative(&f, z, 1e-8);
        let expected = z.exp(); // d/dz e^z = e^z
        assert!((deriv - expected).norm() < 1e-6);
    }

    #[test]
    fn test_complex_derivative_z2() {
        let f = |z: Complex64| z * z;
        let z = c(2.0, 3.0);
        let deriv = complex_derivative(&f, z, 1e-8);
        let expected = 2.0 * z; // d/dz z^2 = 2z
        assert!((deriv - expected).norm() < 1e-6);
    }

    #[test]
    fn test_power_series_exp() {
        // e^z = 1 + z + z^2/2! + z^3/3! + ...
        let coeffs: Vec<Complex64> = (0..20)
            .map(|n| c(1.0 / (1..=n).map(|k| k as f64).product::<f64>(), 0.0))
            .collect();
        let series = PowerSeries::new(c(0.0, 0.0), coeffs);
        let z = c(1.0, 1.0);
        let result = series.evaluate(z);
        let expected = z.exp();
        assert!((result - expected).norm() < 1e-10);
    }

    #[test]
    fn test_power_series_differentiate() {
        let coeffs = vec![c(1.0, 0.0), c(2.0, 0.0), c(3.0, 0.0)];
        // f(z) = 1 + 2z + 3z^2 => f'(z) = 2 + 6z
        let series = PowerSeries::new(c(0.0, 0.0), coeffs);
        let derived = series.differentiate();
        assert_eq!(derived.coefficients.len(), 2);
        assert!((derived.coefficients[0] - c(2.0, 0.0)).norm() < 1e-10);
        assert!((derived.coefficients[1] - c(6.0, 0.0)).norm() < 1e-10);
    }

    #[test]
    fn test_power_series_integrate() {
        // f(z) = 1 + 2z => F(z) = z + z^2
        let coeffs = vec![c(1.0, 0.0), c(2.0, 0.0)];
        let series = PowerSeries::new(c(0.0, 0.0), coeffs);
        let integrated = series.integrate();
        assert_eq!(integrated.coefficients.len(), 3);
        assert!((integrated.coefficients[0]).norm() < 1e-10); // 0
        assert!((integrated.coefficients[1] - c(1.0, 0.0)).norm() < 1e-10);
        assert!((integrated.coefficients[2] - c(1.0, 0.0)).norm() < 1e-10);
    }

    #[test]
    fn test_radius_of_convergence() {
        // Σ z^n has radius of convergence 1
        let coeffs: Vec<Complex64> = (0..50).map(|_| c(1.0, 0.0)).collect();
        let series = PowerSeries::new(c(0.0, 0.0), coeffs);
        let r = series.radius_of_convergence().unwrap();
        assert!((r - 1.0).abs() < 0.2); // approximate due to finite terms
    }

    #[test]
    fn test_analytic_continuation() {
        let coeffs: Vec<Complex64> = (0..20)
            .map(|n| c(1.0 / (1..=n).map(|k| k as f64).product::<f64>(), 0.0))
            .collect();
        let series = PowerSeries::new(c(0.0, 0.0), coeffs);
        let continued = analytic_continuation(&series, c(0.5, 0.0), 15);
        let z = c(0.6, 0.0);
        let result = continued.evaluate(z);
        let expected = z.exp();
        assert!((result - expected).norm() < 1e-6);
    }

    #[test]
    fn test_is_holomorphic_z_cubed() {
        let f = |z: Complex64| z * z * z;
        assert!(is_holomorphic(&f, c(1.0, -1.0), 1e-7));
    }

    #[test]
    fn test_cauchy_riemann_sin() {
        let f = |z: Complex64| z.sin();
        assert!(is_holomorphic(&f, c(0.5, 0.5), 1e-7));
    }
}
