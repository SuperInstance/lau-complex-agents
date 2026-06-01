//! Entire functions, Liouville theorem, order and type, Weierstrass factorization

use num_complex::Complex64;
use serde::{Deserialize, Serialize};

/// Compute the order of an entire function f.
/// ρ = lim sup_{r→∞} (log log M(r)) / (log r)
/// where M(r) = max_{|z|=r} |f(z)|
pub fn compute_order(
    f: &dyn Fn(Complex64) -> Complex64,
    radii: &[f64],
    n_theta: usize,
) -> f64 {
    let pi = std::f64::consts::PI;
    let mut max_order = 0.0_f64;

    for r in radii {
        let mut max_mod = 0.0_f64;
        for k in 0..n_theta {
            let theta = 2.0 * pi * k as f64 / n_theta as f64;
            let z = Complex64::from_polar(*r, theta);
            let val = f(z).norm();
            if val > max_mod {
                max_mod = val;
            }
        }
        if max_mod > 1.0 && *r > 1.0 {
            let log_log_m = max_mod.ln().ln();
            let log_r = r.ln();
            let order = log_log_m / log_r;
            if order > max_order {
                max_order = order;
            }
        }
    }
    max_order
}

/// Check Liouville's theorem: if f is entire and bounded, then f is constant.
pub fn check_liouville(
    f: &dyn Fn(Complex64) -> Complex64,
    radii: &[f64],
    n_theta: usize,
) -> (bool, f64) {
    let pi = std::f64::consts::PI;
    let mut all_values = Vec::new();

    for r in radii {
        for k in 0..n_theta {
            let theta = 2.0 * pi * k as f64 / n_theta as f64;
            let z = Complex64::from_polar(*r, theta);
            all_values.push(f(z));
        }
    }

    let max_mod = all_values.iter().map(|v| v.norm()).fold(0.0_f64, f64::max);
    let min_mod = all_values.iter().map(|v| v.norm()).fold(f64::MAX, f64::min);

    // If bounded and variation is small, it's essentially constant
    let variation = max_mod - min_mod;
    let is_bounded = max_mod.is_finite() && max_mod < 1e10;
    (is_bounded && variation < 0.01 * max_mod.max(1.0), max_mod)
}

/// Compute the type of an entire function of order ρ.
/// σ = lim sup_{r→∞} (log M(r)) / r^ρ
pub fn compute_type(
    f: &dyn Fn(Complex64) -> Complex64,
    order: f64,
    radii: &[f64],
    n_theta: usize,
) -> f64 {
    let pi = std::f64::consts::PI;
    let mut max_type = 0.0_f64;

    for r in radii {
        let mut max_mod = 0.0_f64;
        for k in 0..n_theta {
            let theta = 2.0 * pi * k as f64 / n_theta as f64;
            let z = Complex64::from_polar(*r, theta);
            let val = f(z).norm();
            if val > max_mod {
                max_mod = val;
            }
        }
        if max_mod > 1.0 && *r > 1.0 && order > 0.0 {
            let sigma = max_mod.ln() / r.powf(order);
            if sigma > max_type {
                max_type = sigma;
            }
        }
    }
    max_type
}

/// Canonical Weierstrass factor.
/// E_p(z) = (1 - z) exp(z + z²/2 + ... + z^p/p)
pub fn weierstrass_factor(z: Complex64, p: usize) -> Complex64 {
    let one_minus_z = Complex64::new(1.0, 0.0) - z;
    let mut exponent = Complex64::new(0.0, 0.0);
    for k in 1..=p {
        exponent += z.powi(k as i32) / k as f64;
    }
    one_minus_z * exponent.exp()
}

/// Weierstrass product for a sequence of zeros.
/// f(z) = z^m * e^{g(z)} * Π E_p(z/a_n)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WeierstrassProduct {
    /// Multiplicity of zero at origin.
    pub origin_order: usize,
    /// Zeros (excluding origin).
    pub zeros: Vec<Complex64>,
    /// Genus of the canonical product.
    pub genus: usize,
    /// Entire function g(z) as coefficients of power series.
    pub g_coefficients: Vec<Complex64>,
}

impl WeierstrassProduct {
    pub fn new(origin_order: usize, zeros: Vec<Complex64>, genus: usize) -> Self {
        Self {
            origin_order,
            zeros,
            genus,
            g_coefficients: vec![Complex64::new(0.0, 0.0)],
        }
    }

    /// Evaluate the Weierstrass product at z.
    pub fn evaluate(&self, z: Complex64) -> Complex64 {
        // z^m
        let mut result = z.powi(self.origin_order as i32);

        // Product of Weierstrass factors
        for a_n in &self.zeros {
            if a_n.norm() > 1e-15 {
                result *= weierstrass_factor(z / a_n, self.genus);
            }
        }

        // exp(g(z))
        let mut g_z = Complex64::new(0.0, 0.0);
        let mut z_power = Complex64::new(1.0, 0.0);
        for coeff in &self.g_coefficients {
            g_z += coeff * z_power;
            z_power *= z;
        }
        result * g_z.exp()
    }

    /// Compute the genus from the convergence exponent.
    /// p = smallest integer such that Σ 1/|a_n|^(p+1) converges.
    pub fn compute_genus(zeros: &[Complex64]) -> usize {
        let mut p = 0;
        loop {
            let sum: f64 = zeros.iter().map(|a| 1.0 / a.norm().powi((p + 1) as i32)).sum();
            if sum.is_finite() || p > 20 {
                return p;
            }
            p += 1;
        }
    }
}

/// Hadamard factorization: given the zeros, construct the Hadamard product.
pub fn hadamard_factorization(
    zeros: Vec<Complex64>,
    origin_order: usize,
    order: f64,
) -> WeierstrassProduct {
    let genus = (order.floor() as usize).max(0);
    WeierstrassProduct::new(origin_order, zeros, genus)
}

/// Fundamental theorem of algebra (Liouville-based proof):
/// Every non-constant polynomial has a root.
pub fn fundamental_theorem_of_algebra(coefficients: &[f64]) -> Vec<Complex64> {
    let n = coefficients.len();
    if n <= 1 {
        return vec![];
    }

    // Use Durand-Kerner method to find all roots
    let degree = n - 1;
    let mut roots: Vec<Complex64> = (0..degree)
        .map(|k| Complex64::from_polar(1.0, 2.0 * std::f64::consts::PI * k as f64 / degree as f64 + 0.1))
        .collect();

    for _ in 0..100 {
        let mut new_roots = roots.clone();
        for i in 0..degree {
            // Evaluate polynomial at roots[i]
            let mut p = Complex64::new(coefficients[0], 0.0);
            let mut z_power = roots[i];
            for &c in &coefficients[1..] {
                p = p * roots[i] + Complex64::new(c, 0.0);
                z_power *= roots[i];
            }

            // Actually recompute p properly
            p = Complex64::new(coefficients[0], 0.0);
            for j in 1..n {
                p = p * roots[i] + Complex64::new(coefficients[j], 0.0);
            }

            let mut denom = Complex64::new(1.0, 0.0);
            for j in 0..degree {
                if j != i {
                    denom *= roots[i] - roots[j];
                }
            }
            if denom.norm() > 1e-15 {
                new_roots[i] = roots[i] - p / denom;
            }
        }
        roots = new_roots;
    }

    roots
}

#[cfg(test)]
mod tests {
    use super::*;

    fn c(re: f64, im: f64) -> Complex64 {
        Complex64::new(re, im)
    }

    #[test]
    fn test_order_exp() {
        // e^z has order 1
        let f = |z: Complex64| z.exp();
        let radii: Vec<f64> = (1..10).map(|r| r as f64).collect();
        let order = compute_order(&f, &radii, 100);
        assert!((order - 1.0).abs() < 0.5);
    }

    #[test]
    fn test_order_exp_z_squared() {
        // e^{z^2} has order 2
        let f = |z: Complex64| (z * z).exp();
        let radii: Vec<f64> = (1..5).map(|r| r as f64).collect();
        let order = compute_order(&f, &radii, 100);
        assert!((order - 2.0).abs() < 1.0);
    }

    #[test]
    fn test_liouville_constant() {
        let f = |_z: Complex64| c(5.0, 3.0);
        let (is_const, _) = check_liouville(&f, &[1.0, 2.0, 5.0], 50);
        assert!(is_const);
    }

    #[test]
    fn test_liouville_not_constant() {
        let f = |z: Complex64| z;
        let (is_const, _) = check_liouville(&f, &[1.0, 5.0, 10.0], 50);
        assert!(!is_const);
    }

    #[test]
    fn test_weierstrass_factor_zero() {
        // E_p(0) = 1
        let result = weierstrass_factor(c(0.0, 0.0), 2);
        assert!((result - c(1.0, 0.0)).norm() < 1e-10);
    }

    #[test]
    fn test_weierstrass_product_simple() {
        // Product for zeros at 1 and -1
        let wp = WeierstrassProduct::new(0, vec![c(1.0, 0.0), c(-1.0, 0.0)], 1);
        // At z=0, should be E_1(0)*E_1(0) = 1
        assert!((wp.evaluate(c(0.0, 0.0)) - c(1.0, 0.0)).norm() < 1e-10);
        // At z=1, should be 0
        assert!(wp.evaluate(c(1.0, 0.0)).norm() < 1e-6);
    }

    #[test]
    fn test_compute_genus() {
        // Finitely many zeros: genus 0
        let zeros = vec![c(1.0, 0.0), c(2.0, 0.0), c(3.0, 0.0)];
        assert_eq!(WeierstrassProduct::compute_genus(&zeros), 0);
    }

    #[test]
    fn test_hadamard_factorization() {
        let zeros = vec![c(1.0, 0.0), c(-1.0, 0.0), c(0.0, 1.0)];
        let wp = hadamard_factorization(zeros.clone(), 0, 1.0);
        for z in &zeros {
            assert!(wp.evaluate(*z).norm() < 1e-4);
        }
    }

    #[test]
    fn test_fundamental_theorem_cubic() {
        // z^3 - 1 = 0 has roots 1, e^{2πi/3}, e^{4πi/3}
        let coeffs = [1.0, 0.0, 0.0, -1.0];
        let roots = fundamental_theorem_of_algebra(&coeffs);
        assert_eq!(roots.len(), 3);
        // Check each root
        for r in &roots {
            let val: Complex64 = r * r * r;
            assert!((val - c(1.0, 0.0)).norm() < 0.1);
        }
    }

    #[test]
    fn test_compute_type_exp() {
        let f = |z: Complex64| z.exp();
        let radii: Vec<f64> = (1..8).map(|r| r as f64).collect();
        let sigma = compute_type(&f, 1.0, &radii, 100);
        assert!(sigma > 0.0);
    }
}
