//! Residue calculus, Laurent series, pole classification, argument principle

use num_complex::Complex64;

/// Classify the type of singularity of f at z0.
#[derive(Debug, Clone, PartialEq)]
pub enum SingularityType {
    Removable,
    Pole { order: usize },
    Essential,
    NotIsolated,
}

/// Compute the residue of f at z0 using numerical differentiation of (z-z0)^n * f(z).
/// Res(f, z0) = (1/(n-1)!) * d^(n-1)/dz^(n-1) [(z-z0)^n * f(z)] evaluated at z0
pub fn residue(f: &dyn Fn(Complex64) -> Complex64, z0: Complex64, order: usize, h: f64) -> Complex64 {
    if order == 1 {
        // Simple pole: Res = lim_{z→z0} (z-z0)*f(z)
        let g = |z: Complex64| (z - z0) * f(z);
        // Use numerical limit with Richardson extrapolation
        let r1 = g(z0 + h);
        let r2 = g(z0 + h / 2.0);
        let r3 = g(z0 + h / 4.0);
        // Richardson: extrapolate
        4.0 * r3 / 3.0 - r2 / 3.0
    } else {
        // Higher order pole
        let n = order;
        let factorial_n_minus_1 = (1..n).fold(1.0, |acc, k| acc * k as f64);
        // g(z) = (z - z0)^n * f(z)
        // We need the (n-1)th derivative of g at z0
        let g = |z: Complex64| (z - z0).powi(n as i32) * f(z);
        let deriv = numerical_nth_derivative(&g, z0, n - 1, h);
        deriv / factorial_n_minus_1
    }
}

/// Numerical nth derivative using finite differences.
fn numerical_nth_derivative(
    f: &dyn Fn(Complex64) -> Complex64,
    z0: Complex64,
    n: usize,
    h: f64,
) -> Complex64 {
    if n == 0 {
        return f(z0);
    }
    // Use Cauchy's integral formula for derivatives with a small circle
    let radius = h * 10.0;
    let m = 200;
    let factorial_n = (1..=n).fold(1.0, |acc, k| acc * k as f64);
    let pi = std::f64::consts::PI;
    let mut sum = Complex64::new(0.0, 0.0);
    for k in 0..m {
        let theta = 2.0 * pi * k as f64 / m as f64;
        let z = z0 + radius * Complex64::from_polar(1.0, theta);
        let dtheta = 2.0 * pi / m as f64;
        sum += f(z) / (z - z0).powi(n as i32 + 1) * Complex64::from_polar(radius, theta) * dtheta;
    }
    factorial_n / (2.0 * pi * Complex64::i()) * sum
}

/// Classify a singularity by analyzing the behavior near z0.
pub fn classify_singularity(
    f: &dyn Fn(Complex64) -> Complex64,
    z0: Complex64,
    h: f64,
) -> SingularityType {
    // Try approaching z0 from multiple directions
    let eps = h;
    let directions = [
        Complex64::new(1.0, 0.0),
        Complex64::new(0.0, 1.0),
        Complex64::new(-1.0, 0.0),
        Complex64::new(0.0, -1.0),
    ];

    let values: Vec<Complex64> = directions
        .iter()
        .map(|d| f(z0 + eps * d))
        .collect();

    // Check if values blow up
    let max_norm = values.iter().map(|v| v.norm()).fold(0.0_f64, f64::max);
    if max_norm > 1e3 {
        // Likely a pole. Try to determine order.
        for order in 1..=10 {
            let g = |z: Complex64| (z - z0).powi(order as i32) * f(z);
            let g_vals: Vec<Complex64> = directions.iter().map(|d| g(z0 + eps * d)).collect();
            let g_norms: Vec<f64> = g_vals.iter().map(|v| v.norm()).collect();
            let ratio = g_norms.iter().cloned().fold(f64::MAX, f64::min)
                / g_norms.iter().cloned().fold(0.0_f64, f64::max).max(1e-20);
            if ratio > 0.5 {
                return SingularityType::Pole { order };
            }
        }
        return SingularityType::Essential;
    }

    // Check if limit exists (removable)
    let mean: Complex64 = values.iter().sum::<Complex64>() / values.len() as f64;
    let variance: f64 = values.iter().map(|v| (v - mean).norm()).sum::<f64>() / values.len() as f64;
    if variance < 1e-4 {
        SingularityType::Removable
    } else {
        SingularityType::Essential
    }
}

/// Laurent series coefficient computation around z0.
/// Returns (negative_coeffs, non_negative_coeffs) where
/// f(z) = Σ a_k (z-z0)^k for k from -n to ∞
pub fn laurent_coefficients(
    f: &dyn Fn(Complex64) -> Complex64,
    z0: Complex64,
    inner_radius: f64,
    outer_radius: f64,
    n_negative: usize,
    n_positive: usize,
) -> (Vec<Complex64>, Vec<Complex64>) {
    let n_total = n_negative + n_positive;
    let m = 200; // integration points
    let pi = std::f64::consts::PI;

    let mut all_coeffs = Vec::with_capacity(n_total);
    for k in -(n_negative as i32)..=(n_positive as i32) {
        let radius = (inner_radius + outer_radius) / 2.0;
        let mut sum = Complex64::new(0.0, 0.0);
        for j in 0..m {
            let theta = 2.0 * pi * j as f64 / m as f64;
            let z = z0 + radius * Complex64::from_polar(1.0, theta);
            let dz = radius * Complex64::from_polar(1.0, theta) * (2.0 * pi / m as f64);
            sum += f(z) * (z - z0).powi(-k) * dz;
        }
        all_coeffs.push(sum / (2.0 * pi * Complex64::i()));
    }

    let negative: Vec<Complex64> = all_coeffs[..n_negative].to_vec();
    let positive: Vec<Complex64> = all_coeffs[n_negative..].to_vec();
    (negative, positive)
}

/// Argument principle: number of zeros - number of poles inside γ
/// N - P = (1/2πi) ∮ f'(z)/f(z) dz
pub fn argument_principle(
    f: &dyn Fn(Complex64) -> Complex64,
    f_prime: &dyn Fn(Complex64) -> Complex64,
    contour: &dyn crate::contour::Contour,
    n_points: usize,
) -> i64 {
    let integrand = |z: Complex64| f_prime(z) / f(z);
    let integral = crate::contour::contour_integral(&integrand, contour, n_points);
    let result = integral / (2.0 * std::f64::consts::PI * Complex64::i());
    result.re.round() as i64
}

/// Rouche's theorem: if |f(z) - g(z)| < |f(z)| on γ, then f and g have the same number of zeros inside γ.
pub fn rouche_check(
    f: &dyn Fn(Complex64) -> Complex64,
    g: &dyn Fn(Complex64) -> Complex64,
    contour: &dyn crate::contour::Contour,
    n_points: usize,
) -> bool {
    let dt = 1.0 / n_points as f64;
    for k in 0..n_points {
        let t = (k as f64 + 0.5) * dt;
        let z = contour.evaluate(t);
        let fz = f(z);
        let gz = g(z);
        if (fz - gz).norm() >= fz.norm() {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::contour::CircleContour;

    fn c(re: f64, im: f64) -> Complex64 {
        Complex64::new(re, im)
    }

    #[test]
    fn test_residue_simple_pole() {
        // f(z) = 1/z, Res at 0 should be 1
        let f = |z: Complex64| 1.0 / z;
        let res = residue(&f, c(0.0, 0.0), 1, 1e-6);
        assert!((res - c(1.0, 0.0)).norm() < 0.1);
    }

    #[test]
    fn test_residue_double_pole() {
        // f(z) = 1/z^2, Res at 0 should be 0
        let f = |z: Complex64| 1.0 / (z * z);
        let res = residue(&f, c(0.0, 0.0), 2, 1e-6);
        assert!(res.norm() < 0.1);
    }

    #[test]
    fn test_classify_pole() {
        let f = |z: Complex64| 1.0 / z;
        let st = classify_singularity(&f, c(0.0, 0.0), 1e-6);
        match st {
            SingularityType::Pole { .. } => {},
            _ => panic!("Expected pole, got {:?}", st),
        }
    }

    #[test]
    fn test_classify_removable() {
        // f(z) = sin(z)/z at z=0 is removable
        let f = |z: Complex64| z.sin() / z;
        let st = classify_singularity(&f, c(0.0, 0.0), 1e-4);
        assert_eq!(st, SingularityType::Removable);
    }

    #[test]
    fn test_argument_principle() {
        // f(z) = z^2 - 1 has 2 zeros inside |z|=2
        let f = |z: Complex64| z * z - 1.0;
        let fp = |z: Complex64| 2.0 * z;
        let circle = CircleContour::new(c(0.0, 0.0), 2.0);
        let n = argument_principle(&f, &fp, &circle, 10000);
        assert_eq!(n, 2);
    }

    #[test]
    fn test_rouche_check() {
        // f(z) = z^3, g(z) = z^3 + z on |z|=2
        // |g-f| = |z| = 2 < |z^3| = 8 on |z|=2, so same zeros
        let f = |z: Complex64| z * z * z;
        let g = |z: Complex64| z * z * z + z;
        let circle = CircleContour::new(c(0.0, 0.0), 2.0);
        assert!(rouche_check(&f, &g, &circle, 1000));
    }

    #[test]
    fn test_laurent_coefficients() {
        // f(z) = 1/z = z^(-1), should have a_{-1} = 1
        let f = |z: Complex64| 1.0 / z;
        let (neg, _pos) = laurent_coefficients(&f, c(0.0, 0.0), 0.5, 1.5, 3, 3);
        // neg[0]=a_{-3}, neg[1]=a_{-2}, neg[2]=a_{-1} ≈ 1
        assert!((neg[2] - c(1.0, 0.0)).norm() < 1.0);
    }

    #[test]
    fn test_residue_e_z_over_z() {
        // f(z) = e^z / z, Res at 0 = e^0 = 1
        let f = |z: Complex64| z.exp() / z;
        let res = residue(&f, c(0.0, 0.0), 1, 1e-6);
        assert!((res - c(1.0, 0.0)).norm() < 0.2);
    }
}
