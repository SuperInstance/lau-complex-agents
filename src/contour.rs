//! Contour integration, Cauchy integral formula, Cauchy theorem

use num_complex::Complex64;

/// A parameterized contour γ: [a,b] → ℂ
pub trait Contour {
    /// Evaluate the contour at parameter t ∈ [0, 1].
    fn evaluate(&self, t: f64) -> Complex64;
    /// Derivative γ'(t).
    fn derivative(&self, t: f64) -> Complex64;
}

/// Circle contour centered at `center` with radius `r`, traversed counterclockwise.
pub struct CircleContour {
    pub center: Complex64,
    pub radius: f64,
}

impl CircleContour {
    pub fn new(center: Complex64, radius: f64) -> Self {
        Self { center, radius }
    }
}

impl Contour for CircleContour {
    fn evaluate(&self, t: f64) -> Complex64 {
        let theta = 2.0 * std::f64::consts::PI * t;
        self.center + self.radius * Complex64::from_polar(1.0, theta)
    }
    fn derivative(&self, t: f64) -> Complex64 {
        let theta = 2.0 * std::f64::consts::PI * t;
        self.radius * 2.0 * std::f64::consts::PI * Complex64::from_polar(1.0, theta + std::f64::consts::FRAC_PI_2)
    }
}

/// Line segment from `start` to `end`.
pub struct LineSegment {
    pub start: Complex64,
    pub end: Complex64,
}

impl LineSegment {
    pub fn new(start: Complex64, end: Complex64) -> Self {
        Self { start, end }
    }
}

impl Contour for LineSegment {
    fn evaluate(&self, t: f64) -> Complex64 {
        self.start + t * (self.end - self.start)
    }
    fn derivative(&self, _t: f64) -> Complex64 {
        self.end - self.start
    }
}

/// Polygon contour defined by vertices.
pub struct PolygonContour {
    pub vertices: Vec<Complex64>,
}

impl PolygonContour {
    pub fn new(vertices: Vec<Complex64>) -> Self {
        Self { vertices }
    }
}

impl Contour for PolygonContour {
    fn evaluate(&self, t: f64) -> Complex64 {
        let n = self.vertices.len();
        let scaled = t * n as f64;
        let seg = (scaled.floor() as usize).min(n - 1);
        let local_t = scaled - seg as f64;
        let start = self.vertices[seg];
        let end = self.vertices[(seg + 1) % n];
        start + local_t * (end - start)
    }
    fn derivative(&self, _t: f64) -> Complex64 {
        // Simplified: not segment-aware
        Complex64::new(0.0, 0.0)
    }
}

/// Numerically integrate f along a contour using the trapezoidal rule.
/// ∫_γ f(z) dz = Σ f(γ(t_k)) γ'(t_k) Δt
pub fn contour_integral(f: &dyn Fn(Complex64) -> Complex64, contour: &dyn Contour, n: usize) -> Complex64 {
    let dt = 1.0 / n as f64;
    let mut sum = Complex64::new(0.0, 0.0);
    for k in 0..n {
        let t = (k as f64 + 0.5) * dt;
        sum += f(contour.evaluate(t)) * contour.derivative(t) * dt;
    }
    sum
}

/// Cauchy integral formula: f(a) = (1/2πi) ∮_γ f(z)/(z-a) dz
pub fn cauchy_integral_formula(
    f: &dyn Fn(Complex64) -> Complex64,
    a: Complex64,
    contour: &dyn Contour,
    n: usize,
) -> Complex64 {
    let integrand = |z: Complex64| f(z) / (z - a);
    let integral = contour_integral(&integrand, contour, n);
    integral / (2.0 * std::f64::consts::PI * Complex64::i())
}

/// Cauchy's integral formula for the nth derivative.
/// f^(n)(a) = n!/(2πi) ∮ f(z)/(z-a)^(n+1) dz
pub fn cauchy_nth_derivative(
    f: &dyn Fn(Complex64) -> Complex64,
    a: Complex64,
    n_order: usize,
    contour: &dyn Contour,
    n_points: usize,
) -> Complex64 {
    let n_f64 = n_order as f64;
    let integrand = move |z: Complex64| {
        let denom = (z - a).powi(n_order as i32 + 1);
        f(z) / denom
    };
    let integral = contour_integral(&integrand, contour, n_points);
    let factorial = (1..=n_order).fold(1.0, |acc, k| acc * k as f64);
    factorial / (2.0 * std::f64::consts::PI) * integral / Complex64::i()
}

/// Verify Cauchy's theorem: ∮_γ f(z) dz = 0 for holomorphic f in simply connected domain.
pub fn cauchy_theorem_check(
    f: &dyn Fn(Complex64) -> Complex64,
    contour: &dyn Contour,
    n: usize,
    tolerance: f64,
) -> bool {
    let integral = contour_integral(f, contour, n);
    integral.norm() < tolerance
}

/// Winding number of contour γ around point a.
/// n(γ, a) = (1/2πi) ∮_γ 1/(z-a) dz
pub fn winding_number(contour: &dyn Contour, a: Complex64, n: usize) -> i64 {
    let f = |z: Complex64| 1.0 / (z - a);
    let integral = contour_integral(&f, contour, n);
    let wn = integral / (2.0 * std::f64::consts::PI * Complex64::i());
    (wn.re.round()) as i64
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn c(re: f64, im: f64) -> Complex64 {
        Complex64::new(re, im)
    }

    #[test]
    fn test_circle_contour_at_zero() {
        let circle = CircleContour::new(c(0.0, 0.0), 1.0);
        let p0 = circle.evaluate(0.0);
        let p25 = circle.evaluate(0.25);
        assert!((p0.re - 1.0).abs() < 1e-10);
        assert!((p25.im - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_contour_integral_holomorphic() {
        // ∫_circle e^z dz = 0 (Cauchy's theorem)
        let circle = CircleContour::new(c(0.0, 0.0), 1.0);
        let f = |z: Complex64| z.exp();
        let integral = contour_integral(&f, &circle, 10000);
        assert!(integral.norm() < 1e-6);
    }

    #[test]
    fn test_contour_integral_1_over_z() {
        // ∮ 1/z dz = 2πi around origin
        let circle = CircleContour::new(c(0.0, 0.0), 2.0);
        let f = |z: Complex64| 1.0 / z;
        let integral = contour_integral(&f, &circle, 10000);
        let expected = 2.0 * PI * Complex64::i();
        assert!((integral - expected).norm() < 1e-4);
    }

    #[test]
    fn test_cauchy_integral_formula() {
        // f(z) = e^z, evaluate f(0.5) via CIF around circle of radius 2
        let f = |z: Complex64| z.exp();
        let circle = CircleContour::new(c(0.0, 0.0), 2.0);
        let result = cauchy_integral_formula(&f, c(0.5, 0.0), &circle, 10000);
        let expected = c(0.5, 0.0).exp();
        assert!((result - expected).norm() < 1e-4);
    }

    #[test]
    fn test_cauchy_theorem() {
        let circle = CircleContour::new(c(0.0, 0.0), 1.0);
        let f = |z: Complex64| z * z;
        assert!(cauchy_theorem_check(&f, &circle, 10000, 1e-6));
    }

    #[test]
    fn test_winding_number_inside() {
        let circle = CircleContour::new(c(0.0, 0.0), 1.0);
        let wn = winding_number(&circle, c(0.0, 0.0), 10000);
        assert_eq!(wn, 1);
    }

    #[test]
    fn test_winding_number_outside() {
        let circle = CircleContour::new(c(0.0, 0.0), 1.0);
        let wn = winding_number(&circle, c(5.0, 5.0), 10000);
        assert_eq!(wn, 0);
    }

    #[test]
    fn test_cauchy_nth_derivative() {
        // f(z) = e^z, f'(0) = 1
        let f = |z: Complex64| z.exp();
        let circle = CircleContour::new(c(0.0, 0.0), 1.0);
        let result = cauchy_nth_derivative(&f, c(0.0, 0.0), 1, &circle, 10000);
        assert!((result - c(1.0, 0.0)).norm() < 1e-3);
    }

    #[test]
    fn test_line_segment() {
        let seg = LineSegment::new(c(0.0, 0.0), c(1.0, 1.0));
        let mid = seg.evaluate(0.5);
        assert!((mid - c(0.5, 0.5)).norm() < 1e-10);
    }

    #[test]
    fn test_polygon_contour() {
        let square = PolygonContour::new(vec![c(0.0, 0.0), c(1.0, 0.0), c(1.0, 1.0), c(0.0, 1.0)]);
        let p = square.evaluate(0.125); // halfway along first edge
        assert!((p - c(0.5, 0.0)).norm() < 1e-10);
    }
}
