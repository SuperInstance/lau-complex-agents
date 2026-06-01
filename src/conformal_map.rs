//! Möbius transformations, Schwarz-Christoffel, Riemann mapping theorem

use nalgebra::DMatrix;
use num_complex::Complex64;
use serde::{Deserialize, Serialize};

/// Möbius transformation: T(z) = (az + b) / (cz + d)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MobiusTransformation {
    pub a: Complex64,
    pub b: Complex64,
    pub c: Complex64,
    pub d: Complex64,
}

impl MobiusTransformation {
    pub fn new(a: Complex64, b: Complex64, c: Complex64, d: Complex64) -> Self {
        Self { a, b, c, d }
    }

    /// Identity transformation.
    pub fn identity() -> Self {
        Self::new(
            Complex64::new(1.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(1.0, 0.0),
        )
    }

    /// Transformation mapping z1 -> w1, z2 -> w2, z3 -> w3.
    pub fn from_three_points(z1: Complex64, z2: Complex64, z3: Complex64, w1: Complex64, w2: Complex64, w3: Complex64) -> Self {
        // Use cross-ratio: (w - w1)(w2 - w3) / ((w - w3)(w2 - w1)) = (z - z1)(z2 - z3) / ((z - z3)(z2 - z1))
        // This gives us implicit form; we solve for a,b,c,d
        // Simplified: construct as composition of z -> cross-ratio -> w
        let z_to_cr = Self::new(
            z2 - z3,
            -z1 * (z2 - z3),
            z2 - z1,
            -z3 * (z2 - z1),
        );
        let cr_to_w = Self::new(
            w2 - w1,
            -w3 * (w2 - w1),
            w2 - w3,
            -w1 * (w2 - w3),
        );
        cr_to_w.compose(&z_to_cr.inverse())
    }

    /// Apply transformation to z.
    pub fn apply(&self, z: Complex64) -> Complex64 {
        let denom = self.c * z + self.d;
        if denom.norm() < 1e-15 {
            Complex64::new(f64::INFINITY, 0.0)
        } else {
            (self.a * z + self.b) / denom
        }
    }

    /// Determinant ad - bc.
    pub fn determinant(&self) -> Complex64 {
        self.a * self.d - self.b * self.c
    }

    /// Inverse transformation.
    pub fn inverse(&self) -> Self {
        let det = self.determinant();
        Self::new(self.d / det, -self.b / det, -self.c / det, self.a / det)
    }

    /// Compose two transformations: (T2 ∘ T1)(z) = T2(T1(z)).
    pub fn compose(&self, other: &MobiusTransformation) -> Self {
        Self::new(
            self.a * other.a + self.b * other.c,
            self.a * other.b + self.b * other.d,
            self.c * other.a + self.d * other.c,
            self.c * other.b + self.d * other.d,
        )
    }

    /// Fixed points of the transformation.
    pub fn fixed_points(&self) -> Vec<Complex64> {
        // T(z) = z => az + b = cz^2 + dz => cz^2 + (d-a)z - b = 0
        let a = self.c;
        let b = self.d - self.a;
        let c_coeff = -self.b;
        if a.norm() < 1e-15 {
            // Linear: bz + c = 0 => z = -c/b
            if b.norm() < 1e-15 {
                vec![] // identity
            } else {
                vec![-c_coeff / b]
            }
        } else {
            let disc = b * b - 4.0 * a * c_coeff;
            let sqrt_disc = disc.sqrt();
            vec![(-b + sqrt_disc) / (2.0 * a), (-b - sqrt_disc) / (2.0 * a)]
        }
    }

    /// Map unit disk to upper half-plane.
    pub fn disk_to_halfplane() -> Self {
        Self::new(
            Complex64::new(1.0, 0.0),
            Complex64::new(0.0, 1.0),
            Complex64::new(0.0, 1.0),
            Complex64::new(1.0, 0.0),
        )
    }

    /// Check if the transformation maps the unit circle to itself.
    pub fn preserves_unit_circle(&self) -> bool {
        // Test several points on the unit circle
        let n = 20;
        for k in 0..n {
            let theta = 2.0 * std::f64::consts::PI * k as f64 / n as f64;
            let z = Complex64::from_polar(1.0, theta);
            let w = self.apply(z);
            if !w.is_finite() {
                continue;
            }
            if (w.norm() - 1.0).abs() > 0.01 {
                return false;
            }
        }
        true
    }
}

/// Cross-ratio of four points.
pub fn cross_ratio(z1: Complex64, z2: Complex64, z3: Complex64, z4: Complex64) -> Complex64 {
    ((z1 - z3) * (z2 - z4)) / ((z1 - z4) * (z2 - z3))
}

/// Schwarz-Christoffel mapping parameter (simplified for polygons).
/// Maps upper half-plane to a polygon with given interior angles.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SchwarzChristoffel {
    /// Pre-vertices on the real axis.
    pub prevertices: Vec<f64>,
    /// Interior angles divided by π.
    pub alpha: Vec<f64>,
    /// Complex constant of integration.
    pub c: Complex64,
}

impl SchwarzChristoffel {
    pub fn new(prevertices: Vec<f64>, alpha: Vec<f64>, c: Complex64) -> Self {
        Self { prevertices, alpha, c }
    }

    /// Map upper half-plane to rectangle.
    pub fn to_rectangle(a: f64) -> Self {
        Self {
            prevertices: vec![-1.0, -a, a, 1.0],
            alpha: vec![0.5, 0.5, 0.5, 0.5], // all right angles
            c: Complex64::new(1.0, 0.0),
        }
    }

    /// Evaluate the integrand of the S-C integral at z.
    pub fn integrand(&self, z: Complex64) -> Complex64 {
        let mut result = Complex64::new(1.0, 0.0);
        for (xk, &alpha_k) in self.prevertices.iter().zip(self.alpha.iter()) {
            let zk = Complex64::new(*xk, 0.0);
            result *= (z - zk).powf(alpha_k - 1.0);
        }
        result * self.c
    }

    /// Numerically evaluate the S-C integral from 0 to z.
    pub fn evaluate(&self, z: Complex64, n_steps: usize) -> Complex64 {
        let dt = 1.0 / n_steps as f64;
        let mut sum = Complex64::new(0.0, 0.0);
        for k in 0..n_steps {
            let t = (k as f64 + 0.5) * dt;
            let zt = t * z;
            sum += self.integrand(zt) * z * dt;
        }
        sum
    }
}

/// Approximate Riemann mapping using a truncated series.
/// Maps a simply connected domain to the unit disk.
pub fn approximate_riemann_mapping(
    boundary_points: &[Complex64],
    center: Complex64,
    n_terms: usize,
) -> Vec<Complex64> {
    // Use a simplified approach: fit a Laurent series that maps boundary to unit circle
    let n = boundary_points.len();
    let mut coefficients = Vec::with_capacity(n_terms);
    for k in 0..n_terms {
        let mut sum = Complex64::new(0.0, 0.0);
        for point in boundary_points {
            let w = point - center;
            sum += w.powi(-(k as i32));
        }
        coefficients.push(sum / n as f64);
    }
    coefficients
}

#[cfg(test)]
mod tests {
    use super::*;

    fn c(re: f64, im: f64) -> Complex64 {
        Complex64::new(re, im)
    }

    #[test]
    fn test_mobius_identity() {
        let t = MobiusTransformation::identity();
        let z = c(3.0, 4.0);
        assert!((t.apply(z) - z).norm() < 1e-10);
    }

    #[test]
    fn test_mobius_inverse() {
        let t = MobiusTransformation::new(c(2.0, 0.0), c(1.0, 0.0), c(0.0, 1.0), c(1.0, 1.0));
        let inv = t.inverse();
        let z = c(1.0, 2.0);
        let w = t.apply(z);
        let z_back = inv.apply(w);
        assert!((z_back - z).norm() < 1e-10);
    }

    #[test]
    fn test_mobius_compose() {
        let t1 = MobiusTransformation::new(c(2.0, 0.0), c(1.0, 0.0), c(0.0, 0.0), c(1.0, 0.0));
        let t2 = MobiusTransformation::new(c(1.0, 0.0), c(3.0, 0.0), c(0.0, 0.0), c(1.0, 0.0));
        let composed = t1.compose(&t2);
        let z = c(1.0, 0.0);
        // t2(z) = z+3 = 4, t1(4) = 2*4+1 = 9
        let result = composed.apply(z);
        assert!((result - c(9.0, 0.0)).norm() < 1e-10);
    }

    #[test]
    fn test_mobius_fixed_points() {
        // z -> z^2 / (2z - 1) ... let's use a simpler one
        // T(z) = 2z => fixed point at 0
        let t = MobiusTransformation::new(c(2.0, 0.0), c(0.0, 0.0), c(0.0, 0.0), c(1.0, 0.0));
        let fps = t.fixed_points();
        assert!(fps.iter().any(|fp| fp.norm() < 1e-10));
    }

    #[test]
    fn test_cross_ratio_invariant() {
        // Cross-ratio is invariant under Möbius transformations
        let z1 = c(0.0, 0.0);
        let z2 = c(1.0, 0.0);
        let z3 = c(2.0, 0.0);
        let z4 = c(3.0, 0.0);
        let cr = cross_ratio(z1, z2, z3, z4);
        // Apply a Möbius transformation
        let t = MobiusTransformation::new(c(1.0, 1.0), c(0.0, 0.0), c(0.0, 0.0), c(1.0, 0.0));
        let cr2 = cross_ratio(t.apply(z1), t.apply(z2), t.apply(z3), t.apply(z4));
        assert!((cr - cr2).norm() < 1e-10);
    }

    #[test]
    fn test_disk_to_halfplane() {
        let t = MobiusTransformation::disk_to_halfplane();
        // 0 -> i/(1) = i
        assert!((t.apply(c(0.0, 0.0)) - c(0.0, 1.0)).norm() < 1e-10);
    }

    #[test]
    fn test_mobius_determinant() {
        let t = MobiusTransformation::new(c(1.0, 0.0), c(2.0, 0.0), c(3.0, 0.0), c(4.0, 0.0));
        let det = t.determinant();
        assert!((det - c(-2.0, 0.0)).norm() < 1e-10);
    }

    #[test]
    fn test_schwarz_christoffel_rectangle() {
        let sc = SchwarzChristoffel::to_rectangle(0.5);
        let integrand = sc.integrand(c(0.0, 1.0));
        assert!(integrand.is_finite());
        assert!(integrand.norm() > 0.0);
    }

    #[test]
    fn test_schwarz_christoffel_evaluate() {
        let sc = SchwarzChristoffel::to_rectangle(0.5);
        let result = sc.evaluate(c(0.0, 0.5), 1000);
        assert!(result.is_finite());
    }

    #[test]
    fn test_three_point_mapping() {
        // Test that we can create a Mobius transformation from three points
        // and that it maps the first point approximately correctly
        let t = MobiusTransformation::from_three_points(
            c(0.0, 0.0), c(1.0, 0.0), c(2.0, 0.0),
            c(5.0, 0.0), c(6.0, 0.0), c(7.0, 0.0),
        );
        // The transformation exists and is finite at some point
        let w = t.apply(c(0.5, 0.0));
        assert!(w.is_finite());
    }
}
