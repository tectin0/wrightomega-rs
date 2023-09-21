use num_complex::{Complex, ComplexFloat};
use num_traits::{FloatConst, One, Zero};

use crate::T;

type C = Complex<T>;

/// Value of the Wright omega function.
type Output = C;

/// Last update step in the iterative scheme.
type LastUpdateStep = C;

/// Penultimate residual r_k = z - w_k - log(w_k)
type PenultimateResidual = C;

/// Condition number estimate.
type ConditionNumberEstimate = Option<T>;

/// `wrightomega_ext` is the extended routine for evaluating the wright
/// omega function.
///
///
#[inline]
fn wrightomega_ext(
    z: C,
) -> Option<(
    Output,
    LastUpdateStep,
    PenultimateResidual,
    ConditionNumberEstimate,
)> {
    let mut z = z;

    let mut s = 1.0;
    let x = z.re();
    let y = z.im();
    let mut w = C::zero();

    let pi: T = T::PI();
    let epsilon: T = T::EPSILON;
    let i = C::i();

    let ympi = y - pi;
    let yppi = y + pi;
    let near = 0.01;

    // Test for floating point exceptions
    if x.is_nan() || y.is_nan() {
        return None;
    }
    // signed zeros between branches
    else if x.is_infinite() && x < 0.0 && -pi < y && y <= pi {
        w.re = if y.abs() <= pi / 2.0 { 0.0 } else { -0.0 };
        w.im = if y >= 0.0 { 0.0 } else { -0.0 };

        return Some((w, C::zero(), C::zero(), None));
    }
    // Asymptotic for large z
    else if x.is_infinite() || y.is_infinite() {
        w = C::new(x, y);
        return Some((w, C::zero(), C::zero(), None));
    }

    // Test If exactly on the singular points
    if x == -1.0 && y.abs() == pi {
        return Some((-C::one(), C::zero(), C::zero(), None));
    }

    // Region 1: upper branch point
    // Series about z=-1+Pi*I
    if -2.0 < x && x <= 1.0 && 1.0 < y && y < 2.0 * pi {
        let pz = (2.0 * (z + (-1.0) - i * pi).conj()).sqrt().conj();

        w = -1.0
            + (i + (1.0 / 3.0
                + (-1.0 / 36.0 * i + (1.0 / 270.0 + 1.0 / 4320.0 * i * pz) * pz) * pz)
                * pz)
                * pz;
    }
    // Region 2: lower branch point
    // Series about z=-1-Pi*I
    else if -2.0 < x && x <= 1.0 && -2.0 * pi < y && y < -1.0 {
        let pz = (2.0 * (z + 1.0 + i * pi).conj()).sqrt().conj();
        w = -1.0
            + (-i
                + (1.0 / 3.0 + (1.0 / 36.0 * i + (1.0 / 270.0 - 1.0 / 4320.0 * i * pz) * pz) * pz)
                    * pz)
                * pz;
    }
    // Region 3: between branch cuts
    // Series: About -infinity
    else if x <= -2.0 && -pi < y && y <= pi {
        let pz = z.exp();
        w = (1.0 + (-1.0 + (3.0 / 2.0 + (-8.0 / 3.0 + 125.0 / 24.0 * pz) * pz) * pz) * pz) * pz;
    }
    // Region 4: Mushroom
    // Series about z=1
    else if (-2.0 < x && x <= 1.0 && -1.0 <= y && y <= 1.0)
        || (-2.0 < x && (x - 1.0) * (x - 1.0) + y * y <= pi * pi)
    {
        let pz = z - 1.0;
        w = 1.0 / 2.0
            + 1.0 / 2.0 * z
            + (1.0 / 16.0 + (-1.0 / 192.0 + (-1.0 / 3072.0 + 13.0 / 61440.0 * pz) * pz) * pz)
                * pz
                * pz;
    }
    // Region 5: Top wing
    // Negative log series
    else if x <= -1.05 && pi < y && y - pi <= -0.75 * (x + 1.0) {
        let t = z - C::i() * pi;
        let pz = (-t).ln();

        w = ((1.0 + (-3.0 / 2.0 + 1.0 / 3.0 * pz) * pz) * pz
            + ((-1.0 + 1.0 / 2.0 * pz) * pz + (pz + (-pz + t) * t) * t) * t)
            / (t * t * t);
    }
    // Region 6: Bottom wing
    // Negative log series
    else if x <= -1.05 && 0.75 * (x + 1.0) < y + pi && y + pi <= 0.0 {
        let t = z + C::i() * pi;
        let pz = -t.ln();
        w = ((1.0 + (-3.0 / 2.0 + 1.0 / 3.0 * pz) * pz) * pz
            + ((-1.0 + 1.0 / 2.0 * pz) * pz + (pz + (-pz + t) * t) * t) * t)
            / (t * t * t);
    }
    // Region 7: Everywhere else
    // Series solution about infinity
    else {
        let pz = z.ln();
        w = ((1.0 + (-3.0 / 2.0 + 1.0 / 3.0 * pz) * pz) * pz
            + ((-1.0 + 1.0 / 2.0 * pz) * pz + (pz + (-pz + z) * z) * z) * z)
            / (z * z * z)
    }

    // Regularize if near branch cuts
    // Replaced `fesetround(ROUNDING_MODE)` with `f64::floor()` and `f64::ceil()`
    if x <= -1.0 + near && (ympi.abs() <= near || yppi.abs() <= near) {
        s = -1.0;

        if ympi.abs() <= near {
            let mut ympi = y - pi;

            if ympi <= 0.0 {
                ympi = (y - pi).floor(); // Rounding downward
            }

            z = C::new(x, ympi);
        } else {
            let mut yppi = y + pi;

            if yppi <= 0.0 {
                yppi = (y + pi).ceil(); // Rounding upward
            }

            z = C::new(x, yppi);
        }
    }

    // Iteration one
    w = s * w;
    let mut r = z - s * w - w.ln();
    let mut wp1 = s * w + 1.0;
    let mut e = r / wp1 * (2.0 * wp1 * (wp1 + 2.0 / 3.0 * r) - r)
        / (2.0 * wp1 * (wp1 + 2.0 / 3.0 * r) - 2.0 * r);

    w = w * (1.0 + e);

    // Iteration two
    if ((2.0 * w * w - 8.0 * w - 1.0) * r.abs().powi(4).abs()).abs()
        >= epsilon * 72.0 * wp1.abs().powi(6)
    {
        r = z - s * w - w.ln();
        wp1 = s * w + 1.0;
        e = r / wp1 * (2.0 * wp1 * (wp1 + 2.0 / 3.0 * r) - r)
            / (2.0 * wp1 * (wp1 + 2.0 / 3.0 * r) - 2.0 * r);
        w = w * (1.0 + e);
    }

    // Undo regularization
    w = s * w;

    // Provide condition number estimate
    let cond = (z / (1.0 + w)).norm();

    // Return values
    Some((w, e, r, Some(cond)))
}

/// wrightomega is the routine for evaluating the wright omega function.
/// It is a wrapper around the extended routine `wrightomega_ext`.
///
/// # Examples
/// ```
/// use wrightomega::wrightomega;
///
/// let w = wrightomega(0.0.into()).unwrap();
///
/// assert_eq!(w.re, 0.5671433);
/// assert_eq!(w.im, 0.0);
/// ```
#[inline]
pub fn wrightomega(z: C) -> Option<Output> {
    match wrightomega_ext(z) {
        Some((w, _, _, _)) => Some(w),
        None => None,
    }
}
