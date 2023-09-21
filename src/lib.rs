//! This crate provides an implementation of the Wright Omega function.
//! The Wright Omega function is defined as the solution to the equation
//! `ω + log(ω) = z`. This crate provides two functions for evaluating the
//! Wright Omega function: `wright_omega` and `wright_omega_ext`.
//! The `wright_omega` function returns the value of the Wright Omega function
//! for a given complex number `z`. The `wright_omega_ext` function returns
//! the value of the Wright Omega function, as well as the last update step,
//! the penultimate residual, and an estimate of the condition number.

pub mod wright_omega;

pub use wright_omega::wright_omega;
pub use wright_omega::wright_omega_ext;

#[cfg(feature = "f32")]
type T = f32;
#[cfg(feature = "f64")]
type T = f64;

#[cfg(test)]
mod test {
    use fastapprox::fast::lambertw;
    use num_traits::FloatConst;
    use num_traits::One;
    use num_traits::Zero;

    use crate::wright_omega;
    use crate::T;

    use num_complex::Complex;

    /// Compare the result of the Wright Omega function with the result of the
    /// Lambert W function. As the Wright Omega function is a special case of the
    /// Lambert W function, the results should be identical for these inputs.
    #[test]
    fn compare_with_lambertw() {
        for x in 0..=100 {
            let z = Complex::<T>::new(x as T / 5.0, 0.0);
            let w = wright_omega(z).unwrap();
            let l = lambertw((x as f32 / 5.0).exp()) as T;

            assert!(
                (w.re - l).abs() < 1e-3,
                "Failed for x = {} with w = {} and l = {}",
                x,
                w.re,
                l
            );
        }
    }

    /// Compare the result of the Wright Omega function with known values.
    /// [Source](https://en.wikipedia.org/wiki/Wright_omega_function#Values)
    #[test]
    fn check_known_values() {
        let values = [
            Complex::zero(),
            Complex::one(),
            Complex::new(-1.0, -<T>::PI()),
            Complex::new(-1.0, <T>::PI()),
            Complex::new(-1.0 / 3.0 + T::from(1.0 / 3.0).ln(), <T>::PI()),
            Complex::new(-1.0 / 3.0 + T::from(1.0 / 3.0).ln(), -<T>::PI()),
        ];

        let expected = [
            Complex::new(0.56714, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(-1.0, 0.0),
            Complex::new(-1.0, 0.0),
            Complex::new(-2.23715, 0.0),
            Complex::new(-1.0 / 3.0, 0.0),
        ];

        for (z, w) in values.iter().zip(expected.iter()) {
            let omega = wright_omega(*z).unwrap();
            assert!(
                (omega.re - w.re).abs() < 1e-3,
                "Failed for z = {} with w = {} and expected = {}",
                z,
                omega.re,
                w.re
            );
        }
    }
}
