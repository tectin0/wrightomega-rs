pub mod wright;

pub use wright::wrightomega;

#[cfg(feature = "f32")]
type T = f32;
#[cfg(feature = "f64")]
type T = f64;

#[cfg(test)]
mod test {
    use fastapprox::fast::lambertw;

    use crate::wrightomega;
    use crate::T;

    #[test]
    fn compare_with_lambertw() {
        for x in 0..=100 {
            let z = num_complex::Complex::<T>::new(x as T / 5.0, 0.0);
            let w = wrightomega(z).unwrap();
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
}
