use wright_omega::wright_omega;
use wright_omega::Complex;

fn main() {
    let z = Complex::new(0.0, 0.0);
    let omega = wright_omega(z).unwrap();
    println!("Wright Omega of {} is {}", z, omega);
}
