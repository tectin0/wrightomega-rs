# WrightOmega

This library is a rewrite of the `wright.c` file from Lawrence, Corless, and Jeffrey [[1]](https://dl.acm.org/doi/10.1145/2168773.2168779) in rust.

The [Wright Omega function](https://en.wikipedia.org/wiki/Wright_omega_function) is defined as the solution to the equation `ω + log(ω) = z`. It is a special case of the [Lambert W function](https://en.wikipedia.org/wiki/Lambert_W_function).

## Usage

The library is available on [crates.io](https://crates.io/crates/wright_omega) and can be included in your project by adding the following to your `Cargo.toml` file:

```toml
[dependencies]
wright_omega = "0.1.1"
```

## Example

```rust
use wright_omega::wright_omega;
use wright_omega::Complex;

fn main() {
    let z = Complex::new(0.0, 0.0);
    let omega = wright_omega(z);
    println!("Wright Omega of {} is {}", z, omega);
}
```

> \> Wright Omega of 0+0i is 0.5671433+0i

## Features

The library supports both `f32` and `f64` types. The `f32` feature is enabled by default. To use `f64` instead, add the following to your `Cargo.toml` file:

```toml
[dependencies]
wright_omega = { version = "0.1.0", features = ["f64"] }
```

## Other Implementations

[Scipy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.wrightomega.html) uses the same implementation for Python.
[Matlab](https://www.mathworks.com/help/symbolic/sym.wrightomega.html) also has a built in implementation.

# References

[1] Lawrence, Corless, and Jeffrey, “Algorithm 917: Complex Double-Precision Evaluation of the Wright Function.” ACM Transactions on Mathematical Software, 2012. [DOI:10.1145/2168773.2168779.](https://dl.acm.org/doi/10.1145/2168773.2168779)

```

```
