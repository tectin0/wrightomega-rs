use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rand::Rng;

use wright_omega::wright_omega;

/// Benchmark the Wright Omega function implementation.
pub fn benchmark(c: &mut Criterion) {
    let mut rng = rand::thread_rng();

    c.bench_function("wrightomega", |b| {
        b.iter(|| {
            black_box(wright_omega(rng.gen_range(0.0..100.0).into()));
        })
    });
}

criterion_group!(benches, benchmark);
criterion_main!(benches);
