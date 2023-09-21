// benchmark lambertw
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use fastapprox::fast::lambertw;
use rand::Rng;

use wrightomega::wright::wrightomega;

pub fn benchmark_lambertw(c: &mut Criterion) {
    let mut rng = rand::thread_rng();

    c.bench_function("lambertw", |b| {
        b.iter(|| {
            black_box(lambertw(rng.gen_range(0.0..100.0)));
        })
    });

    c.bench_function("wrightomega", |b| {
        b.iter(|| {
            black_box(wrightomega(rng.gen_range(0.0..100.0).into()));
        })
    });
}

criterion_group!(benches, benchmark_lambertw);
criterion_main!(benches);
