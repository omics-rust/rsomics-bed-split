use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_bed_split(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-bed-split");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let bed = manifest.join("tests/golden/input.bed");
    let tmp = tempfile::TempDir::new().unwrap();
    let prefix = tmp.path().join("bench").to_string_lossy().into_owned();
    c.bench_function("rsomics-bed-split golden", |b| {
        b.iter(|| {
            let status = Command::new(black_box(bin))
                .args([bed.to_str().unwrap(), "-n", "4", "-p", &prefix])
                .status()
                .unwrap();
            assert!(status.success());
        });
    });
}

criterion_group!(benches, bench_bed_split);
criterion_main!(benches);
