//! Compatibility tests: output identical to bedtools split v2.31.1.
//!
//! Golden fixtures (tests/golden/input.bed) were generated with a seed that
//! produces 200 records with no duplicate interval lengths.  For the `size`
//! algorithm, equal-length records are broken differently by C++ introsort
//! (bedtools) and Rust pdqsort; having no ties makes per-file byte comparison
//! valid.
//!
//! Golden fixtures generated with:
//!   bedtools split -i tests/golden/input.bed -n 4 -p /tmp/size
//!   bedtools split -i tests/golden/input.bed -n 4 -p /tmp/simple -a simple
//!
//! Live bedtools check skipped when bedtools is absent or not v2.31.x.

use std::path::{Path, PathBuf};
use std::process::Command;

use rsomics_bed_split::{Algorithm, split};

fn golden(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn run_split(input: &str, n: usize, algorithm: Algorithm, prefix: &str) {
    let input_path = golden(input);
    split(&input_path, n, prefix, algorithm).expect("split failed");
}

fn read_file(path: &str) -> String {
    std::fs::read_to_string(path).unwrap_or_else(|_| panic!("missing file: {path}"))
}

#[test]
fn size_n4_golden() {
    let dir = tempfile::Builder::new()
        .prefix("rsomics-bed-split-size")
        .tempdir()
        .expect("tempdir");
    let prefix = dir.path().join("out").to_string_lossy().into_owned();
    run_split("input.bed", 4, Algorithm::Size, &prefix);

    for i in 1..=4usize {
        let expected_path = golden(&format!("size_n4.{i:05}.bed"));
        let actual_path = format!("{prefix}.{i:05}.bed");
        let expected = std::fs::read_to_string(&expected_path)
            .unwrap_or_else(|_| panic!("golden missing: {}", expected_path.display()));
        let actual = read_file(&actual_path);
        assert_eq!(actual, expected, "size n=4 file {i:05} mismatch");
    }
}

#[test]
fn simple_n4_golden() {
    let dir = tempfile::Builder::new()
        .prefix("rsomics-bed-split-simple")
        .tempdir()
        .expect("tempdir");
    let prefix = dir.path().join("out").to_string_lossy().into_owned();
    run_split("input.bed", 4, Algorithm::Simple, &prefix);

    for i in 1..=4usize {
        let expected_path = golden(&format!("simple_n4.{i:05}.bed"));
        let actual_path = format!("{prefix}.{i:05}.bed");
        let expected = std::fs::read_to_string(&expected_path)
            .unwrap_or_else(|_| panic!("golden missing: {}", expected_path.display()));
        let actual = read_file(&actual_path);
        assert_eq!(actual, expected, "simple n=4 file {i:05} mismatch");
    }
}

/// Live check against the installed bedtools binary.
///
/// Invariant verified: per-file base-count and record-count match bedtools
/// exactly, and the union of all output records across both tools is
/// identical.  Per-file record order is not required to match: when many
/// records share the same length the `size` algorithm's tie-breaking order
/// is implementation-defined (C++ introsort vs Rust pdqsort differ).
///
/// Skipped gracefully when bedtools is absent or not v2.31.x.
#[test]
fn live_bedtools_compat_size() {
    let bt = match which_bedtools() {
        Some(p) => p,
        None => {
            eprintln!("SKIP live_bedtools_compat_size: bedtools not found");
            return;
        }
    };
    if !bedtools_version_ok(&bt) {
        eprintln!("SKIP live_bedtools_compat_size: bedtools version mismatch (want v2.31.x)");
        return;
    }

    let input = golden("input.bed");

    let dir_bt = tempfile::Builder::new()
        .prefix("rsomics-bed-split-live-bt")
        .tempdir()
        .expect("tempdir");
    let bt_prefix = dir_bt.path().join("bt").to_string_lossy().into_owned();

    let status = Command::new(&bt)
        .args(["split", "-i"])
        .arg(&input)
        .args(["-n", "4", "-p", &bt_prefix])
        .status()
        .expect("bedtools split failed");
    assert!(status.success(), "bedtools split exited non-zero");

    let dir_ours = tempfile::Builder::new()
        .prefix("rsomics-bed-split-live-ours")
        .tempdir()
        .expect("tempdir");
    let our_prefix = dir_ours.path().join("out").to_string_lossy().into_owned();

    let our_summary = split(&input, 4, &our_prefix, Algorithm::Size).expect("our split failed");

    // Verify per-file base counts and record counts match exactly.
    for i in 1..=4usize {
        let bt_path = format!("{bt_prefix}.{i:05}.bed");
        let bt_content = read_file(&bt_path);
        let bt_bases: u64 = bt_content
            .lines()
            .map(|l| {
                let mut c = l.splitn(4, '\t');
                let s: u64 = c.nth(1).unwrap_or("0").parse().unwrap_or(0);
                let e: u64 = c.next().unwrap_or("0").parse().unwrap_or(0);
                e.saturating_sub(s)
            })
            .sum();
        let bt_records = bt_content.lines().count() as u64;
        let (_, our_bases, our_records) = &our_summary[i - 1];
        assert_eq!(
            *our_bases, bt_bases,
            "size file {i:05}: base count mismatch (ours={our_bases}, bedtools={bt_bases})"
        );
        assert_eq!(
            *our_records, bt_records,
            "size file {i:05}: record count mismatch"
        );
    }

    // Verify the union of all records across all files is identical.
    let mut bt_all: Vec<String> = (1..=4)
        .flat_map(|i| {
            read_file(&format!("{bt_prefix}.{i:05}.bed"))
                .lines()
                .map(String::from)
                .collect::<Vec<_>>()
        })
        .collect();
    let mut our_all: Vec<String> = (1..=4)
        .flat_map(|i| {
            read_file(&format!("{our_prefix}.{i:05}.bed"))
                .lines()
                .map(String::from)
                .collect::<Vec<_>>()
        })
        .collect();
    bt_all.sort_unstable();
    our_all.sort_unstable();
    assert_eq!(our_all, bt_all, "live size: combined record sets differ");
}

#[test]
fn live_bedtools_compat_simple() {
    let bt = match which_bedtools() {
        Some(p) => p,
        None => {
            eprintln!("SKIP live_bedtools_compat_simple: bedtools not found");
            return;
        }
    };
    if !bedtools_version_ok(&bt) {
        eprintln!("SKIP live_bedtools_compat_simple: bedtools version mismatch (want v2.31.x)");
        return;
    }

    let input = golden("input.bed");

    let dir_bt = tempfile::Builder::new()
        .prefix("rsomics-bed-split-live-bt-simple")
        .tempdir()
        .expect("tempdir");
    let bt_prefix = dir_bt.path().join("bt").to_string_lossy().into_owned();

    let status = Command::new(&bt)
        .args(["split", "-i"])
        .arg(&input)
        .args(["-n", "4", "-p", &bt_prefix, "-a", "simple"])
        .status()
        .expect("bedtools split -a simple failed");
    assert!(status.success(), "bedtools split -a simple exited non-zero");

    let dir_ours = tempfile::Builder::new()
        .prefix("rsomics-bed-split-live-ours-simple")
        .tempdir()
        .expect("tempdir");
    let our_prefix = dir_ours.path().join("out").to_string_lossy().into_owned();

    // simple is deterministic (round-robin in input order): byte-identical is required.
    split(&input, 4, &our_prefix, Algorithm::Simple).expect("our split failed");

    for i in 1..=4usize {
        let bt_path = format!("{bt_prefix}.{i:05}.bed");
        let our_path = format!("{our_prefix}.{i:05}.bed");
        let bt_content = read_file(&bt_path);
        let our_content = read_file(&our_path);
        assert_eq!(our_content, bt_content, "live simple file {i:05} mismatch");
    }
}

fn which_bedtools() -> Option<String> {
    let out = Command::new("which").arg("bedtools").output().ok()?;
    if out.status.success() {
        let p = String::from_utf8(out.stdout).ok()?.trim().to_string();
        if p.is_empty() { None } else { Some(p) }
    } else {
        None
    }
}

fn bedtools_version_ok(bt: &str) -> bool {
    let out = Command::new(bt)
        .arg("--version")
        .output()
        .expect("bedtools --version failed");
    let v = String::from_utf8_lossy(&out.stdout);
    v.contains("v2.31")
}
