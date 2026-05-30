//! Split a BED file into N approximately equal parts — bedtools split equivalent.
//!
//! Two algorithms:
//! - `size` (default): greedy bin-packing by interval length, balanced total bp
//! - `simple`: round-robin by input order, balanced record count
//!
//! Output: `<prefix>.NNNNN.bed`; summary `<file>\t<bases>\t<records>` to stderr.

use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Algorithm {
    Size,
    Simple,
}

pub fn split(
    input: &Path,
    n: usize,
    prefix: &str,
    algorithm: Algorithm,
) -> Result<Vec<(String, u64, u64)>> {
    if n == 0 {
        return Err(RsomicsError::InvalidInput(
            "number of output files must be >= 1".into(),
        ));
    }

    let records = read_records(input)?;

    let write_order = match algorithm {
        Algorithm::Size => assign_size(&records, n),
        Algorithm::Simple => assign_simple(records.len(), n),
    };

    write_files(&records, &write_order, n, prefix)
}

struct Record {
    line: String,
    len: u64,
}

fn read_records(path: &Path) -> Result<Vec<Record>> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let reader = BufReader::new(file);
    let mut records = Vec::new();

    for raw in reader.lines() {
        let line = raw.map_err(RsomicsError::Io)?;
        let trimmed = line.trim_end_matches('\r');
        if trimmed.is_empty()
            || trimmed.starts_with('#')
            || trimmed.starts_with("track")
            || trimmed.starts_with("browser")
        {
            continue;
        }
        let len = parse_len(trimmed)?;
        records.push(Record {
            line: trimmed.to_string(),
            len,
        });
    }

    Ok(records)
}

fn parse_len(line: &str) -> Result<u64> {
    let mut cols = line.splitn(4, '\t');
    let _chrom = cols
        .next()
        .ok_or_else(|| RsomicsError::InvalidInput(format!("bad BED line: {line}")))?;
    let start_s = cols
        .next()
        .ok_or_else(|| RsomicsError::InvalidInput(format!("bad BED line (no start): {line}")))?;
    let end_s = cols
        .next()
        .ok_or_else(|| RsomicsError::InvalidInput(format!("bad BED line (no end): {line}")))?;
    let start: u64 = start_s
        .parse()
        .map_err(|e| RsomicsError::InvalidInput(format!("start parse '{start_s}': {e}")))?;
    let end: u64 = end_s
        .parse()
        .map_err(|e| RsomicsError::InvalidInput(format!("end parse '{end_s}': {e}")))?;
    Ok(end.saturating_sub(start))
}

/// Records appear in descending-length order within each bin (bedtools' processing order).
fn assign_size(records: &[Record], n: usize) -> Vec<Vec<usize>> {
    let mut order: Vec<usize> = (0..records.len()).collect();
    order.sort_unstable_by(|&a, &b| records[b].len.cmp(&records[a].len));

    let mut heap: BinaryHeap<Reverse<(u64, usize)>> = (0..n).map(|i| Reverse((0, i))).collect();

    let mut bin_records: Vec<Vec<usize>> = vec![Vec::new(); n];

    for rec_idx in order {
        let Reverse((bases, bin)) = heap.pop().unwrap();
        bin_records[bin].push(rec_idx);
        heap.push(Reverse((bases + records[rec_idx].len, bin)));
    }

    bin_records
}

fn assign_simple(record_count: usize, n: usize) -> Vec<Vec<usize>> {
    let mut bins: Vec<Vec<usize>> = vec![Vec::new(); n];
    for i in 0..record_count {
        bins[i % n].push(i);
    }
    bins
}

fn write_files(
    records: &[Record],
    write_order: &[Vec<usize>],
    n: usize,
    prefix: &str,
) -> Result<Vec<(String, u64, u64)>> {
    let filenames: Vec<String> = (1..=n).map(|i| format!("{prefix}.{i:05}.bed")).collect();

    let mut summary: Vec<(String, u64, u64)> = Vec::with_capacity(n);

    for (bin, name) in filenames.iter().enumerate() {
        let f = File::create(name)
            .map_err(|e| RsomicsError::InvalidInput(format!("Cannot open \"{name}\". {e}")))?;
        let mut w = BufWriter::new(f);
        let mut bases: u64 = 0;
        let mut count: u64 = 0;
        for &idx in &write_order[bin] {
            writeln!(w, "{}", records[idx].line).map_err(RsomicsError::Io)?;
            bases += records[idx].len;
            count += 1;
        }
        w.flush().map_err(RsomicsError::Io)?;
        summary.push((name.clone(), bases, count));
    }

    let stderr = std::io::stderr();
    let mut err = stderr.lock();
    for (name, b, c) in &summary {
        writeln!(err, "{name}\t{b}\t{c}").map_err(RsomicsError::Io)?;
    }

    Ok(summary)
}
