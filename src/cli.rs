use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_bed_split::{Algorithm, split};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(
    name = "rsomics-bed-split",
    version,
    about = "Split a BED file into N approximately equal parts — Rust port of bedtools split",
    long_about = None,
    disable_help_flag = true
)]
pub struct Cli {
    /// Input BED file.
    #[arg(short = 'i', long = "input", value_name = "FILE")]
    pub input: PathBuf,

    /// Number of output files.
    #[arg(short = 'n', long = "number", value_name = "INT")]
    pub number: usize,

    /// Output file prefix (files named <prefix>.NNNNN.bed).
    #[arg(
        short = 'p',
        long = "prefix",
        value_name = "PREFIX",
        default_value = "split"
    )]
    pub prefix: String,

    /// Splitting algorithm: `size` (equal bases, default) or `simple` (equal records).
    #[arg(
        short = 'a',
        long = "algorithm",
        value_name = "ALG",
        default_value = "size"
    )]
    pub algorithm: String,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }

    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let alg = match self.algorithm.as_str() {
            "size" => Algorithm::Size,
            "simple" => Algorithm::Simple,
            other => {
                return Err(RsomicsError::InvalidInput(format!(
                    "unknown algorithm '{other}'; expected 'size' or 'simple'"
                )));
            }
        };

        split(&self.input, self.number, &self.prefix, alg)?;
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
    tagline: "Split a BED file into N approximately equal parts — Rust port of bedtools split.",
    origin: Some(Origin {
        upstream: "bedtools split",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &["-i <in.bed> -n <N> [-p prefix] [-a size|simple]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('i'),
                long: "input",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("PathBuf"),
                required: true,
                default: None,
                description: "Input BED file.",
                why_default: None,
            },
            FlagSpec {
                short: Some('n'),
                long: "number",
                aliases: &[],
                value: Some("<int>"),
                type_hint: Some("usize"),
                required: true,
                default: None,
                description: "Number of output files to create.",
                why_default: None,
            },
            FlagSpec {
                short: Some('p'),
                long: "prefix",
                aliases: &[],
                value: Some("<str>"),
                type_hint: Some("String"),
                required: false,
                default: Some("split"),
                description: "Output file prefix; files are named <prefix>.NNNNN.bed.",
                why_default: None,
            },
            FlagSpec {
                short: Some('a'),
                long: "algorithm",
                aliases: &[],
                value: Some("<str>"),
                type_hint: Some("String"),
                required: false,
                default: Some("size"),
                description: "size (equal bases, default) or simple (equal records, round-robin).",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "Split peaks.bed into 4 equal-base files",
            command: "rsomics-bed-split -i peaks.bed -n 4 -p out/chunk",
        },
        Example {
            description: "Split by record count (round-robin)",
            command: "rsomics-bed-split -i peaks.bed -n 4 -p out/chunk -a simple",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
