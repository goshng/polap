// Version: v0.1.0
use crate::cli::{AssembleArgs, Cli};
use crate::runner::CommandRunner;
use anyhow::Result;
use std::path::Path;
use tracing::info;

pub fn run(cli: &Cli, a: &AssembleArgs) -> Result<()> {
    if a.help {
        // Build the here-doc style text (you can mirror your Bash snippet)
        let txt = format!(
            r#"Name:
  polap assemble — assemble organelle genome sequences

Synopsis:
  polap assemble [options]

Description:
  polap assemble runs the organelle-genome assembly using Flye (stage1) and seed-based refinement (stage2).
  This mirrors the legacy Bash menus (assemble1/assemble2) and will be migrated over time.

Options:
  -i INT
    index of the source of an organelle-genome assembly [default: {inum}]
  
  -j INT
    index of the target organelle-genome assembly [default: {jnum}]
  
  -w INT
    minimum mapping length for read selection [default: {single_min}]

  -c INT
    the coverage option [default: {coverage}]

  -t INT
    the number of CPU cores [default: {threads}]

  --polap-reads [default: {polap_reads}]
    uses the POLAP read selection (instead of ptGAUL)

  --coverage-check [default: {coverage_check}]
    verify depth uniformity before proceeding

  -l FASTQ
    reads data file [default: {lr}]

  -o OUTDIR
    output folder [default: {outdir}]

  -m INT
    minimum read length [default: {min_read_length}]

Inputs:
  {inputs}

Outputs:
  {outputs}

Examples:
  Get organelle genome sequences using minimum mapping length 6000:
    polap assemble -l l.fq -o outdir -w 6000

  Get organelle genome sequences using seed contig index 0 and target organelle-genome assembly index 1:
    polap assemble -i 0 -j 1 -l l.fq -o outdir

  Get organelle genome sequences using polap read selection:
    polap assemble --polap-reads -l l.fq -o outdir

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024–2025)

Author:
  Sang Chul Choi
"#,
            inum = a.inum,
            jnum = a.jnum,
            single_min = a.single_min,
            coverage = a.coverage,
            threads = a.threads,
            polap_reads = if a.polap_reads { "on" } else { "off" },
            coverage_check = if a.coverage_check { "on" } else { "off" },
            lr = a
                .long_reads
                .as_ref()
                .map(|p| p.to_string_lossy())
                .unwrap_or_else(|| "<none>".into()),
            outdir = a.outdir.to_string_lossy(),
            min_read_length = a.min_read_length,
            inputs = "(list your expected input files here)",
            outputs = "(list your produced outputs here)",
        );

        // Page it as a man page
        crate::man::page_text_as_man("polap-assemble", "1", &txt, &cli.version_string())?;
        return Ok(());
    }

    // Normal execution:
    let _runner = CommandRunner::new(cli);
    // TODO: call your existing scripts
    Ok(())
}
