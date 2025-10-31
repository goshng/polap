// crates/polap/src/subcommands/test.rs
// Version: v0.4.2 — logging sampler + conda script runners (bash/R/python)

use anyhow::{Context, Result};
use std::fs;
use std::path::{Path, PathBuf};
use tracing::{debug, error, info, trace, warn};

use crate::cli::Cli;
use crate::runner::CommandRunner;

pub fn run_log(_cli: &Cli) -> Result<()> {
    // Emit logs at all levels to verify -q / -v and file logging
    trace!("trace message (level 3)");
    debug!("debug message (level 2)");
    info!("info message (level 1)");
    warn!("warn message");
    error!("error message");

    // Lines mimicking your desired style
    info!("[_polap_lib_readassemble-miniasm@polap-lib-readassemble.sh:1153] mtDNA assembly: Anthoceros_angustus-0/mtseed/mt.3.gfa");
    info!("[_run_polap_miniassemble@polap-cmd-miniassemble.sh:218] output plastid assembly graph: SRR9696346.pt.gfa");
    info!("[_run_polap_miniassemble@polap-cmd-miniassemble.sh:222] output mitochondrial assembly graph: SRR9696346.mt.gfa");
    info!("[_run_polap_extract@polap-cmd-extract.sh:153] ptDNA: Anthoceros_angustus-0/extract/oatk.pltd.ctg.fasta");
    info!("[_run_polap_extract@polap-cmd-extract.sh:154] mtDNA: Anthoceros_angustus-0/extract/oatk.mito.ctg.fasta");

    info!("Log test complete. Check screen (-q/-v...) and file log.");
    Ok(())
}

pub fn run_bash(cli: &Cli) -> Result<()> {
    let p = ensure_sample(cli, "test/sample.sh", SAMPLE_BASH, 0o755)?;
    let r = CommandRunner::new(cli);
    r.run_conda("polap", "bash", &[p.to_str().unwrap()])
}

pub fn run_r(cli: &Cli) -> Result<()> {
    let p = ensure_sample(cli, "test/sample.R", SAMPLE_R, 0o755)?;
    let r = CommandRunner::new(cli);
    r.run_conda("polap", "Rscript", &[p.to_str().unwrap()])
}

pub fn run_python(cli: &Cli) -> Result<()> {
    let p = ensure_sample(cli, "test/sample.py", SAMPLE_PY, 0o755)?;
    let r = CommandRunner::new(cli);
    r.run_conda("polap", "python", &[p.to_str().unwrap()])
}

fn ensure_sample(_cli: &Cli, rel: &str, content: &str, mode: u32) -> Result<PathBuf> {
    let lib = std::env::var("_POLAPLIB_DIR").unwrap_or_else(|_| ".".to_string());
    let p = PathBuf::from(lib).join(rel);
    if let Some(parent) = p.parent() {
        fs::create_dir_all(parent)?;
    }
    if !p.exists() {
        fs::write(&p, content.as_bytes()).with_context(|| format!("write {}", p.display()))?;
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            let mut perm = fs::metadata(&p)?.permissions();
            perm.set_mode(mode);
            fs::set_permissions(&p, perm)?;
        }
        info!("Created sample script: {}", p.display());
    }
    Ok(p)
}

// ── Sample scripts ────────────────────────────────────────────────────────────
const SAMPLE_BASH: &str = r#"#!/usr/bin/env bash
set -euo pipefail
echo "[BASH] hello from sample.sh"
echo "[BASH] _POLAPLIB_DIR=${_POLAPLIB_DIR:-}"
echo "[BASH] date: $(date)"
"#;

const SAMPLE_R: &str = r#"#!/usr/bin/env Rscript
cat("[R] hello from sample.R\n")
cat(sprintf("[R] _POLAPLIB_DIR=%s\n", Sys.getenv("_POLAPLIB_DIR")))
cat(sprintf("[R] date: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
"#;

const SAMPLE_PY: &str = r#"#!/usr/bin/env python
import os, datetime
print("[PY] hello from sample.py")
print("[PY] _POLAPLIB_DIR={}".format(os.environ.get("_POLAPLIB_DIR","")))
print("[PY] date: {}".format(datetime.datetime.utcnow().isoformat()))
"#;
