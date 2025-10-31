// crates/polap/src/cli.rs
// Version: v0.4.0 — Presets + Logging controls + Test subcommands.

use anyhow::{anyhow, Result};
use clap::{ArgAction, Args, Parser, Subcommand, ValueEnum, ValueHint};
use std::path::{Path, PathBuf};
use std::time::Instant;

use crate::preset::{finalize_options, save_preset, FinalOpts};

/// Timestamp format for log lines
#[derive(Copy, Clone, Debug, ValueEnum)]
pub enum LogTs {
    Off,
    Date,
    Time,
    Datetime,
}

/// Top-level CLI
#[derive(Parser, Debug)]
#[command(
    name = "polap",
    about = "POLAP – Plant organelle DNA long-read assembly pipeline (Rust CLI)",
    version = env!("CARGO_PKG_VERSION"),
    disable_help_subcommand = true
)]
pub struct Top {
    // ── Output / IO ────────────────────────────────────────────────────────────
    #[arg(short='o', long, value_hint=ValueHint::DirPath, default_value="o", global=true)]
    pub outdir: PathBuf,
    #[arg(short='l', long="long-reads", value_hint=ValueHint::FilePath, global=true)]
    pub long_reads: Option<PathBuf>,
    #[arg(short='a', long="short-read1", value_hint=ValueHint::FilePath, global=true)]
    pub sr1: Option<PathBuf>,
    #[arg(short='b', long="short-read2", value_hint=ValueHint::FilePath, global=true)]
    pub sr2: Option<PathBuf>,

    // ── Tuning / resources ─────────────────────────────────────────────────────
    #[arg(short = 't', long, global = true)]
    pub threads: Option<usize>,
    #[arg(short = 'c', long, global = true)]
    pub coverage: Option<u64>,
    #[arg(short = 'w', long = "single-min", global = true)]
    pub single_min: Option<u64>,
    #[arg(short = 'r', long = "pair-min", global = true)]
    pub pair_min: Option<u64>,
    #[arg(short = 'x', long = "bridge-min", global = true)]
    pub bridge_min: Option<u64>,
    #[arg(short = 'g', long, global = true)]
    pub genomesize: Option<String>,
    #[arg(long = "min-read-length", global = true)]
    pub min_read_length: Option<u64>,
    #[arg(long = "polap-reads", global = true)]
    pub polap_reads: Option<bool>,
    #[arg(long = "coverage-check", global = true)]
    pub coverage_check: Option<bool>,

    // ── Preset system ─────────────────────────────────────────────────────────
    #[arg(long = "preset", global = true)]
    pub preset: Option<String>,
    #[arg(long="preset-path", value_hint=ValueHint::FilePath, global=true)]
    pub preset_path: Option<PathBuf>,
    #[arg(long="preset-dir", value_hint=ValueHint::DirPath, default_value_os_t=default_preset_dir(), global=true)]
    pub preset_dir: PathBuf,
    #[arg(long = "preset-save", default_value_t = false, global = true)]
    pub preset_save: bool,

    // ── Logging controls ──────────────────────────────────────────────────────
    /// Increase screen verbosity (-v, -vv, -vvv). File always logs at TRACE.
    #[arg(short='v', long, action=ArgAction::Count)]
    pub verbose: u8,
    /// Quiet: no screen output at all (file logging still on).
    #[arg(short = 'q', long, default_value_t = false)]
    pub quiet: bool,
    /// Timestamp on each log line (screen + file)
    #[arg(long="log-ts", value_enum, default_value_t=LogTs::Datetime)]
    pub log_ts: LogTs,
    /// Directly set the log file path (otherwise: <outdir>/polap.log)
    #[arg(long="log-file", value_hint=ValueHint::FilePath)]
    pub log_file: Option<PathBuf>,

    // ── Misc ──────────────────────────────────────────────────────────────────
    #[arg(long, default_value_t = false)]
    pub clock: bool,
    #[arg(long, default_value_t = false, global = true)]
    pub dry: bool,
    #[arg(long, default_value_t = false, global = true)]
    pub plastid: bool,
    #[arg(long, default_value_t = false, global = true)]
    pub noncoding: bool,

    // Subcommands
    #[command(subcommand)]
    pub cmd: Option<Cmd>,
}

#[derive(Subcommand, Debug)]
pub enum Cmd {
    #[command(
        name = "assemble",
        about = "Assemble organelle genomes from long/short reads",
        disable_help_flag = true,
        disable_help_subcommand = true
    )]
    Assemble(AssembleArgs),

    PreparePolishing(GenericIo),
    Polish(GenericIo),

    /// Manage presets
    #[command(name = "preset")]
    Preset {
        #[command(subcommand)]
        action: PresetAction,
    },

    /// Test helpers
    #[command(name = "test")]
    Test {
        #[command(subcommand)]
        action: TestAction,
    },

    #[command(hide = true)]
    Plugin {
        name: String,
        #[arg(trailing_var_arg = true)]
        args: Vec<String>,
    },
}

// ── Preset actions ────────────────────────────────────────────────────────────
#[derive(Subcommand, Debug)]
pub enum PresetAction {
    List,
    Show {
        name_or_path: String,
    },
    Edit {
        name_or_path: String,
    },
    /// Save a preset (Defaults < existing preset < CLI); `--print` prints to stdout only.
    Save {
        name_or_path: String,
        #[arg(long)]
        print: bool,
    },
}

// ── Test actions ──────────────────────────────────────────────────────────────
#[derive(Subcommand, Debug)]
pub enum TestAction {
    /// Emit sample logs at all levels to screen/file
    Log,
    /// Run a Bash script under conda env "polap" (auto-creates sample if missing)
    Bash,
    /// Run an R script under conda env "polap" (auto-creates sample if missing)
    R,
    /// Run a Python script under conda env "polap" (auto-creates sample if missing)
    Python,
}

#[derive(Args, Debug)]
pub struct AssembleArgs {
    #[arg(long = "help", default_value_t = false)]
    pub help: bool,

    #[arg(long, value_hint=ValueHint::FilePath)]
    pub long_reads: Option<PathBuf>,
    #[arg(long="short-read1", value_hint=ValueHint::FilePath)]
    pub sr1: Option<PathBuf>,
    #[arg(long="short-read2", value_hint=ValueHint::FilePath)]
    pub sr2: Option<PathBuf>,

    #[arg(long, default_value_t = 3000)]
    pub single_min: u64,
    #[arg(long, default_value_t = 50)]
    pub coverage: u64,
    #[arg(long, default_value_t = 3000)]
    pub min_read_length: u64,
    #[arg(long, default_value_t=num_cpus::get())]
    pub threads: usize,
    #[arg(long, default_value_t = false)]
    pub polap_reads: bool,
    #[arg(long, default_value_t = true)]
    pub coverage_check: bool,

    #[arg(short='o', long, default_value="o", value_hint=ValueHint::DirPath)]
    pub outdir: PathBuf,

    #[arg(long, default_value_t = 0)]
    pub inum: i32,
    #[arg(long, default_value_t = 1)]
    pub jnum: i32,
}

#[derive(Args, Debug)]
pub struct GenericIo {
    #[arg(long, value_hint=ValueHint::FilePath)]
    pub infile: Option<PathBuf>,
    #[arg(long, value_hint=ValueHint::FilePath)]
    pub outfile: Option<PathBuf>,
}

/// Resolved runtime context
pub struct Cli {
    pub args: Top,
    pub opts: FinalOpts,
    pub preset_path: Option<PathBuf>,
    pub out_log: PathBuf,
    pub start: Instant,
    pub version: &'static str,
}

impl Cli {
    pub fn parse_and_init() -> Result<Self> {
        let args = Top::parse();

        // Resolve preset path
        let preset_path = resolve_preset_path(
            args.preset.as_deref(),
            args.preset_path.as_ref(),
            &args.preset_dir,
        );

        // Merge: Defaults < Preset < CLI
        let opts = finalize_options(&args, preset_path.as_deref())?;

        // Direct log file resolution
        let log_file = if let Some(p) = &args.log_file {
            p.clone()
        } else {
            opts.outdir.join("polap.log")
        };
        // Ensure dirs (outdir/log/tmp) + log file dir
        std::fs::create_dir_all(&opts.outdir.join("log"))?;
        std::fs::create_dir_all(&opts.outdir.join("tmp"))?;
        if let Some(parent) = log_file.parent() {
            std::fs::create_dir_all(parent)?;
        } else {
            return Err(anyhow!(
                "Invalid --log-file (no parent dir): {}",
                log_file.display()
            ));
        }

        Ok(Self {
            out_log: log_file,
            preset_path,
            args,
            opts,
            start: Instant::now(),
            version: env!("CARGO_PKG_VERSION"),
        })
    }

    pub fn version_string(&self) -> String {
        format!("v{}", self.version)
    }
}

pub fn default_preset_dir() -> PathBuf {
    dirs::home_dir()
        .unwrap_or_else(|| Path::new(".").to_path_buf())
        .join(".polap")
        .join("presets")
}

fn resolve_preset_path(
    preset: Option<&str>,
    preset_path: Option<&PathBuf>,
    preset_dir: &Path,
) -> Option<PathBuf> {
    if let Some(pp) = preset_path {
        return Some(pp.clone());
    }
    if let Some(name) = preset {
        return Some(preset_dir.join(format!("{name}.yaml")));
    }
    None
}
