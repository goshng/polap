// crates/polap/src/preset.rs
// Version: v0.3.0 – Preset system + commented YAML render + list/show/edit helpers.

use anyhow::{anyhow, Context, Result};
use serde::{Deserialize, Serialize};
use std::{
    fs,
    io::Write,
    path::{Path, PathBuf},
};

use crate::cli::Top;

// ---------- 1) Hard-coded defaults (source-of-truth) ----------
#[derive(Debug, Clone)]
pub struct Defaults;
impl Defaults {
    pub const OUTDIR: &str = "o";
    pub const THREADS: usize = 8;
    pub const COVERAGE: u64 = 50;
    pub const SINGLE_MIN: u64 = 3000;
    pub const PAIR_MIN: u64 = 3000;
    pub const BRIDGE_MIN: u64 = 0;
    pub const MIN_READ_LENGTH: u64 = 3000;
    pub const POLAP_READS: bool = false;
    pub const COVERAGE_CHECK: bool = true;
    pub const GENOMESIZE: &'static str = ""; // empty = unspecified
}

// ---------- 2) Serializable preset file (simple key–value YAML) ----------
#[derive(Debug, Default, Clone, Serialize, Deserialize)]
pub struct Preset {
    pub outdir: Option<PathBuf>,
    pub long_reads: Option<PathBuf>,
    pub sr1: Option<PathBuf>,
    pub sr2: Option<PathBuf>,

    pub threads: Option<usize>,
    pub coverage: Option<u64>,
    pub single_min: Option<u64>,
    pub pair_min: Option<u64>,
    pub bridge_min: Option<u64>,
    pub min_read_length: Option<u64>,
    pub polap_reads: Option<bool>,
    pub coverage_check: Option<bool>,
    pub genomesize: Option<String>,
}

// ---------- 3) Final merged options (concrete types) ----------
#[derive(Debug, Clone, Serialize)]
pub struct FinalOpts {
    pub outdir: PathBuf,
    pub long_reads: Option<PathBuf>,
    pub sr1: Option<PathBuf>,
    pub sr2: Option<PathBuf>,

    pub threads: usize,
    pub coverage: u64,
    pub single_min: u64,
    pub pair_min: u64,
    pub bridge_min: u64,
    pub min_read_length: u64,
    pub polap_reads: bool,
    pub coverage_check: bool,
    pub genomesize: Option<String>,
}

// ---------- 4) Load / Save preset YAML ----------
pub fn load_preset(path: &Path) -> Result<Preset> {
    if !path.exists() || fs::metadata(path)?.len() == 0 {
        return Ok(Preset::default());
    }
    let txt = fs::read_to_string(path).with_context(|| format!("reading preset {:?}", path))?;
    let p: Preset =
        serde_yaml::from_str(&txt).with_context(|| format!("parsing YAML {:?}", path))?;
    Ok(p)
}

pub fn save_preset(path: &Path, opts: &FinalOpts) -> Result<()> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent).with_context(|| format!("create dir {:?}", parent))?;
    }
    let p = Preset {
        outdir: Some(opts.outdir.clone()),
        long_reads: opts.long_reads.clone(),
        sr1: opts.sr1.clone(),
        sr2: opts.sr2.clone(),
        threads: Some(opts.threads),
        coverage: Some(opts.coverage),
        single_min: Some(opts.single_min),
        pair_min: Some(opts.pair_min),
        bridge_min: Some(opts.bridge_min),
        min_read_length: Some(opts.min_read_length),
        polap_reads: Some(opts.polap_reads),
        coverage_check: Some(opts.coverage_check),
        genomesize: opts.genomesize.clone(),
    };
    let yaml = serde_yaml::to_string(&p)?;
    fs::write(path, yaml).with_context(|| format!("writing preset {:?}", path))?;
    Ok(())
}

// ---------- 5) Precedence merge: Defaults < Preset(YAML?) < CLI ----------
pub fn finalize_options(args: &Top, preset_path: Option<&Path>) -> Result<FinalOpts> {
    let preset = if let Some(pp) = preset_path {
        load_preset(pp)?
    } else {
        Preset::default()
    };

    let outdir = args.outdir.clone(); // CLI default "o" always wins (keeps Bash parity)
    let long_reads = pick_opt(args.long_reads.clone(), preset.long_reads);
    let sr1 = pick_opt(args.sr1.clone(), preset.sr1);
    let sr2 = pick_opt(args.sr2.clone(), preset.sr2);

    let threads = pick_val(args.threads, preset.threads, Defaults::THREADS);
    let coverage = pick_val(args.coverage, preset.coverage, Defaults::COVERAGE);
    let single_min = pick_val(args.single_min, preset.single_min, Defaults::SINGLE_MIN);
    let pair_min = pick_val(args.pair_min, preset.pair_min, Defaults::PAIR_MIN);
    let bridge_min = pick_val(args.bridge_min, preset.bridge_min, Defaults::BRIDGE_MIN);
    let min_read_length = pick_val(
        args.min_read_length,
        preset.min_read_length,
        Defaults::MIN_READ_LENGTH,
    );
    let polap_reads = pick_val(args.polap_reads, preset.polap_reads, Defaults::POLAP_READS);
    let coverage_check = pick_val(
        args.coverage_check,
        preset.coverage_check,
        Defaults::COVERAGE_CHECK,
    );

    let genomesize = pick_opt(
        args.genomesize.clone(),
        preset.genomesize.or_else(|| {
            let s = Defaults::GENOMESIZE.trim();
            if s.is_empty() {
                None
            } else {
                Some(s.to_string())
            }
        }),
    );

    Ok(FinalOpts {
        outdir,
        long_reads,
        sr1,
        sr2,
        threads,
        coverage,
        single_min,
        pair_min,
        bridge_min,
        min_read_length,
        polap_reads,
        coverage_check,
        genomesize,
    })
}

fn pick_opt<T>(cli: Option<T>, preset: Option<T>) -> Option<T> {
    match (cli, preset) {
        (Some(c), _) => Some(c),
        (None, Some(p)) => Some(p),
        _ => None,
    }
}
fn pick_val<T: Copy>(cli: Option<T>, preset: Option<T>, default_: T) -> T {
    cli.or(preset).unwrap_or(default_)
}

// ---------- 6) Commented YAML rendering (for show/edit scaffolding) ----------
pub fn render_commented_yaml(preset: &Preset) -> String {
    // Compose a human-documented YAML (comments + current values).
    // We emit assignments only for Some(..) so it remains concise.
    // Comments also show defaults for context.
    let mut s = String::new();
    s.push_str("# POLAP preset (YAML)\n");
    s.push_str("# Lines starting with '#' are comments.\n");
    s.push_str("# Values not present here fall back to program defaults (shown below).\n\n");
    s.push_str(&format!("# Defaults: threads={}, coverage={}, single_min={}, pair_min={}, bridge_min={}, min_read_length={}, polap_reads={}, coverage_check={}, genomesize=\"{}\"\n\n",
        Defaults::THREADS, Defaults::COVERAGE, Defaults::SINGLE_MIN, Defaults::PAIR_MIN, Defaults::BRIDGE_MIN,
        Defaults::MIN_READ_LENGTH, Defaults::POLAP_READS, Defaults::COVERAGE_CHECK, Defaults::GENOMESIZE
    ));

    s.push_str("# Paths\n");
    if let Some(v) = &preset.outdir {
        s.push_str(&format!("outdir: {}\n", quote_path(v)));
    } else {
        s.push_str("# outdir: o\n");
    }
    if let Some(v) = &preset.long_reads {
        s.push_str(&format!("long_reads: {}\n", quote_path(v)));
    } else {
        s.push_str("# long_reads: <path>\n");
    }
    if let Some(v) = &preset.sr1 {
        s.push_str(&format!("sr1: {}\n", quote_path(v)));
    } else {
        s.push_str("# sr1: <path>\n");
    }
    if let Some(v) = &preset.sr2 {
        s.push_str(&format!("sr2: {}\n", quote_path(v)));
    } else {
        s.push_str("# sr2: <path>\n");
    }
    s.push('\n');

    s.push_str("# Numeric knobs\n");
    if let Some(v) = preset.threads {
        s.push_str(&format!("threads: {v}\n"));
    } else {
        s.push_str(&format!("# threads: {}\n", Defaults::THREADS));
    }
    if let Some(v) = preset.coverage {
        s.push_str(&format!("coverage: {v}\n"));
    } else {
        s.push_str(&format!("# coverage: {}\n", Defaults::COVERAGE));
    }
    if let Some(v) = preset.single_min {
        s.push_str(&format!("single_min: {v}\n"));
    } else {
        s.push_str(&format!("# single_min: {}\n", Defaults::SINGLE_MIN));
    }
    if let Some(v) = preset.pair_min {
        s.push_str(&format!("pair_min: {v}\n"));
    } else {
        s.push_str(&format!("# pair_min: {}\n", Defaults::PAIR_MIN));
    }
    if let Some(v) = preset.bridge_min {
        s.push_str(&format!("bridge_min: {v}\n"));
    } else {
        s.push_str(&format!("# bridge_min: {}\n", Defaults::BRIDGE_MIN));
    }
    if let Some(v) = preset.min_read_length {
        s.push_str(&format!("min_read_length: {v}\n"));
    } else {
        s.push_str(&format!(
            "# min_read_length: {}\n",
            Defaults::MIN_READ_LENGTH
        ));
    }
    s.push('\n');

    s.push_str("# Booleans\n");
    if let Some(v) = preset.polap_reads {
        s.push_str(&format!("polap_reads: {}\n", v));
    } else {
        s.push_str(&format!("# polap_reads: {}\n", Defaults::POLAP_READS));
    }
    if let Some(v) = preset.coverage_check {
        s.push_str(&format!("coverage_check: {}\n", v));
    } else {
        s.push_str(&format!("# coverage_check: {}\n", Defaults::COVERAGE_CHECK));
    }
    s.push('\n');

    s.push_str("# Strings\n");
    if let Some(v) = &preset.genomesize {
        s.push_str(&format!("genomesize: {:?}\n", v));
    } else {
        s.push_str(&format!("# genomesize: {:?}\n", Defaults::GENOMESIZE));
    }
    s
}

fn quote_path(p: &PathBuf) -> String {
    let s = p.to_string_lossy();
    if s.contains([' ', ':', '#']) {
        format!("{:?}", s)
    } else {
        s.into_owned()
    }
}

// ---------- 7) Preset directory utilities ----------
pub fn list_presets(dir: &Path) -> Result<Vec<PathBuf>> {
    let mut v = Vec::new();
    if dir.exists() {
        for e in fs::read_dir(dir).with_context(|| format!("read_dir {:?}", dir))? {
            let e = e?;
            if e.file_type()?.is_file() {
                let p = e.path();
                if p.extension().map(|x| x == "yaml").unwrap_or(false) {
                    v.push(p);
                }
            }
        }
        v.sort();
    }
    Ok(v)
}

pub fn show_preset(preset_path: &Path) -> Result<String> {
    let preset = load_preset(preset_path)?;
    Ok(render_commented_yaml(&preset))
}

pub fn write_commented_if_absent(preset_path: &Path) -> Result<()> {
    if preset_path.exists() {
        return Ok(());
    }
    if let Some(parent) = preset_path.parent() {
        fs::create_dir_all(parent)?;
    }
    let txt = render_commented_yaml(&Preset::default());
    let mut f = fs::File::create(preset_path)?;
    f.write_all(txt.as_bytes())?;
    Ok(())
}

pub fn edit_preset_in_editor(preset_path: &Path) -> Result<()> {
    write_commented_if_absent(preset_path)?;
    let editor = std::env::var("EDITOR")
        .ok()
        .or_else(|| std::env::var("VISUAL").ok())
        .unwrap_or_else(|| "vi".to_string());
    let status = std::process::Command::new(editor)
        .arg(preset_path)
        .status()
        .with_context(|| "spawn editor")?;
    if !status.success() {
        return Err(anyhow!("editor exited with non-zero status"));
    }
    Ok(())
}

// Add to crates/polap/src/preset.rs
pub fn to_yaml(opts: &FinalOpts) -> anyhow::Result<String> {
    // Persist as simple key–value YAML (materialize all fields)
    let p = Preset {
        outdir: Some(opts.outdir.clone()),
        long_reads: opts.long_reads.clone(),
        sr1: opts.sr1.clone(),
        sr2: opts.sr2.clone(),
        threads: Some(opts.threads),
        coverage: Some(opts.coverage),
        single_min: Some(opts.single_min),
        pair_min: Some(opts.pair_min),
        bridge_min: Some(opts.bridge_min),
        min_read_length: Some(opts.min_read_length),
        polap_reads: Some(opts.polap_reads),
        coverage_check: Some(opts.coverage_check),
        genomesize: opts.genomesize.clone(),
    };
    Ok(serde_yaml::to_string(&p)?)
}
