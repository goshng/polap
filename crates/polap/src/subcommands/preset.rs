// crates/polap/src/subcommands/preset.rs
// Version: v0.3.4 — adds `preset save <name> --print`.
use anyhow::Result;
use std::path::{Path, PathBuf};

use crate::cli::Cli;
use crate::preset::{edit_preset_in_editor, list_presets, show_preset};

pub fn run_list(cli: &Cli) -> Result<()> {
    ensure_preset_dir(&cli.args.preset_dir)?;
    let dir = &cli.args.preset_dir;
    let items = list_presets(dir)?;
    if items.is_empty() {
        println!("(no presets)  dir={}", dir.display());
        return Ok(());
    }
    for p in items {
        println!("{}", p.display());
    }
    Ok(())
}

pub fn run_show(cli: &Cli, name_or_path: &str) -> Result<()> {
    ensure_preset_dir(&cli.args.preset_dir)?;
    let path = resolve_name_or_path(name_or_path, &cli.args.preset_dir);
    let txt = show_preset(&path)?;
    println!("{txt}");
    Ok(())
}

pub fn run_edit(cli: &Cli, name_or_path: &str) -> Result<()> {
    ensure_preset_dir(&cli.args.preset_dir)?;
    let path = resolve_name_or_path(name_or_path, &cli.args.preset_dir);
    edit_preset_in_editor(&path)
}

pub fn run_save(cli: &Cli, name_or_path: &str, print: bool) -> Result<()> {
    ensure_preset_dir(&cli.args.preset_dir)?;
    let out = resolve_name_or_path(name_or_path, &cli.args.preset_dir);

    // MERGE: Defaults < existing preset at `out` < CLI
    let opts = crate::preset::finalize_options(&cli.args, Some(&out))?;

    if print {
        let yaml = crate::preset::to_yaml(&opts)?;
        println!("{yaml}");
        return Ok(());
    }

    let existed = out.exists();
    crate::preset::save_preset(&out, &opts)?;
    println!(
        "Preset {} {} at {}",
        out.file_stem().and_then(|s| s.to_str()).unwrap_or("preset"),
        if existed { "updated" } else { "created" },
        out.display()
    );
    Ok(())
}

/// Create preset dir if missing; inform user on first creation.
fn ensure_preset_dir(dir: &Path) -> Result<()> {
    if !dir.exists() {
        std::fs::create_dir_all(dir)?;
        println!(
            "Created preset directory: {}\nYou can now create presets with `polap preset edit <name>`.",
            dir.display()
        );
    }
    Ok(())
}

/// Resolve a preset name or path to a concrete YAML file.
/// - If `s` is an existing path → use it directly.
/// - If `s` ends with `.yaml` (even if not yet existing) → treat as full path.
/// - Otherwise → create under preset dir with `.yaml` extension.
fn resolve_name_or_path(s: &str, dir: &Path) -> PathBuf {
    let p = Path::new(s);
    // 1. Already an existing file or absolute/relative path with slash
    if p.exists() {
        return p.to_path_buf();
    }
    // 2. If user explicitly typed a .yaml filename, treat it as a full path
    if s.ends_with(".yaml") {
        return PathBuf::from(s);
    }
    // 3. Otherwise, join with preset dir and append .yaml
    dir.join(format!("{s}.yaml"))
}
