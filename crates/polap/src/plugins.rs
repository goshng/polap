// Version: v0.1.0
use crate::cli::Cli;
use crate::runner::CommandRunner;
use anyhow::{Context, Result};
use std::path::PathBuf;
use walkdir::WalkDir;

pub fn run_plugin(cli: &Cli, name: &str, args: &[String]) -> Result<()> {
    let base = std::env::var("_POLAPLIB_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|_| {
            std::env::current_exe()
                .unwrap()
                .parent()
                .unwrap()
                .join("polaplib")
        });

    let script = find_script(&base, name)
        .with_context(|| format!("plugin '{name}' not found in {}", base.display()))?;

    let argv: Vec<String> = std::iter::once(script.to_string_lossy().to_string())
        .chain(args.iter().cloned())
        .collect();

    let argv_refs: Vec<&str> = argv.iter().map(|s| s.as_str()).collect();
    let runner = CommandRunner::new(cli);

    // Decide interpreter
    let program_args = if script.ends_with(".py") {
        ("python3", &argv_refs[..])
    } else if script.ends_with(".R") || script.ends_with(".r") {
        ("Rscript", &argv_refs[..])
    } else {
        (argv_refs[0], &argv_refs[1..])
    };

    runner.run(program_args.0, program_args.1)
}

fn find_script(base: &PathBuf, name: &str) -> Option<PathBuf> {
    for entry in WalkDir::new(base).min_depth(1).max_depth(2) {
        let e = entry.ok()?;
        if !e.file_type().is_file() {
            continue;
        }
        let p = e.path();
        let fname = p.file_stem()?.to_string_lossy();
        if fname == name || p.file_name()?.to_string_lossy() == name {
            return Some(p.to_path_buf());
        }
    }
    None
}
