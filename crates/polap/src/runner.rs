// crates/polap/src/runner.rs
// Version: v0.4.1 — command runner + conda wrappers.

use anyhow::{Context, Result};
use std::process::Command;
use tracing::{error, info, warn};

use crate::cli::Cli;

pub struct CommandRunner<'a> {
    cli: &'a Cli,
}

impl<'a> CommandRunner<'a> {
    pub fn new(cli: &'a Cli) -> Self {
        Self { cli }
    }

    pub fn run(&self, program: &str, args: &[&str]) -> Result<()> {
        if self.cli.args.dry {
            info!("[DRY] {program} {}", shlex_join(args));
            return Ok(());
        }

        info!("[RUN] {program} {}", shlex_join(args));

        let mut child = Command::new(program)
            .args(args)
            .envs(self.env_overrides())
            .stdout(std::process::Stdio::piped())
            .stderr(std::process::Stdio::piped())
            .spawn()
            .with_context(|| format!("spawn {program}"))?;

        {
            let stdout = child.stdout.take().unwrap();
            let stderr = child.stderr.take().unwrap();

            let t_out = std::thread::spawn(move || {
                use std::io::{BufRead, BufReader};
                for line in BufReader::new(stdout).lines().flatten() {
                    tracing::info!("{}", line); // file log
                }
            });
            let t_err = std::thread::spawn(move || {
                use std::io::{BufRead, BufReader};
                for line in BufReader::new(stderr).lines().flatten() {
                    eprintln!("{}", line); // screen
                    tracing::warn!("{}", line); // file
                }
            });

            let _ = t_out.join();
            let _ = t_err.join();
        }

        let status = child.wait()?;
        if !status.success() {
            error!("Command failed: {program} (code: {:?})", status.code());
            anyhow::bail!("command failed");
        }
        Ok(())
    }

    /// Run a program via `conda run -n <env> ...`. Falls back if conda missing.
    pub fn run_conda(&self, env: &str, program: &str, args: &[&str]) -> Result<()> {
        if which::which("conda").is_err() {
            warn!("conda not found; falling back to `{}` directly", program);
            return self.run(program, args);
        }
        let mut argv: Vec<&str> = vec!["run", "-n", env, program];
        argv.extend_from_slice(args);
        self.run("conda", &argv)
    }

    fn env_overrides(&self) -> Vec<(String, String)> {
        vec![
            (
                "_POLAP_DEBUG".into(),
                std::env::var("_POLAP_DEBUG").unwrap_or_else(|_| "0".into()),
            ),
            (
                "_POLAP_RELEASE".into(),
                std::env::var("_POLAP_RELEASE").unwrap_or_else(|_| "0".into()),
            ),
            (
                "_POLAPLIB_DIR".into(),
                std::env::var("_POLAPLIB_DIR").unwrap_or_else(|_| ".".into()),
            ),
        ]
    }
}

fn shlex_join(args: &[&str]) -> String {
    args.iter()
        .map(|s| shell_words::quote(s))
        .collect::<Vec<_>>()
        .join(" ")
}

pub fn print_clock_to_screen() {
    // Use UTC ISO-like format for simplicity
    let now = chrono::Utc::now();
    eprintln!("{}", now.format("%Y-%m-%dT%H:%M:%SZ"));
}

pub fn print_elapsed(cli: &crate::cli::Cli) {
    let elapsed = cli.start.elapsed();
    let host = std::env::var("HOSTNAME")
        .or_else(|_| std::env::var("COMPUTERNAME"))
        .unwrap_or_else(|_| "unknown-host".into());
    eprintln!(
        "Time at {}: {}hrs {}min {}sec",
        host,
        elapsed.as_secs() / 3600,
        (elapsed.as_secs() / 60) % 60,
        elapsed.as_secs() % 60
    );
}
