// Version: v0.1.0
use anyhow::Result;

mod cli;
mod fastq;
mod logging;
mod man;
mod plugins;
mod preset;
mod profiles;
mod runner;
mod subcommands;
mod units;

fn main() -> Result<()> {
    color_eyre::install().ok(); // optional if you add color-eyre
    let mut cli = cli::Cli::parse_and_init()?; // loads YAML, merges, sets outdirs
    let _guard = logging::init(&cli)?; // tee to file + stderr

    tracing::info!("POLAP {}", cli.version_string());

    // Clock flag parity with Bash
    if cli.args.clock {
        runner::print_clock_to_screen();
    }

    subcommands::dispatch(&mut cli)?;

    if cli.args.clock {
        runner::print_clock_to_screen();
    }
    runner::print_elapsed(&cli);

    Ok(())
}
