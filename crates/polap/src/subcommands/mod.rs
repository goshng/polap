// crates/polap/src/subcommands/mod.rs
use crate::cli::{Cli, Cmd, PresetAction, TestAction};
use anyhow::Result;

mod preset;
mod test;

pub fn dispatch(cli: &mut Cli) -> Result<()> {
    match &cli.args.cmd {
        None => {
            if cli.args.preset.is_some() || cli.args.preset_path.is_some() || cli.args.preset_save {
                return Ok(());
            }
            print_help_like_argv0();
            Ok(())
        }

        Some(Cmd::Preset { action }) => match action {
            PresetAction::List => preset::run_list(cli),
            PresetAction::Show { name_or_path } => preset::run_show(cli, name_or_path),
            PresetAction::Edit { name_or_path } => preset::run_edit(cli, name_or_path),
            PresetAction::Save {
                name_or_path,
                print,
            } => preset::run_save(cli, name_or_path, *print),
        },

        Some(Cmd::Test { action }) => match action {
            TestAction::Log => test::run_log(cli),
            TestAction::Bash => test::run_bash(cli),
            TestAction::R => test::run_r(cli),
            TestAction::Python => test::run_python(cli),
        },

        // TEMP stubs so the match is exhaustive
        Some(Cmd::Assemble(_)) => Ok(()),
        Some(Cmd::PreparePolishing(_)) => Ok(()),
        Some(Cmd::Polish(_)) => Ok(()),

        Some(Cmd::Plugin { name, args }) => crate::plugins::run_plugin(cli, name, args),
    }
}

fn print_help_like_argv0() {
    println!("POLAP – use `polap --help` or a subcommand like `polap preset --help` / `polap test --help`.");
}
