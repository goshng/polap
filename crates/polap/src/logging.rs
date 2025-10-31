// crates/polap/src/logging.rs
// Version: v0.4.5 — event_format + ChronoUtc timer (correct type params)

use anyhow::Result;
use tracing_appender::rolling;
use tracing_subscriber::prelude::*;
use tracing_subscriber::{fmt, EnvFilter};

use crate::cli::{Cli, LogTs};

fn level_from(cli: &Cli) -> &'static str {
    if cli.args.quiet {
        "error"
    } else {
        match cli.args.verbose {
            0 => "info",
            1 => "debug",
            _ => "trace", // -vv / -vvv => trace to screen
        }
    }
}

// Build an event formatter that uses the `Full` formatter kind (the default),
// and a ChronoUtc timer (so `.event_format(...)` satisfies bounds).
fn make_event_format(ts: LogTs) -> fmt::format::Format<fmt::format::Full, fmt::time::ChronoUtc> {
    use fmt::time::ChronoUtc;
    let timer = match ts {
        LogTs::Off => ChronoUtc::new("".into()), // emits no timestamp
        LogTs::Date => ChronoUtc::new("%Y-%m-%d".into()),
        LogTs::Time => ChronoUtc::new("%H:%M:%S".into()),
        LogTs::Datetime => ChronoUtc::new("%Y-%m-%d %H:%M:%S".into()),
    };
    fmt::format().with_timer(timer) // NOTE: returns Format<Full, ChronoUtc>
}

pub fn init(cli: &Cli) -> Result<tracing_appender::non_blocking::WorkerGuard> {
    // File appender
    let file_app = rolling::never(
        cli.out_log
            .parent()
            .unwrap_or_else(|| std::path::Path::new(".")),
        cli.out_log.file_name().unwrap(),
    );
    let (file_writer, guard) = tracing_appender::non_blocking(file_app);

    // Event formatters (screen/file) — same timer config on both
    let ev_file = make_event_format(cli.args.log_ts);
    let ev_term = make_event_format(cli.args.log_ts);

    // File layer: full detail, with file:line
    let file_layer = fmt::layer()
        .event_format(ev_file)
        .with_writer(file_writer)
        .with_ansi(false)
        .with_target(false)
        .with_level(false)
        .with_file(true)
        .with_line_number(true)
        .with_thread_names(false);

    // Screen layer: respect -q/-v; file:line only in debug builds
    let term_layer = fmt::layer()
        .event_format(ev_term)
        .with_writer(std::io::stderr)
        .with_ansi(atty::is(atty::Stream::Stderr))
        .with_target(false)
        .with_level(false)
        .with_file(cfg!(debug_assertions))
        .with_line_number(cfg!(debug_assertions))
        .with_thread_names(false);

    // Screen filter from -q/-v
    let screen_filter = EnvFilter::builder().parse(level_from(cli))?;

    // Registry: file logs everything; screen logs per level_from
    let base = tracing_subscriber::registry()
        .with(EnvFilter::builder().parse("trace")?) // file layer sees all
        .with(file_layer);

    if cli.args.quiet {
        base.init();
    } else {
        base.with(screen_filter).with(term_layer).init();
    }

    Ok(guard)
}
