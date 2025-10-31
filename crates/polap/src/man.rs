// crates/polap/src/man.rs
// Version: v0.1.0
use anyhow::{Context, Result};
use std::{io::Write, path::Path};
use tempfile::NamedTempFile;

/// Convert a simple "section: ...\n" text into a minimal roff manpage and open it with `man -l`.
///
/// `title` becomes the `.TH` title (e.g. "polap-assemble2").
/// `section_no` is typically "1".
/// `body` is your here-doc style content (Headings like "Name:", "Synopsis:", ...).
/// If `man` is missing, we fall back to printing to stdout.
pub fn page_text_as_man(title: &str, section_no: &str, body: &str, version: &str) -> Result<()> {
    let roff = text_to_man_roff(title, section_no, body, version);
    let mut tf = NamedTempFile::new()?;
    tf.write_all(roff.as_bytes())?;
    let path = tf.into_temp_path(); // temp file that deletes on drop (we keep it alive)
                                    // Try `man -l` with a less pager; if not available, just print.
    let status = std::process::Command::new("man")
        .args(["-l", path.to_string_lossy().as_ref()])
        .status();

    match status {
        Ok(s) if s.success() => Ok(()),
        _ => {
            // Fallback: print to stdout
            println!("{}", roff);
            Ok(())
        }
    }
}

/// Very small parser: turns "Name:\n..", "Synopsis:\n..", "Description:\n..", "Options:\n..",
/// "Inputs:", "Outputs:", "Examples:", "Author:", "Copyright:"
/// into .SH sections; everything else becomes paragraph text.
///
/// NOTE: This is intentionally simple & robust. You can extend the `match` as needed.
fn text_to_man_roff(title: &str, section_no: &str, body: &str, version: &str) -> String {
    use std::borrow::Cow;

    fn escape_roff(line: &str) -> String {
        // Minimal escapes so leading dots/hyphens don't become directives.
        if line.starts_with('.') || line.starts_with('\'') {
            format!("\\&{}", line)
        } else {
            line.replace("\\", "\\\\")
        }
    }

    let date = chrono::Local::now().format("%Y-%m-%d").to_string();
    let mut out = String::new();
    // header
    out.push_str(&format!(
        ".TH \"{}\" {} \"{}\" \"{}\" \"POLAP\"\n",
        title, section_no, date, version
    ));

    let mut section_open = false;
    for raw in body.lines() {
        let line = raw.trim_end();

        // Match headings of form "Name:" etc.
        if let Some(h) = line.strip_suffix(':') {
            let sh = match h.trim().to_ascii_lowercase().as_str() {
                "name" => "NAME",
                "synopsis" => "SYNOPSIS",
                "description" => "DESCRIPTION",
                "options" => "OPTIONS",
                "inputs" => "INPUTS",
                "outputs" => "OUTPUTS",
                "examples" => "EXAMPLES",
                "author" => "AUTHOR",
                "authors" => "AUTHORS",
                "copyright" => "COPYRIGHT",
                _ => "DESCRIPTION",
            };
            out.push_str(&format!(".SH {}\n", sh));
            section_open = true;
            continue;
        }

        // blank line -> paragraph break
        if line.trim().is_empty() {
            out.push_str("\n.PP\n");
            continue;
        }

        // normal text line
        if !section_open {
            out.push_str(".SH DESCRIPTION\n");
            section_open = true;
        }
        out.push_str(&escape_roff(line));
        out.push('\n');
    }

    out
}
