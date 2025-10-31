// Version: v0.1.0
use anyhow::{bail, Result};

pub fn parse_human(s: &str) -> Result<u64> {
    // "5k", "10M", "1.5G" → bytes/bases
    let s = s.trim();
    let (num, suf) = s.split_at(s.len() - 1);
    let factor = match suf.to_ascii_lowercase().as_str() {
        "k" => 1_000u64,
        "m" => 1_000_000,
        "g" => 1_000_000_000,
        _ => {
            // maybe pure number
            return s.parse::<u64>().map_err(|e| e.into());
        }
    };
    let val: f64 = num
        .parse()
        .map_err(|_| anyhow::anyhow!("bad number: {num}"))?;
    Ok((val * factor as f64) as u64)
}

// jellyfish -s sanitizer: "5g|5m|5k" → G/M/k normalized or "0" if invalid
pub fn jellyfish_s_normalize(s: &str) -> String {
    let s = s.trim();
    if s.len() < 2 {
        return "0".into();
    }
    let (num, suf) = s.split_at(s.len() - 1);
    if !num.chars().all(|c| c.is_ascii_digit()) {
        return "0".into();
    }
    match suf {
        "g" | "G" => format!("{num}G"),
        "m" | "M" => format!("{num}M"),
        "k" | "K" => format!("{num}k"),
        _ => "0".into(),
    }
}
