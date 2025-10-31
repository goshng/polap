// Version: v0.1.0
use std::path::Path;

pub enum SeqType {
    Illumina,
    PacbioHifi,
    NanoRaw,
    Unknown,
}

pub fn detect_type(path: &Path) -> SeqType {
    let s = path.to_string_lossy().to_ascii_lowercase();
    if s.ends_with(".fq")
        || s.ends_with(".fastq")
        || s.ends_with(".fq.gz")
        || s.ends_with(".fastq.gz")
    {
        // filename hints
        if s.contains("hifi") {
            return SeqType::PacbioHifi;
        }
        if s.contains("ont") || s.contains("nano") {
            return SeqType::NanoRaw;
        }
        // default fallback to illumina for s1/s2 typical names
        if s.contains("s1") || s.contains("s2") || s.contains("_r1") || s.contains("_r2") {
            return SeqType::Illumina;
        }
    }
    SeqType::Unknown
}
