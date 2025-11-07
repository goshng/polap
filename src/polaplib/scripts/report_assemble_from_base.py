#!/usr/bin/env python3
# report_assemble_from_base.py
# Version: v0.3.0
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Build a structured JSON report by scanning a polap-assemble BASE DIR:
#   <Species>/v6/<inum>/polap-assemble
#
# Sections detected (if files exist):
#   1. Seed PT and MT reads
#   2. ptDNA assembly
#   3. Filter out PT reads
#   4. Compute read overlapness
#   5. Filter out NT reads
#   6. Assemble seed contigs using miniasm or mtDNA assembly mt0
#   7. mtDNA assembly  (auto-subsections: 7.1 mt1, 7.2 mt2, 7.3 mt3)
#   8. Oatk's pathfinder
#   9. Polishing
#
# Each entry records:
#   - datetime/date/time from file mtime (UTC, with ns)
#   - size (bytes)
#   - path (displayed as '<Species>/...')
#   - comment (context)
#   - ext
#
import os
import re
import json
import argparse
import datetime
from typing import Dict, Any, List, Optional, Tuple
from glob import glob

# ---------------------------- helpers ----------------------------


def stat_datetime_ns(p: str) -> Tuple[str, str, str, int]:
    """
    Return (datetime, date, time_frac, size) from mtime ns (UTC).
    time_frac like HH:MM:SS.nnnnnnnnn
    """
    st = os.stat(p)
    ns = getattr(st, "st_mtime_ns", int(st.st_mtime * 1_000_000_000))
    sec = ns // 1_000_000_000
    nsec = ns % 1_000_000_000
    dt_utc = datetime.datetime.utcfromtimestamp(sec)
    date = dt_utc.strftime("%Y-%m-%d")
    time_frac = dt_utc.strftime("%H:%M:%S") + f".{nsec:09d}"
    dt_str = f"{date} {time_frac}"
    return dt_str, date, time_frac, int(st.st_size)


def detect_species_from_base(
    base_dir: str,
) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """
    Walk upward from base_dir and return (species_name, species_root_abs, species_parent_abs)
    when we find a 'v<digits>' directory; species is its parent basename.
    """
    cur = os.path.abspath(base_dir)
    up_seen = set()
    while True:
        base = os.path.basename(cur)
        if re.fullmatch(r"v\d+", base):
            species_root = os.path.dirname(cur)  # .../<Species>
            species = os.path.basename(species_root)
            species_parent = os.path.dirname(species_root)
            return species, species_root, species_parent
        up = os.path.dirname(cur)
        if up == cur or up in up_seen:
            break
        up_seen.add(up)
        cur = up
    return None, None, None


def rel_species_display(
    abs_path: str, species_parent: Optional[str], species: Optional[str]
) -> str:
    """
    Turn an absolute path into '<Species>/...' (display path).
    If species_parent known, compute relpath from it.
    Otherwise, fallback to basename surgery.
    """
    if species_parent and species:
        try:
            rel = os.path.relpath(abs_path, start=species_parent)
            # Ensure it starts with <Species>/
            if not rel.startswith(species + os.sep):
                # Fallback: find '/<Species>/' segment inside.
                i = abs_path.find(os.sep + species + os.sep)
                if i != -1:
                    return abs_path[i + 1 :]  # drop leading '/'
            return rel
        except Exception:
            pass
    # Fallback: search for '<Species>/' fragment
    if species:
        i = abs_path.find(os.sep + species + os.sep)
        if i != -1:
            return abs_path[i + 1 :]
    return os.path.basename(abs_path)


def add_file_entry(
    entries: List[Dict[str, Any]],
    abs_path: str,
    species_parent: Optional[str],
    species: Optional[str],
    comment: str,
) -> None:
    if not abs_path or not os.path.exists(abs_path):
        return
    dt_str, date, time_frac, size = stat_datetime_ns(abs_path)
    display = rel_species_display(abs_path, species_parent, species)
    _, ext = os.path.splitext(abs_path)
    entries.append(
        {
            "datetime": dt_str,
            "date": date,
            "time": time_frac,
            "size": size,
            "path": display.replace("\\", "/"),
            "comment": comment,
            "ext": ext.lower(),
        }
    )


def add_dir_entry(
    entries: List[Dict[str, Any]],
    abs_dir: str,
    species_parent: Optional[str],
    species: Optional[str],
    comment: str,
) -> None:
    if not abs_dir or not os.path.isdir(abs_dir):
        return
    dt_str, date, time_frac, size = stat_datetime_ns(abs_dir)
    display = rel_species_display(abs_dir, species_parent, species)
    entries.append(
        {
            "datetime": dt_str,
            "date": date,
            "time": time_frac,
            "size": size,
            "path": display.replace("\\", "/"),
            "comment": comment,
            "ext": "/",
        }
    )


def aggregate_stats(entries: List[Dict[str, Any]]) -> Dict[str, Any]:
    entries_sorted = sorted(entries, key=lambda e: (e["date"], e["time"]))
    total = sum(int(e.get("size", 0) or 0) for e in entries_sorted)
    t0 = entries_sorted[0]["datetime"] if entries_sorted else None
    t1 = entries_sorted[-1]["datetime"] if entries_sorted else None
    return {
        "count": len(entries_sorted),
        "total_bytes": total,
        "time_start": t0,
        "time_end": t1,
    }, entries_sorted


def mk_section(
    sec_id: str,
    title: str,
    entries: List[Dict[str, Any]],
    subsections: Optional[List[Dict[str, Any]]] = None,
) -> Dict[str, Any]:
    stats, ents_sorted = aggregate_stats(entries)
    node: Dict[str, Any] = {
        "id": sec_id,
        "title": title,
        "entries": ents_sorted,
        "stats": stats,
    }
    if subsections:
        node["subsections"] = subsections
    return node


# ---------------------------- collectors per section ----------------------------


def collect_section1(
    base_dir: str, species_parent: Optional[str], species: Optional[str]
) -> Dict[str, Any]:
    """1. Seed PT and MT reads"""
    ents: List[Dict[str, Any]] = []
    base = os.path.join(base_dir, "annotate-read-pt")
    add_file_entry(
        ents,
        os.path.join(base, "assembly_info_organelle_annotation_count-all.txt"),
        species_parent,
        species,
        "Annotation MT PT for all",
    )
    add_file_entry(
        ents,
        os.path.join(base, "contig-annotation-depth-table.txt"),
        species_parent,
        species,
        "Annotation MT",
    )
    add_file_entry(
        ents,
        os.path.join(base, "pt-contig-annotation-depth-table.txt"),
        species_parent,
        species,
        "Annotation PT",
    )
    add_file_entry(
        ents, os.path.join(base, "pt.id.txt"), species_parent, species, "PT reads"
    )
    # common variants for seqkit stats
    for cand in (
        "pt0.fq.seqkit.stats.ta.txt",
        "pt.fq.seqkit.stats.ta.txt",
        "pt0.fq.seqkit.stats.T.txt",
    ):
        add_file_entry(
            ents,
            os.path.join(base, cand),
            species_parent,
            species,
            "PT reads seqkit stats",
        )
    return mk_section("1", "Seed PT and MT reads", ents)


def collect_section2(
    base_dir: str, species_parent: Optional[str], species: Optional[str]
) -> Dict[str, Any]:
    """2. ptDNA assembly"""
    ents: List[Dict[str, Any]] = []
    pt0 = os.path.join(base_dir, "annotate-read-pt", "pt")
    add_file_entry(
        ents,
        os.path.join(pt0, "30-contigger", "graph_final.gfa"),
        species_parent,
        species,
        "Flye ptDNA assembly pt0 from PT reads alone",
    )
    add_file_entry(
        ents,
        os.path.join(pt0, "30-contigger", "graph_final.png"),
        species_parent,
        species,
        "pt0 assembly graph png",
    )
    add_file_entry(
        ents,
        os.path.join(pt0, "contig-annotation-depth-table.txt"),
        species_parent,
        species,
        "Annotation MT from assembly pt0",
    )
    add_file_entry(
        ents,
        os.path.join(pt0, "pt-contig-annotation-depth-table.txt"),
        species_parent,
        species,
        "Annotation PT from assembly pt0",
    )
    add_file_entry(
        ents,
        os.path.join(pt0, "assembly_info_organelle_annotation_count-all.txt"),
        species_parent,
        species,
        "Assembly info (pt0)",
    )
    add_file_entry(
        ents,
        os.path.join(pt0, "mt.contig.name-pt1"),
        species_parent,
        species,
        "Seed contigs in pt0 assembly for pt1",
    )

    pt1 = os.path.join(base_dir, "annotate-read-pt", "pt1")
    add_file_entry(
        ents,
        os.path.join(pt1, "01-contig", "contig.fa"),
        species_parent,
        species,
        "Seed contig sequence from pt0 for pt1",
    )
    add_file_entry(
        ents,
        os.path.join(pt1, "assembly_graph.gfa"),
        species_parent,
        species,
        "Flye assembly graph of pt1",
    )
    add_file_entry(
        ents,
        os.path.join(pt1, "assembly.fasta"),
        species_parent,
        species,
        "Flye assembly pt1 sequence",
    )
    add_file_entry(
        ents,
        os.path.join(pt1, "ptdna.metrics.csv"),
        species_parent,
        species,
        "ptDNA reads info in pt1",
    )
    add_file_entry(
        ents,
        os.path.join(pt1, "assembly_graph.png"),
        species_parent,
        species,
        "Assembly pt1 graph png",
    )
    add_file_entry(
        ents,
        os.path.join(pt1, "contig-annotation-depth-table.txt"),
        species_parent,
        species,
        "Assembly pt1 annotation table for MT",
    )
    add_file_entry(
        ents,
        os.path.join(pt1, "pt-contig-annotation-depth-table.txt"),
        species_parent,
        species,
        "Assembly pt1 annotation table for PT",
    )

    # Scatter plots from annotate-read-mtseed (useful context)
    arm = os.path.join(base_dir, "annotate-read-mtseed")
    add_file_entry(
        ents,
        os.path.join(arm, "pt", "pt-contig-annotation-depth-table.txt.scatter.pdf"),
        species_parent,
        species,
        "PT read scatter plot",
    )
    add_file_entry(
        ents, os.path.join(arm, "pt.id.all.txt"), species_parent, species, "PT reads"
    )
    add_file_entry(
        ents,
        os.path.join(arm, "mt", "contig-annotation-depth-table.txt.scatter.pdf"),
        species_parent,
        species,
        "MT read scatter plot",
    )
    add_file_entry(
        ents, os.path.join(arm, "mt.id.all.txt"), species_parent, species, "MT reads"
    )

    return mk_section("2", "ptDNA assembly", ents)


def collect_section3(
    base_dir: str, species_parent: Optional[str], species: Optional[str]
) -> Dict[str, Any]:
    """3. Filter out PT reads"""
    ents: List[Dict[str, Any]] = []
    mtseed = os.path.join(base_dir, "mtseed")
    panel = os.path.join(mtseed, "01-panel")
    add_file_entry(
        ents,
        os.path.join(panel, "pt_isomerA.fa"),
        species_parent,
        species,
        "ptDNA isomer A sequence fasta",
    )
    add_file_entry(
        ents,
        os.path.join(panel, "pt_isomerB.fa"),
        species_parent,
        species,
        "ptDNA isomer B sequence fasta",
    )
    add_file_entry(
        ents,
        os.path.join(mtseed, "pt_thresh.diag.tsv"),
        species_parent,
        species,
        "ptDNA assembly pt0 PT reads selection criteria",
    )
    add_file_entry(
        ents,
        os.path.join(mtseed, "pt.ids"),
        species_parent,
        species,
        "PT reads IDs selected for ptDNA assembly pt0",
    )
    add_file_entry(
        ents,
        os.path.join(mtseed, "pt_thresh.vars"),
        species_parent,
        species,
        "pt0 reads selection thresholds from pt_thresh.diag.tsv",
    )

    # stats variants
    stat_candidates = [
        "reads.nonpt.fq.gz.seqkit.stats.ta.tsv",
        "reads.nonpt.fq.gz.seqkit.stats.ta.txt",
        "reads.nonpt.fq.gz.seqkit.stats.T.txt",
    ]
    for cand in stat_candidates:
        add_file_entry(
            ents,
            os.path.join(mtseed, cand),
            species_parent,
            species,
            "read seqkit stats after filtering out PT reads",
        )
    return mk_section("3", "Filter out PT reads", ents)


def collect_section4(
    base_dir: str, species_parent: Optional[str], species: Optional[str]
) -> Dict[str, Any]:
    """4. Compute read overlapness"""
    ents: List[Dict[str, Any]] = []
    m3 = os.path.join(base_dir, "mtseed", "03-allvsall")
    add_file_entry(
        ents,
        os.path.join(m3, "edges_loose.tsv.35.txt"),
        species_parent,
        species,
        "Edges above 35 cutoff (summary)",
    )
    qc = os.path.join(m3, "04-qc")
    add_file_entry(
        ents,
        os.path.join(qc, "scan_degree_hist.pdf"),
        species_parent,
        species,
        "Distribution of degree values of all reads",
    )
    add_file_entry(
        ents,
        os.path.join(qc, "scan_wdegree_hist.pdf"),
        species_parent,
        species,
        "Distribution of weighted degree values of all reads",
    )
    add_file_entry(
        ents,
        os.path.join(qc, "scan_cum_wdegree.pdf"),
        species_parent,
        species,
        "Cumulative distribution of weighted degree values (threshold)",
    )
    add_file_entry(
        ents,
        os.path.join(qc, "overlap_qc.vars"),
        species_parent,
        species,
        "Overlapness QC values",
    )
    add_file_entry(
        ents,
        os.path.join(base_dir, "mtseed", "pipeline.log"),
        species_parent,
        species,
        "mtseed procedure log file",
    )
    return mk_section("4", "Compute read overlapness", ents)


def collect_section5(
    base_dir: str, species_parent: Optional[str], species: Optional[str]
) -> Dict[str, Any]:
    """5. Filter out NT reads"""
    ents: List[Dict[str, Any]] = []
    add_file_entry(
        ents,
        os.path.join(base_dir, "mtseed", "04-busco", "nuc.ids.sample"),
        species_parent,
        species,
        "NT reads IDs",
    )
    add_file_entry(
        ents,
        os.path.join(base_dir, "mtseed", "05-round", "threshold_from_nuclear.tsv"),
        species_parent,
        species,
        "Weighted-degree threshold to filter NT reads",
    )
    return mk_section("5", "Filter out NT reads", ents)


def collect_section6(
    base_dir: str, species_parent: Optional[str], species: Optional[str]
) -> Dict[str, Any]:
    """6. Assemble seed contigs using miniasm or mtDNA assembly mt0"""
    ents: List[Dict[str, Any]] = []
    m6 = os.path.join(base_dir, "mtseed", "06-miniasm")
    add_file_entry(
        ents,
        os.path.join(m6, "log", "steps.log"),
        species_parent,
        species,
        "Logs from miniasm executions",
    )
    m7 = os.path.join(base_dir, "mtseed", "07-flye")
    add_dir_entry(ents, m7, species_parent, species, "Folder for assembly mt0")
    add_file_entry(
        ents,
        os.path.join(m7, "30-contigger", "graph_final.png"),
        species_parent,
        species,
        "miniasm assembly graph png",
    )
    # '3-gfa.all.gfa' may live here or be named slightly differently
    add_file_entry(
        ents,
        os.path.join(m7, "30-contigger", "3-gfa.all.gfa"),
        species_parent,
        species,
        "miniasm assembly graph (GFA without sequences)",
    )
    add_file_entry(
        ents,
        os.path.join(m7, "contig-annotation-depth-table.txt"),
        species_parent,
        species,
        "MT annotation table for miniasm assembly (mt0)",
    )
    add_file_entry(
        ents,
        os.path.join(m7, "pt-contig-annotation-depth-table.txt"),
        species_parent,
        species,
        "PT annotation table for miniasm assembly (mt0)",
    )
    add_file_entry(
        ents,
        os.path.join(m7, "mt.contig.name-mt1"),
        species_parent,
        species,
        "Seed contigs selected from mt0 for mt1",
    )
    return mk_section(
        "6", "Assemble seed contigs using miniasm or mtDNA assembly mt0", ents
    )


def collect_section7(
    base_dir: str, species_parent: Optional[str], species: Optional[str]
) -> Dict[str, Any]:
    """7. mtDNA assembly with subsections"""
    subsections: List[Dict[str, Any]] = []

    mtseed = os.path.join(base_dir, "mtseed")
    for name, idx in (("mt1", "7.1"), ("mt2", "7.2"), ("mt3", "7.3")):
        d = os.path.join(mtseed, name)
        if not os.path.isdir(d):
            continue
        ents: List[Dict[str, Any]] = []
        add_dir_entry(ents, d, species_parent, species, f"Folder for assembly {name}")
        # reads runs
        for length_txt in glob(
            os.path.join(d, "02-reads", "ptgaul-reads", "*", "ptgaul.names.length.txt")
        ):
            run = os.path.basename(os.path.dirname(length_txt))
            add_file_entry(
                ents,
                length_txt,
                species_parent,
                species,
                f"Reads total from the ptgaul selection run {run}",
            )
            names_txt = os.path.join(os.path.dirname(length_txt), "ptgaul.names")
            add_file_entry(
                ents,
                names_txt,
                species_parent,
                species,
                f"Read IDs from the ptgaul selection run {run}",
            )

        # selected indices/omega
        add_file_entry(
            ents,
            os.path.join(d, "01-contig", "index.txt"),
            species_parent,
            species,
            "Selected run index",
        )
        add_file_entry(
            ents,
            os.path.join(d, "01-contig", "omega.txt"),
            species_parent,
            species,
            "Selected omega to use",
        )
        add_file_entry(
            ents,
            os.path.join(d, "01-contig", "contig_total_length.txt"),
            species_parent,
            species,
            "Seed contig total length (bp)",
        )

        # Any test flye runs (paths vary)
        for flye_log in glob(
            os.path.join(d, "05-flye", "**", "flye.log"), recursive=True
        ):
            add_file_entry(ents, flye_log, species_parent, species, "Test Flye run log")
        for gfa in glob(os.path.join(d, "06-summary", "**", "*.gfa"), recursive=True):
            add_file_entry(ents, gfa, species_parent, species, "Test Flye run GFA")
        for frags in glob(
            os.path.join(d, "06-summary", "**", "*.fragments"), recursive=True
        ):
            add_file_entry(
                ents, frags, species_parent, species, "Test Flye run fragments"
            )
        for bases in glob(
            os.path.join(d, "06-summary", "**", "*.bases"), recursive=True
        ):
            add_file_entry(
                ents, bases, species_parent, species, "Test Flye run total bases"
            )
        for depth in glob(
            os.path.join(d, "06-summary", "**", "*.depth"), recursive=True
        ):
            add_file_entry(ents, depth, species_parent, species, "Test Flye run depths")

        # main assembly artifacts
        add_file_entry(
            ents,
            os.path.join(d, "assembly_graph.png"),
            species_parent,
            species,
            f"Flye assembly graph png for {name}",
        )
        add_file_entry(
            ents,
            os.path.join(d, "30-contigger", "3-gfa.all.gfa"),
            species_parent,
            species,
            f"Flye assembly {name} graph without sequences",
        )
        add_file_entry(
            ents,
            os.path.join(d, "contig-annotation-depth-table.txt"),
            species_parent,
            species,
            f"MT annotation table for assembly {name}",
        )
        add_file_entry(
            ents,
            os.path.join(d, "pt-contig-annotation-depth-table.txt"),
            species_parent,
            species,
            f"PT annotation table for assembly {name}",
        )

        # next seeds marker
        nxt = {"mt1": "mt2", "mt2": "mt3"}.get(name)
        if nxt:
            add_file_entry(
                ents,
                os.path.join(d, f"mt.contig.name-{nxt}"),
                species_parent,
                species,
                f"Seed contigs from assembly {name} for assembly {nxt}",
            )

        # push subsection
        subsections.append(mk_section(idx, f"mtDNA assembly {name}", ents))

    # parent section (entries empty; only subsections)
    parent = {
        "id": "7",
        "title": "mtDNA assembly",
        "entries": [],
        "stats": {"count": 0, "total_bytes": 0, "time_start": None, "time_end": None},
        "subsections": subsections,
    }
    return parent


def collect_section8(
    base_dir: str, species_parent: Optional[str], species: Optional[str]
) -> Dict[str, Any]:
    """8. Oatk's pathfinder"""
    ents: List[Dict[str, Any]] = []
    ex = os.path.join(base_dir, "extract")
    add_file_entry(
        ents,
        os.path.join(ex, "pt.norm.gfa"),
        species_parent,
        species,
        "Assembly graph PT for Oatk pathfinder",
    )
    add_file_entry(
        ents,
        os.path.join(ex, "mt.norm.gfa"),
        species_parent,
        species,
        "Assembly graph MT for Oatk pathfinder",
    )
    add_file_entry(
        ents,
        os.path.join(ex, "oatk.pltd.ctg.fasta"),
        species_parent,
        species,
        "PT DNA sequence by Oatk pathfinder",
    )
    add_file_entry(
        ents,
        os.path.join(ex, "oatk.mito.ctg.fasta"),
        species_parent,
        species,
        "MT DNA sequence by Oatk pathfinder",
    )
    return mk_section("8", "Oatk's pathfinder", ents)


def collect_section9(
    base_dir: str, species_parent: Optional[str], species: Optional[str]
) -> Dict[str, Any]:
    """9. Polishing"""
    ents: List[Dict[str, Any]] = []
    p = os.path.join(base_dir, "polish-longshort")
    logs = os.path.join(p, "logs")
    add_file_entry(
        ents,
        os.path.join(logs, "r1.minimap2.err"),
        species_parent,
        species,
        "Minimap2 round 1 log (time/memory last line)",
    )
    add_file_entry(
        ents,
        os.path.join(logs, "r1.filter.err"),
        species_parent,
        species,
        "Filter round 1 log (reads used)",
    )
    add_file_entry(
        ents,
        os.path.join(p, "stage1", "polished.r1.fa"),
        species_parent,
        species,
        "Racon round 1 polished sequence",
    )
    add_file_entry(
        ents,
        os.path.join(logs, "r1.racon.err"),
        species_parent,
        species,
        "Racon round 1 log (end time)",
    )
    add_file_entry(
        ents,
        os.path.join(logs, "r2.minimap2.err"),
        species_parent,
        species,
        "Minimap2 round 2 log (time/memory last line)",
    )
    add_file_entry(
        ents,
        os.path.join(logs, "r2.filter.err"),
        species_parent,
        species,
        "Filter round 2 log (reads used)",
    )
    add_file_entry(
        ents,
        os.path.join(p, "stage1", "polished.r2.fa"),
        species_parent,
        species,
        "Racon round 2 polished sequence",
    )
    add_file_entry(
        ents,
        os.path.join(logs, "r2.racon.err"),
        species_parent,
        species,
        "Racon round 2 log (end time)",
    )
    add_file_entry(
        ents,
        os.path.join(logs, "ropebwt2.err"),
        species_parent,
        species,
        "ropebwt2 log (end time)",
    )
    add_file_entry(
        ents,
        os.path.join(logs, "fmlrc2-convert.err"),
        species_parent,
        species,
        "fmlrc2 convert log",
    )
    add_file_entry(
        ents, os.path.join(logs, "fmlrc2.err"), species_parent, species, "fmlrc2 log"
    )
    add_file_entry(
        ents,
        os.path.join(p, "stage2", "segments.fmlrc2.fa"),
        species_parent,
        species,
        "fmlrc2 result sequence",
    )
    add_file_entry(
        ents,
        os.path.join(logs, "bt2-build.out"),
        species_parent,
        species,
        "bowtie2 build stdout",
    )
    add_file_entry(
        ents,
        os.path.join(logs, "bt2-build.err"),
        species_parent,
        species,
        "bowtie2 build stderr",
    )
    add_file_entry(
        ents,
        os.path.join(logs, "polypolish.bt2.err"),
        species_parent,
        species,
        "polypolish bowtie2 log",
    )
    add_file_entry(
        ents,
        os.path.join(p, "stage3", "polished.polypolish.fa"),
        species_parent,
        species,
        "polypolish result (FASTA)",
    )
    add_file_entry(
        ents,
        os.path.join(logs, "polypolish.err"),
        species_parent,
        species,
        "polypolish log",
    )
    add_file_entry(
        ents,
        os.path.join(p, "polished.fa"),
        species_parent,
        species,
        "Racon×2 + fmlrc2 + polypolish final (FASTA)",
    )
    return mk_section("9", "Polishing", ents)


# ---------------------------- main ----------------------------


def main():
    ap = argparse.ArgumentParser(
        description="Scan a polap-assemble base dir and build JSON report."
    )
    ap.add_argument(
        "--base-dir", required=True, help="Path to <Species>/vX/<inum>/polap-assemble"
    )
    ap.add_argument("--out", required=True, help="Output JSON file")
    args = ap.parse_args()

    base_dir = os.path.abspath(args.base_dir)
    if not os.path.isdir(base_dir):
        raise SystemExit(f"[ERR] --base-dir not found: {base_dir}")

    species, species_root, species_parent = detect_species_from_base(base_dir)

    sections: List[Dict[str, Any]] = []
    sections.append(collect_section1(base_dir, species_parent, species))
    sections.append(collect_section2(base_dir, species_parent, species))
    sections.append(collect_section3(base_dir, species_parent, species))
    sections.append(collect_section4(base_dir, species_parent, species))
    sections.append(collect_section5(base_dir, species_parent, species))
    sections.append(collect_section6(base_dir, species_parent, species))
    sections.append(collect_section7(base_dir, species_parent, species))
    sections.append(collect_section8(base_dir, species_parent, species))
    sections.append(collect_section9(base_dir, species_parent, species))

    # Remove empty sections (no entries and no subsections)
    filtered_sections: List[Dict[str, Any]] = []
    for s in sections:
        keep = False
        if s.get("entries"):
            keep = True
        if s.get("subsections"):
            if any(ss.get("entries") for ss in s["subsections"]):
                keep = True
        if keep:
            filtered_sections.append(s)

    doc = {
        "base_dir": base_dir,
        "species": species,
        "generated_at": datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
        "sections": filtered_sections,
    }

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as fo:
        json.dump(doc, fo, indent=2)

    print(f"[OK] Wrote JSON: {os.path.abspath(args.out)}")


if __name__ == "__main__":
    raise SystemExit(main())
