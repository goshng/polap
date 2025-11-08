## `polaplib/tooltest/minimap2/README.md`

# minimap2 — smoke & extended tests

**Folder:** `polaplib/tooltest/minimap2/`

## What this test covers

- **Smoke** (default): `minimap2 --version` (tool available on `PATH`)
- **Extended** (`--extended`): tiny FASTA + FASTQ → SAM; header & alignment presence

## Prerequisites

- `minimap2` on `PATH` (e.g., `mamba install -c bioconda minimap2`)

## How to run

```bash
# smoke (fast)
bash polaplib/tooltest/minimap2/run.sh

# extended (tiny mapping)
bash polaplib/tooltest/minimap2/run.sh --extended

# with logging
LOG_FILE=o/tooltest/minimap2.log bash polaplib/tooltest/minimap2/run.sh --extended
```

```
## Exit codes

| Code | Meaning             |
| ---: | ------------------- |
|    0 | PASS                |
|    1 | FAIL                |
|   99 | SKIP (tool missing) |

## Expected output

- One or more `RUN:` lines and `PASS minimap2` on success.

## Troubleshooting

- `SKIP`: ensure `minimap2` is on `PATH`.
- `FAIL`: copy the **expanded** command after `;` in the `RUN:` line and execute it manually to see the raw error.

## Versions tested

| Tool     | Version | OS/Arch      | Result | Notes                    |
| -------- | ------- | ------------ | ------ | ------------------------ |
| minimap2 | 2.28    | linux/x86_64 | PASS   | smoke; extended tiny map |
```
