## `polaplib/tooltest/samtools/README.md`

````markdown
# samtools — smoke & extended tests

**Folder:** `polaplib/tooltest/samtools/`

## What this test covers

- **Smoke** (default): `samtools --version-only`
- **Extended** (`--extended`): minimal SAM → BAM → sort → index → idxstats

## Prerequisites

- `samtools` on `PATH` (e.g., `mamba install -c bioconda samtools`)

## How to run

```bash
# smoke (fast)
bash polaplib/tooltest/samtools/run.sh

# extended (tiny BAM workflow)
bash polaplib/tooltest/samtools/run.sh --extended

# with logging
LOG_FILE=o/tooltest/samtools.log bash polaplib/tooltest/samtools/run.sh --extended
```
````

## Exit codes

| Code | Meaning |
| ---: | ------- |
|    0 | PASS    |
|    1 | FAIL    |
|   99 | SKIP    |

## Expected output

- `RUN:` lines and `PASS samtools` on success.

## Troubleshooting

- `SKIP`: tool not installed / not on `PATH`.
- `FAIL`: copy **expanded** command and run manually for raw messages.

## Versions tested

| Tool     | Version | OS/Arch      | Result | Notes                    |
| -------- | ------- | ------------ | ------ | ------------------------ |
| samtools | 1.18    | linux/x86_64 | PASS   | smoke; extended pipeline |

```

```
