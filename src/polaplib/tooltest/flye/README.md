## `polaplib/tooltest/flye/README.md`

````markdown
# Flye — smoke & optional extended tests

**Folder:** `polaplib/tooltest/flye/`

## What this test covers

- **Smoke** (default): `flye --version` and `flye --help` first lines
- **Extended** (`--extended`): optional tiny run (likely to fail by design on miniature data, intended to exercise error paths)

## Prerequisites

- `flye` on `PATH` (e.g., `mamba install -c bioconda flye`)

## How to run

```bash
# smoke
bash polaplib/tooltest/flye/run.sh

# extended (may fail by design; check logs)
bash polaplib/tooltest/flye/run.sh --extended
```
````

## Exit codes

| Code | Meaning                 |
| ---: | ----------------------- |
|    0 | PASS                    |
|    1 | FAIL                    |
|   99 | SKIP (tool not present) |

## Expected output

- `RUN:` lines and `PASS flye` for smoke; extended prints rc and leaves traces under `o/tooltest/flye/asm/` if it runs.

## Troubleshooting

- `SKIP`: not installed or not on `PATH`.
- For extended, a non-zero exit is acceptable on tiny data; inspect logs.

## Versions tested

| Tool | Version | OS/Arch      | Result | Notes      |
| ---- | ------- | ------------ | ------ | ---------- |
| flye | 2.9     | linux/x86_64 | PASS   | smoke only |

```

```
