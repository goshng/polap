## `polaplib/tooltest/iqtree2/README.md`

````markdown
# IQ-TREE 2 — smoke & extended tests

**Folder:** `polaplib/tooltest/iqtree2/`

## What this test covers

- **Smoke** (default): `iqtree2 --version` (fallback `iqtree -h`)
- **Extended** (`--extended`): tiny alignment → tree inference (GTR+G, -nt 1, -seed 1)

## Prerequisites

- `iqtree2` (or `iqtree`) on `PATH` (e.g., `mamba install -c bioconda iqtree`)

## How to run

```bash
# smoke
bash polaplib/tooltest/iqtree2/run.sh

# extended
bash polaplib/tooltest/iqtree2/run.sh --extended
```
````

## Exit codes

| Code | Meaning |
| ---: | ------- |
|    0 | PASS    |
|    1 | FAIL    |
|   99 | SKIP    |

## Expected output

- `RUN:` lines; extended mode writes outputs under `o/tooltest/iqtree2/demo.*`.

## Troubleshooting

- `SKIP`: neither `iqtree2` nor `iqtree` found.
- `FAIL`: run the expanded command printed after `;` to inspect.

## Versions tested

| Tool    | Version | OS/Arch      | Result | Notes                        |
| ------- | ------- | ------------ | ------ | ---------------------------- |
| iqtree2 | 2.3.x   | linux/x86_64 | PASS   | smoke; tiny tree in extended |

```

```
