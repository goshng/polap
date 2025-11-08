Here’s a practical, copy-pasteable guide to make your Bash, Python, and R scripts “failsafe-ready” so they’re debuggable with polaplib/polap-lib-failsafe.sh v0.3.0 (the quiet, reliable one you’re using).

⸻

A) Bash scripts (parents and children)

1. Minimal header for every Bash process you want captured

Put this at the top, before any code (and before any source):

#!/usr/bin/env bash

# — re-exec in bash if invoked via /bin/sh —

had_u=0; case $- in *u*) had_u=1;; esac; set +u
if [ -z "${BASH_VERSION+x}" ]; then exec /usr/bin/env bash "$0" "$@"; fi
[ "$had_u" -eq 1 ] && set -u

# — strict + trap propagation —

set -Eeuo pipefail
set -o errtrace
set -o functrace

# — enable failsafe —

source "${\_POLAPLIB_DIR}/polap-lib-failsafe.sh" # v0.3.0 baseline
polap_enable_failsafe

# re-exec in bash if invoked via /bin/sh

had_u=0; case $- in *u*) had_u=1;; esac; set +u
if [ -z "${BASH_VERSION+x}" ]; then exec /usr/bin/env bash "$0" "$@"; fi
[ "$had_u" -eq 1 ] && set -u

# strict + trap propagation

set -Eeuo pipefail
set -o errtrace
set -o functrace

# enable failsafe

source "${\_POLAPLIB_DIR}/polap-lib-failsafe.sh"
polap_enable_failsafe

Why:
• set -Eeuo pipefail makes failures bubble out.
• errtrace/functrace makes ERR fire inside functions/subshells.
• The re-exec shim guarantees we’re in bash (not dash/sh).
• polap_enable_failsafe installs a quiet DEBUG snapshot + the ERR printer.

2. Keep failures visible (so ERR actually fires)

Avoid these patterns (they mask errors or suppress ERR):
• if some_cmd; then ... fi
• some_cmd || true
• ! some_cmd
• ( some_cmd ) (subshell that swallows status)
• pipelines without pipefail (you already set -o pipefail)

Prefer simple commands that fail plainly, e.g.:

bash -c "exit 7" # good: triggers ERR
false # good: triggers ERR

3. Make parent wrappers thin and honest

When launching a child, do just enough:

# Optional: human-friendly raw + expanded one-liner (so you can copy/paste):

printf '[%s:%s] parent: %s\n' "${BASH_SOURCE[0]##*/}" "${LINENO}" \
 'python3 "$child_py" --message "$message" --code "$code"'

expanded=( python3 "$child_py" --message "$message" --code "$code" )
printf '[%s:%s] parent: %s\n' "${BASH_SOURCE[0]##\*/}" "${LINENO}" \
  "$(printf '%q ' "${expanded[@]}")"

# Actual run (let non-zero bubble so ERR prints our call-site)

python3 "$child_py" --message "$message" --code "$code"

If you want the parent crash line to include (conda-env) and an expanded command automatically, add that logic to your failsafe printer or wrapper. (You already experimented with it.)

4. Don’t clear traps (unless you restore them)

Avoid trap - DEBUG / trap - ERR without restoring; the failsafe relies on them.

⸻

B) Python scripts (children)

1. Minimal crash helper template

#!/usr/bin/env python3
import argparse, sys, traceback

class ExitCodeError(RuntimeError):
def **init**(self, msg, code): super().**init**(msg); self.code = int(code)

def \_hook_trace(typ, exc, tb):
traceback.print_exception(typ, exc, tb)
code = int(getattr(exc, "code", 1)) or 1
sys.exit(code)

def \_hook_quiet(typ, exc, tb):
print(str(exc) or typ.**name**, file=sys.stderr)
code = int(getattr(exc, "code", 1)) or 1
sys.exit(code)

p = argparse.ArgumentParser()
p.add_argument("--message","-m", default="intentional test failure (python)")
p.add_argument("--code","-c", type=int, default=7)
p.add_argument("--exit-code", dest="code_alt", type=int)
p.add_argument("--in", dest="infile", default="")
p.add_argument("--out", dest="outfile", default="")
p.add_argument("--no-trace", dest="no_trace", action="store_true")
a = p.parse_args()

code = a.code_alt if a.code_alt is not None else a.code
if code == 0: code = 7
sys.excepthook = \_hook_quiet if a.no_trace else \_hook_trace

print(f"PY test_crash: message='{a.message}' code={code} infile='{a.infile}' outfile='{a.outfile}'",
file=sys.stderr)

if a.outfile:
try:
with open(a.outfile, "w", encoding="utf-8") as fh:
fh.write("polap-py-test-crash wrote before failing\n")
except Exception as e:
print(f"[child] warn: {e}", file=sys.stderr)

def level3(): raise ExitCodeError(a.message, code)
def level2(): level3()
def level1(): level2()
level1()

Why:
• Always raises an exception → non-zero exit.
• Traceback (or quiet) printing is controlled via --no-trace.
• The parent Bash failsafe prints where you called Python; Python prints what failed.

⸻

C) R scripts (children) with optparse (your requirement)

1. Simple, robust optparse + stack on by default

#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
options(show.error.messages = FALSE) # handler controls printing

opt_list <- list(
make_option(c("-m","--message"), type="character", default="intentional test failure (R)"),
make_option(c("-c","--exit-code"), type="integer", default=7, dest="exit_code"),
make_option(c("--in"), type="character", default="", dest="infile"),
make_option(c("--out"), type="character", default="", dest="outfile"),
make_option(c("--no-stack"), action="store_true", default=FALSE, dest="no_stack")
)
parser <- OptionParser(usage = "%prog [options]", option_list = opt_list)
opts <- parse_args(parser)

msg <- opts$message
ecode <- if (!is.null(opts$exit_code) && !is.na(opts$exit_code)) opts$exit_code else 7L
if (ecode == 0L) ecode <- 7L

# short stack printer (file:line when srcref exists; otherwise calls only)

stack_handler <- function() {
calls <- sys.calls(); nc <- length(calls)
if (!nc) { cat("Stack (empty)\n"); quit(save="no", status=ecode) }
cat("Stack (oldest → newest):\n")
for (i in seq_len(nc)) {
cl <- calls[[i]]; sr <- attr(cl,"srcref")
txt <- paste(deparse(cl), collapse=" ")
if (!is.null(sr)) {
sf <- attr(sr,"srcfile"); fn <- if (!is.null(sf) && !is.null(sf$filename)) basename(sf$filename) else "<src>"
cat(sprintf("#%d %s:%d:%d %s\n", i, fn, as.integer(sr[[1]]), as.integer(sr[[2]]), txt))
} else {
cat(sprintf("#%d %s\n", i, txt))
}
}
quit(save="no", status=ecode)
}
quiet_handler <- function() quit(save="no", status=ecode)
options(error = if (!isTRUE(opts$no_stack)) stack_handler else quiet_handler)

cat(sprintf("R test_crash: message='%s' code=%d infile='%s' outfile='%s'\n",
msg, as.integer(ecode), opts$infile, opts$outfile))

if (nzchar(opts$outfile)) {
  con <- file(opts$outfile, open="w")
writeLines("polap-r-test-crash wrote before failing", con)
close(con)
}

third <- function(m,c) stop(m, call.=TRUE)
second <- function(m,c) third(m,c)
first <- function(m,c) second(m,c)
first(msg, as.integer(ecode))

Why:
• optparse is kept and used correctly (no nested $options objects, no conflicting dest).
• Stack is ON by default; pass --no-stack to silence.
• Under plain Rscript, srcref can be missing → the stack prints calls without line numbers (as you accepted). That’s fine; the Bash failsafe still prints the parent call site accurately.

If you must have file:line from R under Rscript, use the earlier self-resourcing trick and accept the complexity, or run via Rscript -e 'options(keep.source=TRUE); source("file.R")'.

⸻

D) Parent vs Child: who enables failsafe?
• Every Bash process whose call-site you want printed must enable the failsafe.
Parents: print the line that launched Python/R.
Children (Bash): print their own failing line.
• Language children (Python/R) do not use the Bash failsafe; they print their own traceback/stack. The Bash failsafe prints the call site in the wrapper.

⸻

E) Quick testing recipe 1. Bash child

bash polaplib/polap-bash-test-crash-bash.sh --mode cmd --code 7

    2.	Python child via Bash parent

bash polaplib/polap-bash-test-crash-python.sh --message "boom" --code 7

    3.	R child via Bash parent (stack on by default)

bash polaplib/polap-bash-test-crash-r.sh --message "boom" --code 7

(Add --no-stack if you want to silence R’s stack.)

⸻

TL;DR
• Bash: add the strict header + polap_enable_failsafe in every process; avoid masking failures; use thin wrappers.
• Python: raise exceptions; choose traceback vs quiet with a flag; let non-zero exit reach the parent.
• R (optparse): stack on by default (use --no-stack to silence); srcref optional; exit with code.

With those edits, your scripts are fully “failsafe-aware”: Bash immediately tells you where and what failed in the shell, and Python/R tell you why inside their own stacks—exactly the division of labor you want.
