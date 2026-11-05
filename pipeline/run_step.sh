#!/usr/bin/env bash
# Run a pipeline step and keep its stderr in the run's own Error folder (never loose in data/).
#   pipeline/run_step.sh pipeline/steps/05_filter.R --db=data/<prefix>/<prefix>.db --action=preview
#   pipeline/run_step.sh pipeline/steps/01_load.R --design=... --data=... --prefix=<prefix> [--outdir=data]
# stderr is shown as usual and also appended to <run_dir>/Error/<step>_stderr.log; the exit code is the step's.
set -o pipefail
script="$1"; shift || { echo "usage: $0 <step.R> [args...]" >&2; exit 2; }
run_dir=""; prefix=""; outdir="${DPDP_DATA:-$(cd "$(dirname "$0")/.." && pwd)/data}"
for a in "$@"; do
  case "$a" in
    --db=*) run_dir="$(dirname "${a#--db=}")" ;;
    --prefix=*) prefix="${a#--prefix=}" ;;
    --outdir=*) outdir="${a#--outdir=}" ;;
  esac
done
[ -z "$run_dir" ] && [ -n "$prefix" ] && run_dir="$outdir/$prefix"
step="$(basename "$script" .R)"
if [ -n "$run_dir" ]; then
  mkdir -p "$run_dir/Error"
  log="$run_dir/Error/${step}_stderr.log"
  printf '\n==== %s %s %s\n' "$(date '+%F %T')" "$step" "$*" >> "$log"
  Rscript "$script" "$@" 2> >(tee -a "$log" >&2)
  code=$?
  wait
  exit $code
else
  Rscript "$script" "$@"   # no run folder known (e.g. a tool without --db/--prefix): stderr only to the console
fi
