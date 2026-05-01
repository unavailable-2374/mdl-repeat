#!/bin/bash
# Explicit regeneration of synthetic test fixtures.
#
# tests/data/*.{fa,truth.bed} are committed to git as frozen fixtures so
# that test results stay comparable across runs.  Regenerating them is a
# deliberate act — it can change discovered family counts by changing the
# input data, even with deterministic Python seeds (Python's random module
# is not guaranteed bit-stable across major versions).
#
# Run this script ONLY when:
#   - You want to refresh fixtures after intentionally changing
#     gen_synthetic.py
#   - You are setting up a new test scenario
#
# After regenerating, commit the new fixtures and update baseline numbers
# in tests/baseline/*.txt as part of the same commit.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="$SCRIPT_DIR/data"

mkdir -p "$DATA_DIR"
python3 "$SCRIPT_DIR/gen_synthetic.py" -o "$DATA_DIR"

echo
echo "Done. Review changes with:"
echo "  git diff --stat tests/data/"
echo
echo "If they look right, commit them along with updated baseline files."
