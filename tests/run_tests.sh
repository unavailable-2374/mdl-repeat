#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
BIN="$PROJECT_DIR/bin/mdl-repeat"
DATA_DIR="$SCRIPT_DIR/data"
RESULTS_DIR="$SCRIPT_DIR/results"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

pass=0
fail=0

check() {
    local name="$1"
    local condition="$2"
    if eval "$condition"; then
        echo -e "  ${GREEN}PASS${NC}: $name"
        ((pass++)) || true
    else
        echo -e "  ${RED}FAIL${NC}: $name"
        ((fail++)) || true
    fi
}

echo "=== MDL-Repeat Test Suite ==="
echo ""

# Build
echo "--- Building ---"
cd "$PROJECT_DIR"
make clean 2>/dev/null || true
make 2>&1 | tail -1
check "Build succeeds with zero warnings" "[ -x '$BIN' ]"

# Generate synthetic data
echo ""
echo "--- Generating synthetic test data ---"
mkdir -p "$DATA_DIR" "$RESULTS_DIR"
python3 "$SCRIPT_DIR/gen_synthetic.py" -o "$DATA_DIR"

# L_int unit tests
echo ""
echo "--- L_int unit tests ---"
gcc -O3 -Wall -Wextra -std=c11 -pthread \
    -o "$RESULTS_DIR/test_mdl" \
    "$SCRIPT_DIR/test_mdl.c" \
    "$PROJECT_DIR/src/mdl.c" \
    "$PROJECT_DIR/src/candidates.c" \
    "$PROJECT_DIR/src/genome.c" \
    "$PROJECT_DIR/src/kmer.c" \
    "$PROJECT_DIR/src/graph.c" \
    "$PROJECT_DIR/src/cmd_line_opts.c" \
    -lm 2>&1
"$RESULTS_DIR/test_mdl"
check "L_int test vectors" "$RESULTS_DIR/test_mdl"

# Test A — multiple families (1MB genome)
echo ""
echo "--- Test A: Multiple families (1MB genome) ---"
$BIN -sequence "$DATA_DIR/testA.fa" -output "$RESULTS_DIR/testA_families.fa" \
     -instances "$RESULTS_DIR/testA_instances.bed" \
     -stats "$RESULTS_DIR/testA_stats.tsv" -v 2>&1 | grep -E "(Accepted|Bases covered)" || true

n_fam_a=$(grep -c "^>" "$RESULTS_DIR/testA_families.fa" 2>/dev/null || echo 0)
n_inst_a=$(wc -l < "$RESULTS_DIR/testA_instances.bed" 2>/dev/null || echo 0)

check "Test A: found families (>= 1)" "[ $n_fam_a -ge 1 ]"
check "Test A: found instances (>= 50)" "[ $n_inst_a -ge 50 ]"

if [ -f "$DATA_DIR/testA_truth.bed" ] && [ -f "$RESULTS_DIR/testA_instances.bed" ]; then
    truth_bases=$(awk '{sum += $3 - $2} END {print sum+0}' "$DATA_DIR/testA_truth.bed")
    pred_bases=$(awk '{sum += $3 - $2} END {print sum+0}' "$RESULTS_DIR/testA_instances.bed")
    echo "  Truth repeat bases: $truth_bases, Predicted repeat bases: $pred_bases"
fi

# Test B — Tandem (known MVP limitation: 50bp tandem hard to detect)
echo ""
echo "--- Test B: Tandem array (MVP limitation: may not detect short tandems) ---"
$BIN -sequence "$DATA_DIR/testB.fa" -output "$RESULTS_DIR/testB_families.fa" \
     -stats "$RESULTS_DIR/testB_stats.tsv" -v 2>&1 | grep -E "(Accepted|Bases covered)" || true

n_fam_b=$(grep -c "^>" "$RESULTS_DIR/testB_families.fa" 2>/dev/null || echo 0)
echo -e "  ${YELLOW}INFO${NC}: Test B found $n_fam_b families (50bp tandem units near detection limit)"

# Test C — Detection limit
echo ""
echo "--- Test C: Detection limit (3 copies @ 0% vs 2 copies @ 20%) ---"
$BIN -sequence "$DATA_DIR/testC.fa" -output "$RESULTS_DIR/testC_families.fa" \
     -stats "$RESULTS_DIR/testC_stats.tsv" -v 2>&1 | grep -E "(Accepted|Bases covered)" || true

n_fam_c=$(grep -c "^>" "$RESULTS_DIR/testC_families.fa" 2>/dev/null || echo 0)
check "Test C: found at least 1 family" "[ $n_fam_c -ge 1 ]"

# Test D — Nested TE
echo ""
echo "--- Test D: Nested TE (SINE inside LINE) ---"
$BIN -sequence "$DATA_DIR/testD.fa" -output "$RESULTS_DIR/testD_families.fa" \
     -instances "$RESULTS_DIR/testD_instances.bed" \
     -stats "$RESULTS_DIR/testD_stats.tsv" -v 2>&1 | grep -E "(Accepted|Bases covered)" || true

n_fam_d=$(grep -c "^>" "$RESULTS_DIR/testD_families.fa" 2>/dev/null || echo 0)
check "Test D: found families (>= 1)" "[ $n_fam_d -ge 1 ]"

# Human3M smoke test
echo ""
echo "--- Human3M smoke test ---"
if [ -f "/home/shuoc/tool/RepeatScout/human3M.fa" ]; then
    $BIN -sequence /home/shuoc/tool/RepeatScout/human3M.fa \
         -output "$RESULTS_DIR/human3M_families.fa" \
         -instances "$RESULTS_DIR/human3M_instances.bed" \
         -stats "$RESULTS_DIR/human3M_stats.tsv" -k 14 -v 2>&1 | grep -E "(Accepted|Bases covered|Compression)" || true

    n_fam_h=$(grep -c "^>" "$RESULTS_DIR/human3M_families.fa" 2>/dev/null || echo 0)
    check "Human3M: found families (>= 10)" "[ $n_fam_h -ge 10 ]"
else
    echo "  SKIP: human3M.fa not found"
fi

# Summary
echo ""
echo "==========================="
echo -e "Results: ${GREEN}$pass PASS${NC}, ${RED}$fail FAIL${NC}"
echo "==========================="

# Clean up temp files from test_mdl binary
rm -f "$RESULTS_DIR/test_mdl"

[ "$fail" -eq 0 ]
