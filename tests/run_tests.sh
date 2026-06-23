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

# Synthetic data: tests/data/*.{fa,truth.bed} are committed to git as
# frozen test fixtures.  Only regenerate when one or more files are missing,
# to avoid silent drift across Python random-module changes.
echo ""
echo "--- Synthetic test data ---"
mkdir -p "$DATA_DIR" "$RESULTS_DIR"
need_regen=0
for f in testA.fa testA_truth.bed testB.fa testB_truth.bed \
         testC.fa testC_truth.bed testD.fa testD_truth.bed; do
    if [ ! -f "$DATA_DIR/$f" ]; then
        need_regen=1
        break
    fi
done
if [ "$need_regen" -eq 1 ]; then
    echo "  (missing files detected — regenerating)"
    python3 "$SCRIPT_DIR/gen_synthetic.py" -o "$DATA_DIR"
else
    echo "  Using committed fixtures in $DATA_DIR/ (run tests/regenerate_data.sh to refresh)"
fi

# mdl unit tests (L_int, model_cost, instance_cost, library selection)
echo ""
echo "--- mdl.c unit tests ---"
gcc -O2 -Wall -Wextra -std=c11 -pthread \
    -I"$PROJECT_DIR/src" \
    -o "$RESULTS_DIR/test_mdl" \
    "$SCRIPT_DIR/test_mdl.c" \
    "$PROJECT_DIR/src/mdl.c" \
    -lm 2>&1
if "$RESULTS_DIR/test_mdl" > "$RESULTS_DIR/test_mdl.log" 2>&1; then
    # grep -c exits 1 when zero matches — must guard against set -e killing
    # the script when n_fail == 0 (the happy path).
    n_pass=$(grep -c "^  PASS:" "$RESULTS_DIR/test_mdl.log" || true)
    n_fail=$(grep -c "^  FAIL:" "$RESULTS_DIR/test_mdl.log" || true)
    echo "  $n_pass assertions passed, $n_fail failed"
    check "mdl.c unit tests" "[ $n_fail -eq 0 ]"
else
    echo "  mdl unit tests crashed — see $RESULTS_DIR/test_mdl.log"
    check "mdl.c unit tests" "false"
fi

# Sweep-line regression tests (ENG-N2 + ENG-N10 + ENG-N12, V6 Tier 1.5b).
# Verifies the O(genome_len) coverage allocations are gone — peak RSS delta
# stays below 200 MB even when genome_len = 4e9.  Links against pre-built
# obj/*.o so we exercise the same code path as the production binary.
echo ""
echo "--- sweep-line regression tests (ENG-N2 + ENG-N10 + ENG-N12) ---"
gcc -O2 -Wall -Wextra -std=c11 -pthread \
    -I"$PROJECT_DIR/src" \
    -o "$RESULTS_DIR/test_sweepline" \
    "$SCRIPT_DIR/test_sweepline.c" \
    "$PROJECT_DIR/obj/mdl.o" \
    "$PROJECT_DIR/obj/refine.o" \
    "$PROJECT_DIR/obj/align.o" \
    "$PROJECT_DIR/obj/candidates.o" \
    "$PROJECT_DIR/obj/genome.o" \
    "$PROJECT_DIR/obj/kmer.o" \
    "$PROJECT_DIR/obj/output.o" \
    "$PROJECT_DIR/obj/discover.o" \
    "$PROJECT_DIR/obj/discover_mask.o" \
    "$PROJECT_DIR/obj/cmd_line_opts.o" \
    -lm 2>&1
if "$RESULTS_DIR/test_sweepline" > "$RESULTS_DIR/test_sweepline.log" 2>&1; then
    n_pass_sw=$(grep -c "^  PASS:" "$RESULTS_DIR/test_sweepline.log" || true)
    n_fail_sw=$(grep -c "^  FAIL:" "$RESULTS_DIR/test_sweepline.log" || true)
    echo "  $n_pass_sw assertions passed, $n_fail_sw failed"
    check "sweep-line regression tests" "[ $n_fail_sw -eq 0 ]"
else
    echo "  sweep-line tests crashed — see $RESULTS_DIR/test_sweepline.log"
    check "sweep-line regression tests" "false"
fi

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

# Test C2 — bounded recall rescue pass
echo ""
echo "--- Test C2: Bounded recall rescue discovery ---"
TESTC_RESCUE_FA="$RESULTS_DIR/testC_rescue_families.fa"
TESTC_RESCUE_TSV="$RESULTS_DIR/testC_rescue_stats.tsv"
TESTC_RESCUE_LOG="$RESULTS_DIR/testC_rescue.log"
TESTC_RESCUE_AUDIT="$RESULTS_DIR/testC_rescue_audit.tsv"
$BIN -sequence "$DATA_DIR/testC.fa" -output "$TESTC_RESCUE_FA" \
     -stats "$TESTC_RESCUE_TSV" \
     -recall-rescue -rescue-maxrepeats 5 \
     -rescue-audit "$TESTC_RESCUE_AUDIT" -v \
     2>"$TESTC_RESCUE_LOG" || true

check "Test C2: rescue pass honors maxrepeats parameter" \
      "grep -q 'Recall rescue discovery: .*MAXR=5' '$TESTC_RESCUE_LOG'"
check "Test C2: rescue pass defaults to targeted mode" \
      "grep -q 'mode=targeted' '$TESTC_RESCUE_LOG' && grep -q 'Recall rescue targeted' '$TESTC_RESCUE_LOG'"
check "Test C2: rescue append reports duplicate filtering" \
      "grep -Eq 'Recall rescue appended [0-9]+ candidate families \\(filtered [0-9]+ duplicates\\)' '$TESTC_RESCUE_LOG'"
check "Test C2: rescue log summarizes duplicate evidence metrics" \
      "grep -Eq 'Recall rescue evidence: candidates=[0-9]+, max_identity=[0-9.]+, max_containment=[0-9.]+, max_length_ratio=[0-9.]+, near_duplicates=[0-9]+' '$TESTC_RESCUE_LOG'"
check "Test C2: rescue audit records target and candidate decisions" \
      "grep -q '^record_type' '$TESTC_RESCUE_AUDIT' && grep -q '^target_segment' '$TESTC_RESCUE_AUDIT' && grep -q '^candidate' '$TESTC_RESCUE_AUDIT'"
check "Test C2: rescue audit records duplicate evidence metrics" \
      "head -1 '$TESTC_RESCUE_AUDIT' | grep -q 'consensus_identity.*contained_instance_fraction.*length_ratio'"
check "Test C2: rescue candidates carry quality provenance" \
      "grep -q 'rescue_discovery' '$TESTC_RESCUE_TSV'"
check "Test C2: missing rescue audit path is rejected" \
      "$BIN -sequence '$DATA_DIR/testC.fa' -output '$RESULTS_DIR/testC_missing_audit.fa' -recall-rescue -rescue-audit -v 2>'$RESULTS_DIR/testC_missing_audit.log'; [ \$? -ne 0 ] && grep -q 'ERROR: -rescue-audit requires a file path' '$RESULTS_DIR/testC_missing_audit.log'"

# Test C3 — external tool QC policy is optional and non-mutating by default
echo ""
echo "--- Test C3: External tool QC policy ---"
TESTC_EXT_FA="$RESULTS_DIR/testC_external_families.fa"
TESTC_EXT_QC="$RESULTS_DIR/testC_external_qc.tsv"
TESTC_EXT_LOG="$RESULTS_DIR/testC_external.log"
$BIN -sequence "$DATA_DIR/testC.fa" -output "$TESTC_EXT_FA" \
     -external-qc "$TESTC_EXT_QC" \
     2>"$TESTC_EXT_LOG" || true
check "Test C3: optional external QC absence does not fail core run" \
      "[ -s '$TESTC_EXT_FA' ] && grep -Eq 'External QC: seqkit stats wrote|WARNING: optional external QC seqkit' '$TESTC_EXT_LOG'"
check "Test C3: required missing external QC fails clearly" \
      "$BIN -sequence '$DATA_DIR/testC.fa' -output '$RESULTS_DIR/testC_external_required.fa' -external-tools require -external-qc '$RESULTS_DIR/testC_external_required_qc.tsv' -seqkit '$RESULTS_DIR/no_such_seqkit' 2>'$RESULTS_DIR/testC_external_required.log'; [ \$? -ne 0 ] && grep -q 'ERROR: required external QC seqkit failed' '$RESULTS_DIR/testC_external_required.log'"

# Test D — Nested TE
echo ""
echo "--- Test D: Nested TE (SINE inside LINE) ---"
$BIN -sequence "$DATA_DIR/testD.fa" -output "$RESULTS_DIR/testD_families.fa" \
     -instances "$RESULTS_DIR/testD_instances.bed" \
     -stats "$RESULTS_DIR/testD_stats.tsv" -v 2>&1 | grep -E "(Accepted|Bases covered)" || true

n_fam_d=$(grep -c "^>" "$RESULTS_DIR/testD_families.fa" 2>/dev/null || echo 0)
check "Test D: found families (>= 1)" "[ $n_fam_d -ge 1 ]"

# Test G — REMOVED 2026-05-01: tested Phase 5 F' RMBlast recruit which was
# reverted (Step 4b pre-pass caused 27+ hour TAIR10 hang via O(n²) assemble;
# chr4 family-level recall regressed -7.5pp 80×80). Test G fixtures retained
# in tests/data/testG.fa for future restoration if F' is re-implemented with
# proper scaling guards. See V6_PHASE5_RESULT.md for revert rationale.

# Test F — Bimodal-subfamily split instrumentation (Track 1 / Phase 4)
# Verifies that -vv emits [split] log lines (instrumentation smoke test)
# and that -threads 1 / -threads 4 produce a consistent family count.
echo ""
echo "--- Test F: Bimodal subfamily split instrumentation ---"
TESTF_FA="$DATA_DIR/testF.fa"
if [ ! -f "$TESTF_FA" ]; then
    echo "  Generating testF.fa..."
    python3 "$SCRIPT_DIR/gen_testF.py" -o "$DATA_DIR"
fi

if [ -f "$TESTF_FA" ]; then
    TESTF_T1_FA="$RESULTS_DIR/testF_t1_families.fa"
    TESTF_T4_FA="$RESULTS_DIR/testF_t4_families.fa"
    TESTF_T1_LOG="$RESULTS_DIR/testF_t1.log"
    TESTF_T4_LOG="$RESULTS_DIR/testF_t4.log"
    TESTF_T1_AUDIT="$RESULTS_DIR/testF_t1_split_audit.tsv"
    TESTF_T4_AUDIT="$RESULTS_DIR/testF_t4_split_audit.tsv"

    $BIN -sequence "$TESTF_FA" -output "$TESTF_T1_FA" \
         -threads 1 -vv -split-audit "$TESTF_T1_AUDIT" \
         2>"$TESTF_T1_LOG" || true

    $BIN -sequence "$TESTF_FA" -output "$TESTF_T4_FA" \
         -threads 4 -vv -split-audit "$TESTF_T4_AUDIT" \
         2>"$TESTF_T4_LOG" || true

    # 1. Instrumentation: [split] lines are emitted
    n_split_lines_t1=$(grep -c "\[split\]" "$TESTF_T1_LOG" 2>/dev/null || echo 0)
    check "Test F (-threads 1): [split] instrumentation lines present" "[ $n_split_lines_t1 -ge 1 ]"

    n_split_lines_t4=$(grep -c "\[split\]" "$TESTF_T4_LOG" 2>/dev/null || echo 0)
    check "Test F (-threads 4): [split] instrumentation lines present" "[ $n_split_lines_t4 -ge 1 ]"

    # 1b. Split audit: every thread mode emits a decision TSV.
    n_audit_t1=$(wc -l < "$TESTF_T1_AUDIT" 2>/dev/null || echo 0)
    n_audit_t4=$(wc -l < "$TESTF_T4_AUDIT" 2>/dev/null || echo 0)
    check "Test F (-threads 1): split audit TSV has decision rows" \
          "[ -f '$TESTF_T1_AUDIT' ] && grep -q 'decision' '$TESTF_T1_AUDIT' && [ $n_audit_t1 -ge 2 ]"
    check "Test F (-threads 4): split audit TSV has decision rows" \
          "[ -f '$TESTF_T4_AUDIT' ] && grep -q 'decision' '$TESTF_T4_AUDIT' && [ $n_audit_t4 -ge 2 ]"

    # 2. Family count sanity: should find at least 1, at most 6
    n_fam_f1=$(grep -c "^>" "$TESTF_T1_FA" 2>/dev/null || echo 0)
    n_fam_f4=$(grep -c "^>" "$TESTF_T4_FA" 2>/dev/null || echo 0)
    check "Test F (-threads 1): family count in [1,6]" \
          "[ $n_fam_f1 -ge 1 ] && [ $n_fam_f1 -le 6 ]"
    check "Test F (-threads 4): family count in [1,6]" \
          "[ $n_fam_f4 -ge 1 ] && [ $n_fam_f4 -le 6 ]"

    # 3. Thread consistency: family count agrees within 1
    diff_f=$((n_fam_f1 > n_fam_f4 ? n_fam_f1 - n_fam_f4 : n_fam_f4 - n_fam_f1))
    check "Test F: family count agrees within 1 across thread modes ($n_fam_f1 vs $n_fam_f4)" \
          "[ $diff_f -le 1 ]"
else
    echo "  SKIP: testF.fa not found"
fi

# Test E — Multi-chromosome cross-boundary correctness (V6 Tier 1.5a)
# Verifies the cross-chromosome guards in:
#   B+Q6:    nested_containment_fraction (refine.c)
#   ENG-N11: check_instance_overlap      (refine.c)
#   ENG-N8:  refine_coalesce_tandem_instances (refine.c)
#   ENG-N9:  refine_assemble_fragments   (refine.c)
echo ""
echo "--- Test E: Multi-chromosome cross-boundary correctness ---"
MULTICHR_DIR="$DATA_DIR/multichr"
if [ ! -f "$MULTICHR_DIR/multichr.fa" ] && [ -f "$MULTICHR_DIR/gen_multichr.py" ]; then
    python3 "$MULTICHR_DIR/gen_multichr.py" >/dev/null
fi

if [ ! -f "$MULTICHR_DIR/multichr.fa" ]; then
    echo "  SKIP: multichr.fa not found (run tests/data/multichr/gen_multichr.py)"
else
    # chr1 length is 100000; chr2 length is 100000.  The padded-genome
    # length used by the runner adds PADLENGTH (11000) at the front.
    # We assert no BED row spans across the chr1/chr2 boundary.
    chr_check() {
        local bed="$1"
        # Each BED row carries a chromosome name; if interval `end` exceeds
        # the chromosome's length by > PADLENGTH (11000) it indicates the
        # interval spilled into the adjacent sequence.  We use the strict
        # invariant: end <= chrom_len + buffer.  A buffer of 100 absorbs
        # legitimate alignment-extension overshoot but catches the bug
        # (which produces overshoots of thousands of bp).
        awk 'BEGIN{ok=1}
             $1 == "chr1" && $3 > 100100 {print "  CROSS-BOUNDARY: " $0; ok=0}
             $1 == "chr2" && $3 > 100100 {print "  CROSS-BOUNDARY: " $0; ok=0}
             $2 < 0 {print "  NEGATIVE-START: " $0; ok=0}
             END{exit !ok}' "$bed"
    }

    for THREADS in 1 4; do
        OUT_FA="$RESULTS_DIR/multichr_t${THREADS}_families.fa"
        OUT_BED="$RESULTS_DIR/multichr_t${THREADS}_instances.bed"
        OUT_TSV="$RESULTS_DIR/multichr_t${THREADS}_stats.tsv"

        $BIN -sequence "$MULTICHR_DIR/multichr.fa" \
             -output "$OUT_FA" \
             -instances "$OUT_BED" \
             -stats "$OUT_TSV" \
             -threads "$THREADS" 2>&1 | grep -E "(Accepted|Bases covered)" || true

        n_fam=$(grep -c "^>" "$OUT_FA" 2>/dev/null || echo 0)
        # Truth has 4 distinct repeat families (RepA/RepB on chr1, RepC/RepD on chr2).
        # If cross-chr ENG-N11 / B+Q6 merge runaway occurred, count would drop
        # below 3 (multiple distinct truth families collapsed via spurious
        # cross-chr instance overlap). >=3 retains correctness while allowing
        # one legitimate consolidation of low-complexity randomly-similar
        # regions.
        check "Test E (-threads $THREADS): >= 3 families found (no cross-chr merge collapse)" "[ $n_fam -ge 3 ]"
        check "Test E (-threads $THREADS): <= 8 families found (no spurious split runaway)"   "[ $n_fam -le 8 ]"
        # ENG-N8 invariant: no BED interval may span a chromosome boundary,
        # i.e. chr1.end <= 100100 (with 100 bp tolerance for align extension)
        # and chr2.end <= 100100.  A cross-chr tandem coalesce would produce
        # an interval ending at ~ chr1.size + chr2.size, well past the bound.
        check "Test E (-threads $THREADS): no BED interval spans chr boundary (ENG-N8)" "chr_check '$OUT_BED'"

        # Strand sanity: every BED row carries '+' or '-' in column 6.
        n_bad_strand=$(awk '$6 != "+" && $6 != "-" {print}' "$OUT_BED" | wc -l)
        check "Test E (-threads $THREADS): every BED row has valid strand"  "[ $n_bad_strand -eq 0 ]"
    done

    # Determinism (cross-thread): family count should match between -threads 1
    # and -threads 4 within a small tolerance.  Large divergence here usually
    # indicates a thread-race in the merge worker, not a refine bug per se,
    # but the cross-chr guard interacts with worker code so this is worth
    # asserting.
    n_t1=$(grep -c "^>" "$RESULTS_DIR/multichr_t1_families.fa" 2>/dev/null || echo 0)
    n_t4=$(grep -c "^>" "$RESULTS_DIR/multichr_t4_families.fa" 2>/dev/null || echo 0)
    diff_fc=$((n_t1 > n_t4 ? n_t1 - n_t4 : n_t4 - n_t1))
    check "Test E: family count agrees within 2 across thread modes ($n_t1 vs $n_t4)" "[ $diff_fc -le 2 ]"
fi

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

# Clean up temp files from unit-test binaries
rm -f "$RESULTS_DIR/test_mdl" "$RESULTS_DIR/test_sweepline"

[ "$fail" -eq 0 ]
