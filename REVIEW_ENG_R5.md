## Round 5 Engineering Review

### Integration check (R4 items in v5)

**ENG-N12 (mdl.c:247-248 bundled with ENG-N2 + N10): correctly integrated in body.**

The bundled section "ENG-N2 + N10 + N12" correctly names both allocation sites, gives sizes, states "REGRESSION test must cover BOTH allocations", effort updated 1.0 → 1.5 days, silent-discard failure mode at mdl.c:249 called out.

**Citation correction (check_instance_overlap call sites: 496, 532, 696, 787): correctly integrated.**

Verified by `grep -n 'check_instance_overlap' src/refine.c`: definition at line 252, calls at 496, 532, 696, 787.

### Documentation inconsistencies found in v5

These were introduced in the v4 → v5 edit. They do not change implementation requirements but could mislead the implementer.

**Doc-1: Phase 2 execution order still lists "ENG-N2 + ENG-N10 (1 day)" — omits N12 and shows wrong effort.**

Body says ENG-N2+N10+N12 bundled at 1.5 days; execution order step 9 says "ENG-N2 + ENG-N10 (1 day)". Inconsistent.

**Doc-2: ENG-N9.aux appears in Phase 1 execution order (step 4) and test table despite being marked "completed inline / no further audit needed" in body.**

Self-contradictory.

**Doc-3: ENG-N9.aux result text mislabels 7 sites as "5 unsafe."**

`5 unsafe: B+Q6 (lines 225, 230), ENG-N11 (lines 268, 275, 284, 307, 311)` — that's 7 line numbers under the "5 unsafe" label. Should split into "2 covered by B+Q6" + "5 covered by ENG-N11".

### O(genome_len) allocation sweep — no new sites

Confirmed against live grep: only refine.c:1769 and mdl.c:248 are O(genome_len). Other allocations are O(num_instances), O(num_families), or constant. Workspace arrays in discover.c use MAXN constant.

### NEW action items (Round 5 only)

**No new engineering action items in Round 5. This is round 1 of "no new" implementation-correctness items.**

The three documentation issues above (Doc-1/2/3) are not ENG-level items — they affect only the proposal document, not source code. The author should fix them during the next proposal edit; no review round is required to gate that.

### Final assessment

R5 = 0 new code-level items. Doc fixes are housekeeping. Round 6 needed to satisfy 2-consecutive-round criterion under strict reading.
