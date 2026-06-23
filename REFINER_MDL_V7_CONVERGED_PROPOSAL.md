# Refiner_mdl Redesign — v7 Converged Proposal

**Status**: 8-round iterative review complete. **Bio + Eng dual-track convergence**.
**Goal**: reduce false positives in mdl-repeat output AND improve consensus quality. Output quality dominant; no concern for implementation cost.
**Empirical baseline** (TAIR10 nuclear, family-level recall): mdl-repeat 0.821 / 0.661 (80×80 / 90×80). v7 designed to ≥ 0.85 / ≥ 0.70.

---

## 0. Iteration Journey (8 rounds)

| Round | Bio findings | Eng findings | v# produced | Verdict |
|---|---:|---:|---|---|
| R1 | 33 | 30 | v2 | Major architectural rewrite |
| R2 | 19 | 20 | v3 | Layer 0/2 made continuous; layer-cascade redesigned |
| R3 | 12 | 22 (5 Critical) | v4 | Output rewrite + 3-track training + tool batch |
| R4 | 7 | 15 | v5 | Tool selection tightened; calibration genomes; provenance |
| R5 | 4 | 13 | v6 | EC artifact; bit_per_residue; 3-state CI |
| R6 | 1 (Minor) | **0** ✅ | v7 | Eng clean round 1/2; Bio added HELIANO Flag 4 |
| R7 | **0** ✅ | **0** ✅ | v7 (no change) | Eng CONVERGED 2/2; Bio clean round 1/2 |
| R8 | **0** ✅ | **0** ✅ | v7 (locked) | **BOTH CONVERGED 2/2** |

**Total findings collected & integrated**: 76 Bio + 100 Eng = **176 unique design improvements** across 7 versions.

**Convergence**: monotone-decreasing trajectory. Combined Bio+Eng findings: 63 → 39 → 34 → 22 → 17 → 1 → 0 → 0.

---

## 1. v7 ARCHITECTURE — 7-layer evidence-collection pipeline (LOCKED)

```
mdl-repeat output (FASTA + BED + stats)
     │
     ▼
LAYER 0 — Pre-filter EVIDENCE (continuous scores; NO boolean drops)
     0.1  TRF tandem repeat                    → trf_score, trf_fraction
     0.2  rDNA blastn (5S+18S+5.8S+25S+IGS)   → rdna_pident, rdna_coverage
     0.3  Rfam cmscan tRNA/snoRNA/snRNA       → rfam_best_e, rfam_best_family
     0.4  Telomeric motif (TTTAGGG/CCCTAAA)   → telomere_fraction
     0.5  TE-purged SwissProt blastx          → host_gene_evalue, host_gene_coverage
     0.6  NUMT/NUPT (minimap2 vs mt+cp)       → numt_pident, nupt_pident
     0.7  Pack-TYPE umbrella (6 superfamilies + putative Pack-PiggyBac flag):
          Pack-MULE, Pack-CACTA, Pack-hAT, Pack-Harbinger/PIF, Pack-Mariner/Tc1, Helitron-NLR fusion
     ─→ continuous evidence dictionary per family

LAYER 1 — HMM-based primary classification
     1.1  hmmscan vs Dfam 3.8 (Zenodo SHA256, hmmpress'd plant subset)
     1.2  hmmscan vs Pfam TE profiles (PF00078 RT, PF00665 IN, PF03017 LINE_endo,
          PF00872 DDE_Tnp, PF14223 Helitron_helicase, ...)
     1.3  RepeatClassifier (RepeatModeler 2.0.6+)
     1.4  TEsorter REXdb HMMs
     1.5  NeuralTE (2024) / Terrier (2025 PMID 40862518) — auxiliary signals only
     1.6  Dfam-to-Wicker mapping post-hmmscan

LAYER 2 — Class-specific structural validators (parallel cascades, EDTA-style; emit continuous confidence)
     2.1  LTR-RT (LTR_retriever architecture):
          - LTR self-blastn ≥100 bp; canonical TG/CA + non-canonical TGTT/AACA
          - LTR-pair identity ≥85%; TSD 4-6 bp from BED flanks
          - LARD/TRIM sub-classification (internal<500bp → TRIM_candidate;
                                          internal>3000bp + zero Pfam TE → LARD_candidate)
          - tandem-LTR-RT flag
     2.2  LINE / non-LTR retro:
          - Layer 1 ORF2 RT+EN domain hit; allow 5'-truncated copies (no ORF1 required)
          - 3' poly(A) ≥10 A's; length ≥500 bp
     2.3  SINE (PRIMARY: AnnoSINE2 2024; SECONDARY: motif scan):
          - tRNA-head blastn vs tRNAscan-SE (≥70% over ≥30 bp)
          - Degenerate `RVTGG` + `GTTCRA`, 24-45 bp spacing
          - SINE3 5S rRNA promoter (A-box + IE-box + C-box)
          - snRNA-SINEU (taxon-conditional, --kingdom animal)
          - 3' poly(A); length 80-500 bp
     2.4  TIR DNA TE (PRIMARY: MITE-Tracker; HelitronScanner removed):
          - TIR 10-80 bp (raised upper from MITE-Hunter 40)
          - Class-specific TSDs: Tc1/Mariner 2bp TA; Stowaway 2bp TA;
            PIF/Harbinger (incl. Tourist-MITE, Pong) 3bp TAA/TTA;
            hAT 8bp; Mutator/MULE 9-11bp; CACTA 3bp; Helitron NONE (explicit)
     2.5+2.6  Helitron HLE1 + HLE2 (HELIANO 2024 unified):
          4 QC flags:
            Flag 1: N_nonauto < 0.5 × N_auto → underdetection
            Flag 2: N_auto == 0 + Dfam ancestor has Helitron → autonomous_zero
            Flag 3: N_orfonly / (N_auto + N_orfonly) > 0.20 → boundary_unresolved
            Flag 4: HLE2 proto-Helentron OR length > 3× median → 5p_boundary_uncertain
          On Flag 3/4: suppress L4 TSD positional check; "structure_evidence_unavailable" in L5 Track B
     2.7  DIRS-order (4-superfamily disambiguation):
          Step 1: MT domain (PF05869/PF13578) present → DIRS-like
          Step 2 (MT absent): terminal architecture
            ITR + ICR → DIRS-like (some lack MT)
            SDR (A1-B1-A2-B2) → continue
          Step 3: TEsorter REXdb sub-clade phylogeny:
            Ngaro-clade → Ngaro
            PAT-clade → PAT-like
            VIPER-clade → "DIRS-order/VIPER-like (kinetoplastid-associated)"
            Else → "DIRS_order_unresolved"
     2.8  PLE (animal): GIY-YIG endonuclease + telomerase-like RT
     2.9  Polinton/Maverick (animal): pPolB + retroviral IN + 6 bp TSD
     2.10 Crypton (--kingdom non-plant): tyrosine recombinase Pfam PF02899
     2.11 CRM-clade Gypsy: chromodomain PF01108; centromeric_gypsy_candidate tag

LAYER 3 — Copy-number + dispersion EVIDENCE
     3.1  Within-dataset percentile (NOT magic numbers): copy_pctile, mdl_pctile
     3.2  Dispersion (assembly-quality conditional):
          N50 ≥1Mb: ≥2 chrs OR ≥1Mb apart on same chr
          N50 <1Mb: Gini coefficient > 0.4 OR ≥3 distinct scaffolds with ≥2 copies each
                    AND scaffold_length ≥5× consensus_length
     3.3  Class-aware copy ranges (relaxed bounds; flags only, not exclusions):
          LTR-RT 10-3000; LINE 10-1000; SINE 50-10000; MITE 10-10000;
          Helitron 10-3000; CRM 5-1000

LAYER 4 — TSD evidence from BED flanks
     - Per-class TSD validation
     - ≥80% of instances must agree on TSD length/seq
     - Helitron explicitly exempt from TSD requirement
     - Flag 4 suppresses positional check for proto-Helentron candidates

LAYER 5 — Probabilistic classifier (calibrated)
     5.1  4-track training (breaks Pfam leak):
          Track A: Pfam-discovered TE proteins (RT/IN/RH/transposase)
          Track B: Structure-discovered TIR/TSD/hairpin (MITEs, Helitrons)
          Track C: Copy-number-discovered (TAIR10 manually curated low-identity)
          Track D: Pack-TYPE (host-protein Pfam features kept but is_pack_type flagged)
          k-fold stratified by discovery method
          Cross-track validation: train on B+C+D, test on A; etc.

     5.2  Feature equivalence class (EC) ablation:
          - Ship feature_equivalence_classes.tsv via git-lfs
          - ~25-30 ECs: EC_RT, EC_IN, EC_RNase_H, EC_HELITRON_HELICASE, EC_DDE_TNP,
            EC_GAG, EC_TYR_RECOMBINASE, EC_METHYLTRANSFERASE, EC_PRIMASE_POLB, ...
          - Track A k-fold ablates entire EC, not individual feature
          - CI test: every Pfam accession in features.json maps to exactly one EC

     5.3  Features:
          Layer 0-4 outputs + bit_per_residue scores (HMMER null-model length-corrected)
                  + auxiliary (HMM_length, target_length)
          Family-internal: length_cv, gc_delta, inter_copy_distance_entropy,
                           dispersion_chr_count, mdl_pctile, copy_pctile
          consensus_source ∈ {hmmemit, raw_mdl_pre_refine, raw_mdl_post_refine,
                              ltr_pair_consensus, famsa_only_no_hmm, cdhit_centroid}

     5.4  Model: LightGBM 4.5.0 (pinned + model.txt SHA256 git-lfs)
          + SHAP 0.46 (deterministic TreeSHAP)
          + isotonic calibration (NOT Platt; reliability diagram + Brier + ECE)
          + ECE block release if any per-track or per-consensus_source ECE > 0.07

     5.5  TOOL_MISSING handling:
          Required tool MISSING → hard-fail at pipeline start
          Optional tool MISSING → GBT native missingness (np.nan)
          --strict-evidence-contract flag for production safety
          --require-min-classifiers N (default 2)

     5.6  6-genome calibration panel:
          TAIR10 / rice / maize / Drosophila / C. elegans / wheat-chr3B subset
          (Full wheat 17Gb only for release validation, not iteration)

     5.7  Output: P(is_TE) ∈ [0,1] calibrated probability + per-family TreeSHAP top-3 features

LAYER 6 — Consensus polish (HMM-based)
     6.1  Orient sequences via blastn vs consensus FIRST (avoids MAFFT --adjustdirection bug)
     6.2  Recruit copies: BED instances PRIMARY; minimap2 -x asm20 SECONDARY at 70% identity
     6.3  MSA: FAMSA + Kalign 3 fallback for short cases
     6.4  hmmbuild → hmmemit -c (principled match/insert handling)
     6.5  Subfamily detection (NEW per F2-7):
          PRIMARY: COSEG (Goubert et al. 2022 standard); 90/90 fallback
          PRE-CHECK: truncation screen (>30% instances share boundary → cd-hit instead)
          nhmmer self-rescan as sanity check
          Recursion: depth cap 2; min 5 copies/leaf; dot-separated naming (R=42.a.b)
     6.6  Length anchor: ±10% of MEDIAN recruited copy length
     6.7  5'-truncated LINE special path (Castaño-Basques 2022 / Goubert 2022 fallback citation):
          Build consensus from 3' ~900 bp first; extend 5' from minority full-length copies
     6.8  Skip MSA: copies<3 OR length<150 (CLI-overridable)

LAYER 7 — Output
     ONE ranked FASTA with calibrated P + Wicker class + tier label + SHAP top-3 + all classifications + disagreement_flag in header:
          >R=42|wicker=LTR/Copia/Athila|P=0.872|n_copies=14|tier=high|consensus_source=hmmemit|
           HELIANO=NA|RC=Unknown|TEsorter=Athila|NeuralTE=Helitron(0.62)|disagreement_flag=False|
           shap_top1=pfam_RT_bit_per_res(0.42)|shap_top2=ltr_pair_identity(0.31)|shap_top3=copies(0.18)
     + sidecar {prefix}.ranked.tsv (sorted by P descending)
     + provenance.yaml byte-identity (all hashes/versions/seeds/thresholds)
     + tier presets: high (P≥0.95), medium (P≥0.80), low (P≥0.50) — annotations, NOT file boundaries
     + optional --filter-threshold N flag for downstream filtering
```

---

## 2. Engineering replacements (LOCKED)

| Component | v7 final |
|---|---|
| Workflow | Nextflow DSL2 + nf-core template |
| Containers | per-process BioContainers + Apptainer SIF for HPC; **batch-per-tool** |
| Dedup | mmseqs2 linclust --threads 1 --shuffle 0 --cov-mode 0 -c 0.8 (deterministic) |
| Two-tier | --repro-strict (--threads 1) vs default (--threads N); CI asserts ≤0.5pp recall delta |
| MSA | FAMSA + explicit pre-orientation + Kalign 3 fallback for short |
| Recruit | minimap2 -x asm20 |
| Helitron | HELIANO 2024 (NOT HelitronScanner 2014) |
| MITE | MITE-Tracker 2018 + RepeatModeler2 integrated (NOT MITE-Hunter 2010) |
| SINE | AnnoSINE2 2024 (NOT regex motif) |
| HMMER | 3.4 pinned (Aug 2023) |
| Dfam | 3.8 pinned (Oct 2023, Zenodo SHA256, hmmpress'd plant subset) |
| LightGBM | 4.5.0 pinned + model.txt git-lfs SHA256 |
| SHAP | 0.46 pinned |
| Calibration | mandatory isotonic (NOT Platt for imbalanced) |
| Library scan | nhmmer subsampled 10Mb genome + RepeatMasker (full at end) |
| LTR_retriever | separate Apptainer SIF; canonical manifest TSV for cross-SIF hand-off |

---

## 3. Three-state evidence contract (CI-enforced)

Every layer emits one of three states per evidence field:
- `RAN_PASSED`: tool ran, evidence collected
- `RAN_FAILED`: tool ran, no positive evidence (real negative)
- `TOOL_MISSING`: tool absent or failed to start

CI test (tests/eng/test_evidence_contract.py):
- Fixture with each state per feature
- Assert LightGBM encodes TOOL_MISSING as np.nan, NOT 0.0
- Assert RAN_FAILED ≠ TOOL_MISSING in calibration accounting

`--strict-evidence-contract` CLI: hard-fail if required tool has TOOL_MISSING.

Per CRITICAL_DATA_INTEGRITY: NEVER simulate, mock, or impute missing evidence.

---

## 4. All thresholds in YAML config (citation-tagged)

`config/thresholds.yaml`:
- ~25 entries with `value`/`unit`/`source`/`justification`
- CLI overrides exposed via Pydantic-validated CLI args
- CI static check: numeric literals outside config-loader → hard-fail

---

## 5. Reproducibility commitment (LOCKED)

| Field | Commitment |
|---|---|
| Family count + ID set | byte-identical across linux-amd64/linux-arm64 |
| Calibrated P values | agree to 3 decimals across architectures |
| Per-family Wicker class | byte-identical |
| Per-family consensus FASTA | may differ in MSA gap placement ≤1% positions |
| Provenance YAML main block | SHA256-identical across runs (excluding `runtime` block) |

---

## 6. Diagnostic experiments (must run BEFORE deployment)

| ID | Experiment | Question |
|---|---|---|
| **D12** | Feasibility prototype REQUIRED before full build | +3pp/+4pp recall achievable? Go/no-go: 95% CI lower>0pp + point>1pp on TAIR10 + ≥1 other genome |
| D1 | 4-cell ablation on Phase 1 | Which polish component drives recall change |
| D2 | cd-hit dedup ablation | Is dedup the source of recall loss |
| D3 | EDTA cascade benchmark | SOTA comparator |
| D4 | Per-detector FPR on shuffled control | Calibrate per-test alpha |
| D5 | Layer 5 LR vs GBT cross-validation | Choose model |
| D6 | Per-class recall (Helitron / MITE / SINE / LTR-RT / LINE / Helentron) | Detect class-specific regression |
| D7 | Earl Grey 2024 head-to-head benchmark | Contemporary plant TE pipeline comparison |
| D8 | Layer 5 calibration LR vs isotonic-calibrated GBT | Calibrate decision |
| D9 | TE-purged SwissProt loss validation | Confirm TE-purge doesn't lose true LINEs/LTRs |
| D10 | NUMT/NUPT detection sensitivity (TAIR10 + wheat) | Validate Layer 0.6 |
| D11 | Subfamily-split rate (Layer 6.6) on synthetic chimeric | Validate COSEG truncation pre-check |

---

## 7. Tool naming + migration path

**v7 = NEW separate tool**, name TBD (e.g., `te-polish-v7` or `refiner-mdl-v7`):
- Consumes mdl-repeat ≥ v6.1 FASTA + BED + stats
- Backward compat: parses v6.0 (no div=/topo=) and v6.1+ (with div=/topo=) headers
- Produces ranked FASTA + ranked.tsv + provenance.yaml
- Does NOT modify mdl-repeat in any way
- Independent versioning

---

## 8. Empirical context

**Why this redesign was needed** (per `REFINER_MDL_BENCHMARK_RESULT.md`):

The original Refiner_mdl tested on TAIR10 nuclear (mdl-repeat → Refiner_mdl pipeline) destroyed recall:
- mdl-repeat baseline: 0.821 / 0.661 (80×80 / 90×80)
- Refiner_mdl analysis lib (full pipeline): 0.231 / 0.181 (**-59 pp at 80×80**)
- Cause: Phase 2 te_structure_filter dropped 12522→854 sequences (93% loss), incompatible with Helitron/MITE/SINE/novel families

v7 redesign:
- ALL Layer 0/2 evidence is continuous (no boolean drops)
- Layer 5 calibrated probabilistic classifier integrates evidence
- Pack-TYPE and gene-capture cases handled
- HELIANO/AnnoSINE2/MITE-Tracker (modern 2024 tools) replace stale 2010-2014 tools
- Single ranked FASTA output with P(is_TE) per family

**v7 design target**: ≥ 0.85 / ≥ 0.70 (improve over baseline by 3-4 pp).

---

## 9. Reviewer convergence record

After 8 rounds with `bio-interpreter` (TE biology) + `general-purpose` (bioinformatics engineering) agents:

| Reviewer | Final verdict |
|---|---|
| Bio (R8) | "**NO NEW BIOLOGICAL ISSUES** — v7/v8 biologically defensible; bio clean round 2 of 2 → BIO CONVERGED" |
| Eng (R8) | "**NO NEW ENGINEERING ISSUES** — v7/v8 engineering-locked; final confirmation of CONVERGENCE" |

Per user-set stop criterion ("until 2 consecutive rounds yield no new conclusions"):
- Eng: R6 (0) + R7 (0) + R8 (0) = 3 clean rounds (criterion met after R7; R8 confirms)
- Bio: R7 (0) + R8 (0) = 2 clean rounds (criterion met)

**Both tracks converged. v7 is locked.**

---

## 10. Implementation status

This is a **DESIGN PROPOSAL** — NOT a built tool. Implementation NOT YET ATTEMPTED.

**Pre-implementation requirement (per F3-10)**: D12 feasibility prototype on TAIR10 to validate +3pp recall target is achievable before committing to the full v7 build.

**Estimated v7 implementation cost**: 6-12 person-months (Nextflow DSL2 + 11 class cascades + 4-track GBT training + 6-genome calibration + integration tests + CI infrastructure).

**Decision pending**: whether to proceed with v7 implementation or defer pending other priorities.

---

## 11. Files produced during iteration

The 7 versioned proposal documents are at `/tmp/refiner_proposal_v{1..7}.md`. v1 (the initial 5-change R5/R6/R7/R8 proposal that this iteration started from) is reconstructable from session history; v2-v7 saved as artifacts.

Final document (this file): `/home/shuoc/tool/mdl-repeat/REFINER_MDL_V7_CONVERGED_PROPOSAL.md`.

Empirical baseline benchmark: `/home/shuoc/tool/mdl-repeat/REFINER_MDL_BENCHMARK_RESULT.md`.
