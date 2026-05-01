#!/usr/bin/env python3
"""
Diagnostic synthetic genome builder for mdl-repeat failure analysis.
Implements Genomes A, B, C to test hypotheses about missed repeat families.

Usage: python3 tools/diag_setup.py [--step {all,cluster,genomes,check}]
"""

import os
import sys
import re
import random
import argparse
import struct

OUTDIR = "/tmp/ath_bench/diag"
MISSED_FA = "/tmp/ath_bench/missed_long.fa"
CLUSTER_REPS_FA = "/tmp/ath_bench/missed_cluster_reps.fa"
CLUSTER_TSV = "/tmp/ath_bench/missed_clusters.tsv"
CHR4_FA = "/tmp/ath_bench/chr4.fa"
RANDOM_SEED = 12345

# ── helpers ──────────────────────────────────────────────────────────────────

def load_fasta(path):
    """Return OrderedDict {name: seq} preserving order."""
    from collections import OrderedDict
    seqs = OrderedDict()
    curr = None
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                curr = line[1:].split()[0]
                seqs[curr] = []
            else:
                if curr is not None:
                    seqs[curr].append(line.upper())
    return {k: "".join(v) for k, v in seqs.items()}


def write_fasta(path, seqs, line_width=80):
    """seqs: list of (name, seq) tuples."""
    with open(path, "w") as f:
        for name, seq in seqs:
            f.write(f">{name}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i+line_width] + "\n")


def write_bed(path, records):
    """records: list of (chrom, start, end, name, score, strand)."""
    with open(path, "w") as f:
        for r in records:
            f.write("\t".join(str(x) for x in r) + "\n")


def gc_content(seq):
    s = seq.upper()
    gc = s.count("G") + s.count("C")
    return gc / max(len(s), 1)


def is_low_complexity(seq):
    """Return True if the sequence is likely low-complexity."""
    s = seq.upper()
    gc = gc_content(s)
    if gc < 0.25 or gc > 0.75:
        return True
    # Check max poly-A/T/G/C run > 15 bp
    for base in "ACGT":
        runs = re.findall(f"{base}{{15,}}", s)
        if runs:
            return True
    # Check if dominant 4-mer > 5x expected
    n = len(s)
    expected = n / 256
    tetramers = {}
    for i in range(n - 3):
        t = s[i:i+4]
        tetramers[t] = tetramers.get(t, 0) + 1
    if tetramers and max(tetramers.values()) > expected * 5:
        return True
    return False


def mutate_seq(seq, rng, divergence):
    """Introduce substitutions at given divergence rate."""
    bases = list(seq.upper())
    n_mut = int(len(bases) * divergence)
    positions = rng.sample(range(len(bases)), min(n_mut, len(bases)))
    for pos in positions:
        orig = bases[pos]
        alts = [b for b in "ACGT" if b != orig]
        bases[pos] = rng.choice(alts)
    return "".join(bases)


# ── chr4 background model ─────────────────────────────────────────────────────

def build_dinucleotide_model(fa_path, sample_bp=2_000_000):
    """
    Train a 1st-order Markov model on chr4.
    Returns (initial_freq dict, transition_matrix dict-of-dict).
    """
    print(f"  Loading chr4 for background model (first {sample_bp} bp) ...")
    seqs = load_fasta(fa_path)
    chrom = list(seqs.values())[0]

    # Skip leading Ns, take first sample_bp ACGT bp
    clean = re.sub(r"[^ACGT]", "", chrom.upper())
    chunk = clean[:sample_bp]
    if len(chunk) < 100_000:
        raise RuntimeError(f"Not enough non-N chr4 sequence ({len(chunk)} bp)")

    bases = "ACGT"
    # Count transitions
    trans = {b: {b2: 1 for b2 in bases} for b in bases}  # pseudo-count 1
    init = {b: 1 for b in bases}

    for i, c in enumerate(chunk):
        if c not in bases:
            continue
        init[c] = init.get(c, 0) + 1
        if i > 0:
            prev = chunk[i-1]
            if prev in bases:
                trans[prev][c] = trans[prev].get(c, 0) + 1

    # Normalize
    total = sum(init.values())
    init_p = {b: init[b] / total for b in bases}

    trans_p = {}
    for b in bases:
        row_total = sum(trans[b].values())
        trans_p[b] = {b2: trans[b][b2] / row_total for b2 in bases}

    gc = init_p["G"] + init_p["C"]
    print(f"  Background model: GC={gc:.3f}, trained on {len(chunk)} bp")
    return init_p, trans_p


def generate_markov_seq(length, init_p, trans_p, rng):
    """Generate a DNA sequence using a 1st-order Markov model."""
    bases = "ACGT"
    # Sample initial base
    r = rng.random()
    cumulative = 0.0
    cur = bases[0]
    for b in bases:
        cumulative += init_p[b]
        if r <= cumulative:
            cur = b
            break

    result = [cur]
    for _ in range(length - 1):
        r = rng.random()
        cumulative = 0.0
        row = trans_p[cur]
        for b in bases:
            cumulative += row[b]
            if r <= cumulative:
                cur = b
                break
        result.append(cur)
    return "".join(result)


# ── representative selection ──────────────────────────────────────────────────

CLUSTER_INFO = [
    # (cluster_id, n_members, rep_len, gc_pct)
    (0,  7, 1425, 24.1),  # LOW COMPLEXITY – AT-rich; SKIP
    (1,  7, 1700, 38.3),  # ~medium-long, good GC
    (2,  6, 2194, 37.1),  # long, good GC
    (3,  6, 22266, 35.8), # giant (22 kb) – too large for diagnostic; SKIP
    (4,  6, 1547, 30.2),  # medium-long, acceptable GC
    (5,  5, 1611, 43.8),  # medium-long, good GC
    (6,  5, 1556, 40.0),  # medium, good GC
    (7,  5, 959,  25.1),  # AT-rich; SKIP
    (8,  4, 703,  43.8),  # short-medium, good GC
    (9,  4, 2283, 45.6),  # long, good GC
    (10, 4, 752,  41.4),  # short-medium, good GC
    (11, 4, 1373, 35.0),  # medium, acceptable GC
    (12, 4, 790,  40.8),  # short-medium, good GC
    (13, 4, 4783, 33.5),  # very long; skip for size
    (14, 4, 662,  32.3),  # short, acceptable GC
    (15, 3, 580,  35.5),  # short, good GC
    (16, 3, 4915, 33.0),  # very long; skip for size
    (17, 3, 1620, 24.8),  # AT-rich; SKIP
]

# Picked 6 reps: 2 short (~600-800 bp), 2 medium (~1500-1700 bp), 2 long (~2000-2300 bp)
# Skip cluster_0 (AT-rich), cluster_7 (AT-rich), cluster_17 (AT-rich), cluster_3 (22kb), cluster_13/16 (>4kb)
PICKED = [
    # (cluster_id, role, length_class)
    # cluster_15 (580bp): LC due to poly-AAAA runs (7.9x 4-mer) → REPLACED by cluster_8 (703bp, OK)
    # cluster_1 (1700bp): LC due to polyA=16 run → REPLACED by cluster_11 (1373bp, OK)
    # cluster_5 (1611bp): LC due to TTTT 5.7x → keep cluster_6 (1556bp) which tested OK
    # cluster_14 (662bp): borderline LC (TTTT 5.0x) → use cluster_10 (752bp) instead
    (8,  "short_A",  "short"),   # n=4, 703 bp, gc=43.8%, OK
    (10, "short_B",  "short"),   # n=4, 752 bp, gc=41.4%, OK
    (11, "medium_A", "medium"),  # n=4, 1373 bp, gc=35.0%, OK  (replaces cluster_1)
    (6,  "medium_B", "medium"),  # n=5, 1556 bp, gc=40.0%, OK
    (2,  "long_A",   "long"),    # n=6, 2194 bp, gc=37.1%, OK
    (9,  "long_B",   "long"),    # n=4, 2283 bp, gc=45.6%, OK
]

# For Genome C: use long_A (2194 bp) as LINE-like, short_A (703 bp) as LTR-like


def step_cluster(rng):
    """Extract picked representatives and write /tmp/ath_bench/diag/representatives.fa"""
    print("\n=== Step 1: Extract picked representatives ===")
    os.makedirs(OUTDIR, exist_ok=True)

    cluster_seqs = load_fasta(CLUSTER_REPS_FA)
    # Map cluster_id -> seq name
    name_map = {}
    for name in cluster_seqs:
        m = re.match(r"cluster_(\d+)_", name)
        if m:
            name_map[int(m.group(1))] = name

    print(f"  Loaded {len(cluster_seqs)} cluster reps from {CLUSTER_REPS_FA}")
    print(f"  Picking {len(PICKED)} representatives:")

    reps = []
    for cid, role, lclass in PICKED:
        name = name_map.get(cid)
        if name is None:
            print(f"  ERROR: cluster_{cid} not found in {CLUSTER_REPS_FA}")
            sys.exit(1)
        seq = cluster_seqs[name]
        gc = gc_content(seq)
        lc = is_low_complexity(seq)
        status = "LC!" if lc else "OK"
        print(f"  cluster_{cid:2d}  role={role:12s}  len={len(seq):5d}  gc={gc:.3f}  {status}")
        # New informative header
        new_name = f"rep_c{cid}_{role}_len{len(seq)}"
        reps.append((new_name, seq))

    out_path = os.path.join(OUTDIR, "representatives.fa")
    write_fasta(out_path, reps)
    print(f"  Written {len(reps)} reps to {out_path}")
    return {role: (seq, cid) for (cid, role, _), (_, seq) in zip(PICKED, reps)}


# ── insertion helper ──────────────────────────────────────────────────────────

def find_free_pos(used_ranges, length, genome_len, rng, tries=2000):
    """Find a non-overlapping insertion site with 500 bp flanking buffer."""
    buf = 500
    for _ in range(tries):
        pos = rng.randint(buf, genome_len - length - buf - 1)
        ok = True
        for (s, e) in used_ranges:
            if pos < e + 100 and pos + length > s - 100:
                ok = False
                break
        if ok:
            used_ranges.append((pos, pos + length))
            return pos
    return None


# ── Genome A ──────────────────────────────────────────────────────────────────

def build_genome_A(reps, init_p, trans_p, rng):
    """
    Genome A: each picked representative ×8 copies, 5–15% divergence,
    in clean Markov background.
    """
    print("\n=== Step 2a: Building Genome A ===")
    GENOME_LEN = 5_000_000
    N_COPIES = 8
    DIV_MIN, DIV_MAX = 0.05, 0.15

    genome = list(generate_markov_seq(GENOME_LEN, init_p, trans_p, rng))
    used = []
    truth = []
    chrom = "genomA"

    for role, (cons, cid) in reps.items():
        for copy_i in range(N_COPIES):
            div = rng.uniform(DIV_MIN, DIV_MAX)
            copy_seq = mutate_seq(cons, rng, div)
            pos = find_free_pos(used, len(copy_seq), GENOME_LEN, rng)
            if pos is None:
                print(f"  WARNING: could not place copy {copy_i} of {role}")
                continue
            genome[pos:pos+len(copy_seq)] = list(copy_seq)
            truth.append((chrom, pos, pos+len(copy_seq),
                          f"rep_c{cid}_{role}", 0, "+"))

    genome_str = "".join(genome)
    n_inserted = len(truth)
    total_rep_bp = sum(r[2]-r[1] for r in truth)
    print(f"  Genome A: {GENOME_LEN:,} bp, {n_inserted} inserted copies")
    print(f"  Repeat content: {total_rep_bp/GENOME_LEN*100:.1f}%")

    fa_path  = os.path.join(OUTDIR, "genomA.fa")
    bed_path = os.path.join(OUTDIR, "genomA_truth.bed")
    write_fasta(fa_path, [(chrom, genome_str)])
    write_bed(bed_path, sorted(truth, key=lambda r: r[1]))
    print(f"  Written: {fa_path}, {bed_path}")
    return genome_str, truth


# ── Genome B ──────────────────────────────────────────────────────────────────

def build_genome_B(reps, init_p, trans_p, rng):
    """
    Genome B: dominant family (×100) + victim families (×8 each).
    The dominant family has a 14-bp segment copied from victim's consensus
    engineered into it (simulating cross-talk).
    """
    print("\n=== Step 2b: Building Genome B ===")
    GENOME_LEN = 5_000_000
    N_DOMINANT = 100
    N_VICTIM_COPIES = 8
    DIV_DOM = 0.05         # dominant: tight family
    DIV_VICTIM_MIN = 0.05
    DIV_VICTIM_MAX = 0.15

    # Pick dominant: medium_A (1700 bp, n=7) — largest cluster, prominent signal
    dom_role = "medium_A"
    dom_cons, dom_cid = reps[dom_role]

    # Victim: use short_A (580 bp) — shortest, most likely to get shadowed
    victim_role = "short_A"
    victim_cons, victim_cid = reps[victim_role]

    # ── Engineer cross-talk: copy 14-bp segment from victim pos 200..214
    #    into dominant consensus at position 300..314
    shared_14mer = victim_cons[200:214]
    dom_cons_modified = (dom_cons[:300] + shared_14mer + dom_cons[314:])
    print(f"  Shared 14-mer (victim pos 200-214): {shared_14mer}")
    print(f"  Inserted into dominant at pos 300-314")
    print(f"  Dominant: cluster_{dom_cid} ({dom_role}), len={len(dom_cons_modified)}")
    print(f"  Victim:   cluster_{victim_cid} ({victim_role}), len={len(victim_cons)}")

    genome = list(generate_markov_seq(GENOME_LEN, init_p, trans_p, rng))
    used = []
    truth = []
    chrom = "genomB"

    # Insert dominant ×100
    for i in range(N_DOMINANT):
        div = rng.uniform(0.0, DIV_DOM)
        copy_seq = mutate_seq(dom_cons_modified, rng, div)
        pos = find_free_pos(used, len(copy_seq), GENOME_LEN, rng)
        if pos is None:
            print(f"  WARNING: could not place dominant copy {i}")
            continue
        genome[pos:pos+len(copy_seq)] = list(copy_seq)
        truth.append((chrom, pos, pos+len(copy_seq),
                      f"dominant_c{dom_cid}_{dom_role}", 0, "+"))

    # Insert all 6 picked families ×8 each (including victim)
    for role, (cons, cid) in reps.items():
        for copy_i in range(N_VICTIM_COPIES):
            div = rng.uniform(DIV_VICTIM_MIN, DIV_VICTIM_MAX)
            copy_seq = mutate_seq(cons, rng, div)
            pos = find_free_pos(used, len(copy_seq), GENOME_LEN, rng)
            if pos is None:
                print(f"  WARNING: could not place victim copy {copy_i} of {role}")
                continue
            tag = "victim" if role == victim_role else "background"
            genome[pos:pos+len(copy_seq)] = list(copy_seq)
            truth.append((chrom, pos, pos+len(copy_seq),
                          f"{tag}_c{cid}_{role}", 0, "+"))

    genome_str = "".join(genome)
    n_inserted = len(truth)
    total_rep_bp = sum(r[2]-r[1] for r in truth)
    print(f"  Genome B: {GENOME_LEN:,} bp, {n_inserted} inserted copies")
    print(f"  Repeat content: {total_rep_bp/GENOME_LEN*100:.1f}%")
    print(f"  Dominant copies: {sum(1 for t in truth if 'dominant' in t[3])}")
    print(f"  Victim copies:   {sum(1 for t in truth if 'victim' in t[3])}")

    fa_path  = os.path.join(OUTDIR, "genomB.fa")
    bed_path = os.path.join(OUTDIR, "genomB_truth.bed")
    write_fasta(fa_path, [(chrom, genome_str)])
    write_bed(bed_path, sorted(truth, key=lambda r: r[1]))
    print(f"  Written: {fa_path}, {bed_path}")
    return genome_str, truth


# ── Genome C ──────────────────────────────────────────────────────────────────

def build_genome_C(reps, init_p, trans_p, rng):
    """
    Genome C: nested elements.
    - LINE-like (long_A, ~2194 bp) ×10 copies
    - LTR-like (short_A, ~580 bp) nested inside 3 of the LINE copies
    - LTR-like standalone ×5 copies
    """
    print("\n=== Step 2c: Building Genome C ===")
    GENOME_LEN = 5_000_000
    N_LINE = 10          # total LINE copies
    N_LINE_WITH_LTR = 3  # how many LINEs contain a nested LTR
    N_STANDALONE_LTR = 5 # standalone LTR copies
    DIV_LINE = 0.08
    DIV_LTR = 0.05

    line_role = "long_A"
    ltr_role  = "short_A"
    line_cons, line_cid = reps[line_role]
    ltr_cons,  ltr_cid  = reps[ltr_role]

    # LTR insertion offset inside LINE: aim for ~30% into element
    ltr_insert_offset = int(len(line_cons) * 0.30)
    print(f"  LINE: cluster_{line_cid} ({line_role}), len={len(line_cons)}")
    print(f"  LTR:  cluster_{ltr_cid} ({ltr_role}), len={len(ltr_cons)}")
    print(f"  LTR nested at offset {ltr_insert_offset} within LINE")

    # Build nested LINE consensus (LINE with LTR inserted)
    nested_line_cons = (line_cons[:ltr_insert_offset] +
                        ltr_cons +
                        line_cons[ltr_insert_offset:])

    genome = list(generate_markov_seq(GENOME_LEN, init_p, trans_p, rng))
    used = []
    truth = []
    chrom = "genomC"

    # Insert LINE copies — first N_LINE_WITH_LTR get nested LTR, rest are plain
    for i in range(N_LINE):
        if i < N_LINE_WITH_LTR:
            # nested LINE: length = len(line_cons) + len(ltr_cons)
            template = nested_line_cons
            line_len = len(template)
            div = rng.uniform(0.0, DIV_LINE)
            copy_seq = mutate_seq(template, rng, div)
            pos = find_free_pos(used, line_len, GENOME_LEN, rng)
            if pos is None:
                print(f"  WARNING: could not place nested LINE copy {i}")
                continue
            genome[pos:pos+line_len] = list(copy_seq)
            # BED for the whole LINE
            truth.append((chrom, pos, pos + line_len,
                          f"LINE_nested_c{line_cid}", 0, "+"))
            # BED for the embedded LTR
            ltr_start = pos + ltr_insert_offset
            ltr_end   = ltr_start + len(ltr_cons)
            truth.append((chrom, ltr_start, ltr_end,
                          f"LTR_nested_c{ltr_cid}", 0, "+"))
        else:
            # plain LINE
            template = line_cons
            div = rng.uniform(0.0, DIV_LINE)
            copy_seq = mutate_seq(template, rng, div)
            pos = find_free_pos(used, len(copy_seq), GENOME_LEN, rng)
            if pos is None:
                print(f"  WARNING: could not place plain LINE copy {i}")
                continue
            genome[pos:pos+len(copy_seq)] = list(copy_seq)
            truth.append((chrom, pos, pos + len(copy_seq),
                          f"LINE_plain_c{line_cid}", 0, "+"))

    # Insert standalone LTR copies
    for i in range(N_STANDALONE_LTR):
        div = rng.uniform(0.0, DIV_LTR)
        copy_seq = mutate_seq(ltr_cons, rng, div)
        pos = find_free_pos(used, len(copy_seq), GENOME_LEN, rng)
        if pos is None:
            print(f"  WARNING: could not place standalone LTR copy {i}")
            continue
        genome[pos:pos+len(copy_seq)] = list(copy_seq)
        truth.append((chrom, pos, pos + len(copy_seq),
                      f"LTR_solo_c{ltr_cid}", 0, "+"))

    genome_str = "".join(genome)
    n_inserted = len(truth)
    total_rep_bp = sum(r[2]-r[1] for r in truth)
    print(f"  Genome C: {GENOME_LEN:,} bp, {n_inserted} BED records")
    n_line = sum(1 for t in truth if "LINE" in t[3])
    n_ltr  = sum(1 for t in truth if "LTR" in t[3])
    print(f"  LINE records: {n_line}, LTR records: {n_ltr}")
    print(f"    LTR nested: {sum(1 for t in truth if 'LTR_nested' in t[3])}")
    print(f"    LTR solo:   {sum(1 for t in truth if 'LTR_solo' in t[3])}")
    print(f"  Repeat content: {total_rep_bp/GENOME_LEN*100:.1f}%")

    fa_path  = os.path.join(OUTDIR, "genomC.fa")
    bed_path = os.path.join(OUTDIR, "genomC_truth.bed")
    write_fasta(fa_path, [(chrom, genome_str)])
    write_bed(bed_path, sorted(truth, key=lambda r: r[1]))
    print(f"  Written: {fa_path}, {bed_path}")
    return genome_str, truth


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--step", choices=["all", "cluster", "genomes", "check"],
                        default="all")
    args = parser.parse_args()

    os.makedirs(OUTDIR, exist_ok=True)
    rng = random.Random(RANDOM_SEED)

    print(f"Output directory: {OUTDIR}")
    print(f"Random seed: {RANDOM_SEED}")

    # Step 1: extract representatives
    reps = step_cluster(rng)

    if args.step == "cluster":
        return

    # Build background model
    print("\n=== Building chr4 background model ===")
    init_p, trans_p = build_dinucleotide_model(CHR4_FA)

    # Step 2: build genomes
    build_genome_A(reps, init_p, trans_p, rng)
    build_genome_B(reps, init_p, trans_p, rng)
    build_genome_C(reps, init_p, trans_p, rng)

    print("\n=== All genomes built successfully ===")
    print(f"Output files in {OUTDIR}:")
    for f in sorted(os.listdir(OUTDIR)):
        p = os.path.join(OUTDIR, f)
        if os.path.isfile(p):
            sz = os.path.getsize(p)
            print(f"  {f}  ({sz:,} bytes)")


if __name__ == "__main__":
    main()
