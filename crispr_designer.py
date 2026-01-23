# This project was made by Maheen Abbasi on Jan 2026.

# ------------------------------------------------------------------------------------ #
# Basic File and DNA setup
# ------------------------------------------------------------------------------------ #

# Step 1.1: Accept DNA from text
def get_dna_sequence():
    seq = input("Enter DNA sequence: ").upper().strip()
    return seq

# Step 1:2 Accept DNA from FASTA file by normalizing
def read_fasta(filename):
    sequence = ""

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            sequence += line.upper()
    
    return sequence

# Step 2: Validate that all bases in the DNA sequence are one of: A, T, C, G
def validate_dna(seq):
    valid_bases = {"A", "T", "C", "G"}

    for base in seq:
        if base not in valid_bases:
            raise ValueError("Invalid DNA sequence: only A, T, C, G are allowed. Please check your sequence.")
    
# Step 3: Get the reverse complement of each base
def reverse_complement(seq):
    complement = {
        "A":"T",
        "T":"A",
        "C":"G",
        "G":"C"
    }

    rev_comp = ""
    for base in reversed(seq):
        rev_comp += complement[base]

    return rev_comp

# Step 4: Find PAM sites (NGG, where N is any base)
def find_grnas(seq):
    guides = []

    # This is a sliding window scan; move the DNA one base at a time. For each base:
    #   --> look at the first 20 bases for a possible guide
    #   --> look at the next 3 bases for a possible PAM (check if PAM = NGG)

    # Note: CRISPR needs a 20bp guide + a 3bp PAM, so the sequence needs to be >= 23bp
    for i in range(len(seq) - 23 + 1):
        guide = seq[i:i+20]
        pam = seq[i+20:i+23]

        # guides are stored in a dictionary
        if pam[1:] == "GG":
            guides.append({
                "guide": guide,
                "pam": pam,
                "position": i,
                "strand": "+"
            })
    
    return guides

# ------------------------------------------------------------------------------------ #
# Methods for Scoring and Ranking gRNA
#   --> Rule based CRISPR Scoring heuristic
# ------------------------------------------------------------------------------------ #

# GC Content Calculation
def gc_content(seq):
    gc = (seq.count("G")) + (seq.count("C"))
    return gc/len(seq)

def has_bad_repeats(seq):
    bad_patterns = ["AAAA", "TTTT", "CCCC", "GGGG"]
    return any(p in seq for p in bad_patterns)

def position_penalty(guide):
    # Avoid having a T at position 1 (it is a known weak transcription start)
    if guide[0] == "T":
        return 10
    return 0

def score_grna(guide):
    # By default, score starts at 100
    score = 100
    gc = gc_content(guide)

    # Weighted scoring for gc content
    score -= abs(gc - 0.5) * 40

    # pos penalty
    score -= (position_penalty(guide))

    # bad repeat penalty
    if has_bad_repeats(guide):
        score -= 30
    
    # off-target penalty
    if dna:
        off_targets = count_off_targets(guide, dna)
        score -= off_targets * 10

    return round(score, 2)

def rank_grnas(guides, dna):
    for g in guides:
        # add entries to the guides dictionary for gc content and score
        g["gc"] = round(gc_content(g["guide"]), 2)
        g["off_targets"] = count_off_targets(g["guide"], dna)
        g["score"] = score_grna(g["guide"])
    
    # Sort the dictionary by score
    return sorted(guides, key=lambda x: x["score"], reverse=True)

# ------------------------------------------------------------------------------------ #
# Off-target checking
# ------------------------------------------------------------------------------------ #

# Step 1: Compute Hamming distance
def hamming_distance(a, b):
    return sum(x != y for x, y in zip(a,b))

def count_off_targets(guide, dna, max_mismatches=2):
    count = 0

    for i in range(len(dna) - len(guide)):
        window = dna[i:i+len(guide)]
        if hamming_distance(guide, window) <= max_mismatches:
            count += 1
    
    # When returnign, subtract the perfect match, aka the guide we are referring to
    return count - 1

# ------------------------------------------------------------------------------------ #
# Main
# ------------------------------------------------------------------------------------ #
if __name__ == "__main__":
    # let the user choose how to give input
    mode = input("Input type (raw / fasta): ").lower()

    if mode == "fasta":
        filename = input("FASTA filename: ")
        dna = read_fasta(filename)
        validate_dna(dna)
    else:
        dna = get_dna_sequence()
        
    validate_dna(dna)

    guides = find_grnas(dna)
    ranked = rank_grnas(guides, dna)

    print(f"Total candidate guides: {len(ranked)}\n")

    for g in ranked[:5]:
        print(
            f"Guide: {g['guide']} | "
            f"PAM: {g['pam']} | "
            f"Off-targets: {g['off_targets']} | "
            f"GC: {g['gc']} | "
            f"Score: {g['score']} | "
            f"Pos: {g['position']} | "
            f"Strand: {g['strand']}"
        )