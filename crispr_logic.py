# This project was made by Maheen Abbasi on Jan 2026.

CAS_ENZYMES = {
    "SpCas9": {
        "guide_length": 20,
        "pam_patterns": ["NGG", "NAG"],
        "pam_length": 3,
        "pam_side": "downstream",
        "relaxed_penalty": 15
    },
    "SaCas9": {
        "guide_length": 21,
        "pam_patterns": ["NNGRRT"],
        "pam_length": 6,
        "pam_side": "downstream",
        "relaxed_penalty": 0
    },
    "Cas12a": {
        "guide_length": 23,
        "pam_patterns": ["TTTV"],
        "pam_length": 4,
        "pam_side": "upstream",
        "relaxed_penalty": 0
    }
}

# ------------------------------------------------------------------------------------ #
# Basic File and DNA setup
# ------------------------------------------------------------------------------------ #

# Step 1.1: Accept DNA from text
def get_dna_sequence():
    seq = input("Enter DNA sequence: ").upper().strip()
    return seq

# Step 1.2: Accept DNA from FASTA file by normalizing
def read_fasta(filename):
    sequence = ""

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            sequence += line.upper()
    
    return sequence

# Step 1.3: Normalize DNA data
def normalize(seq):
    lines = seq.splitlines()
    cleaned = ""

    for line in lines:
        if line.startswith(">"):
            continue
        cleaned += line.strip()

    return cleaned.upper()

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

# Step 4: Check for valid PAMs
def pam_matches(pam, pattern):
    if len(pam) != len(pattern):
        return False
    
    # Compare PAM sequence against a PAM pattern (NGG, NAG, NNGRRT)
    iupac = {
        "N": {"A", "T", "C", "G"},
        "R": {"A", "G"},
        "Y": {"C", "T"},
        "V": {"A", "C", "G"}
    }

    for base, rule in zip(pam, pattern):
        if rule in iupac:
            if base not in iupac[rule]:
                return False
        else:
            if base != rule:
                return False
    
    return True

def classify_pam(pam, enzyme = "SpCas9"):
    patterns = CAS_ENZYMES[enzyme]["pam_patterns"]

    for pattern in patterns:
        if pam_matches(pam, pattern):
            if enzyme == "SpCas9" and pattern != "NGG":
                return "non-canonical"
            return "canonical"
    return None

    
# Step 5: Find PAM sites
def find_grnas(seq, enzyme = "SpCas9"):
    guides = []
    cfg = CAS_ENZYMES[enzyme]

    g_len = cfg["guide_length"]
    pam_len = cfg["pam_length"]
    pam_side = cfg["pam_side"]

    total_len = g_len + pam_len

    # Slide a window along the DNA, and at each position, check:
    #   -> if a guide + PAM started here (or vice versa), would it be valid for THIS enzyme
    for i in range(len(seq) - total_len + 1):

        if pam_side == "downstream":
            # downstream PAM (Cas12a)
            guide = seq[i:(i + g_len)]
            pam = seq[i+g_len: (i + g_len + pam_len)]
        else:
            # upstream PAM (Cas12a)
            pam = seq[i:(i + pam_len)]
            guide = seq[i+pam_len: (i + pam_len + g_len)]

        if len(guide) != g_len or len(pam) != pam_len:
            continue

        pam_class = classify_pam(pam, enzyme)

        # guides are stored in a dictionary
        if pam_class:
            guides.append({
                "guide": guide,
                "pam": pam,
                "pam_class": pam_class,
                "position": i,
                "strand": "+",
                "enzyme": enzyme
            })
    
    return guides

def find_all_grnas(seq, enzyme = "SpCas9"):
    # Forward strand
    forward = find_grnas(seq, enzyme)

    # Reverse
    rev_seq = reverse_complement(seq)
    reverse = find_grnas(rev_seq, enzyme)

    # Update strand info for reverse hits since find_grnas default to +
    for g in reverse:
        g["strand"] = "-"

    return forward + reverse

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

def score_grna(guide, dna, pam_class, enzyme = "SpCas9"):
    # By default, score starts at 100
    score = 100
    gc = gc_content(guide)

    # Weighted scoring for gc content with guide length normalization
    score -= abs(gc - 0.5) * (40 * (20/len(guide)))

    # pos penalty
    score -= (position_penalty(guide))

    # bad repeat penalty
    if has_bad_repeats(guide):
        score -= 30
    
    # non-canonical/relaxed PAM penalty
    if pam_class == "non-canonical":
        score -= CAS_ENZYMES[enzyme]["relaxed_penalty"]

    # off-target penalty
    off_targets = count_off_targets(guide, dna)
    score -= off_targets * 10

    return round(score, 2)

def rank_grnas(guides, dna):
    for g in guides:
        # add entries to the guides dictionary for gc content and score
        g["gc"] = round(gc_content(g["guide"]), 2)
        g["off_targets"] = count_off_targets(g["guide"], dna)
        g["score"] = score_grna(
            g["guide"], 
            dna, 
            g["pam_class"],
            g["enzyme"]
            )
    
    # Sort the dictionary by score
    return sorted(guides, key=lambda x: x["score"], reverse=True)

# ------------------------------------------------------------------------------------ #
# Off-target checking
# ------------------------------------------------------------------------------------ #

# Step 1: Compute Hamming distance
def hamming_distance(a, b):
    return sum(x != y for x, y in zip(a,b))

# Step 2: Check for off-targets
def count_off_targets(guide, dna, max_mismatches=2):
    count = 0
    rev_dna = reverse_complement(dna)

    for genome in (dna, rev_dna):
        for i in range(len(genome) - len(guide) + 1):
            window = genome[i:i+len(guide)]
            if hamming_distance(guide, window) <= max_mismatches:
                count += 1
    
    # When returnign, subtract the perfect match, aka the guide we are referring to
    return count - 1