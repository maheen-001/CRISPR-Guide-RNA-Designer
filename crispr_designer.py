# This project was made by Maheen Abbasi on Jan. 22, 2026.

# Step 1: Accept DNA from user + ensure it is valid
def get_dna_sequence():
    seq = input("Enter DNA sequence: ").upper().strip()

    valid_bases = {"A", "T", "C", "G"}
    for base in seq:
        if base not in valid_bases:
            raise ValueError("Invalid DNA sequence: only A, T, C, G are allowed. Please check your sequence.")
    
    return seq

# Step 2: Get the reverse complement of each base
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

# Step 3: Find PAM sites (NGG, where N is any base)
def find_grnas(seq):
    guides = []

    # This is a sliding window scan; move the DNA one base at a time. For each base:
    #   --> look at the first 20 bases for a possible guide
    #   --> look at the next 3 bases for a possible PAM (check if PAM = NGG)

    # Note: CRISPR needs a 20bp guide + a 3bp PAM, so the sequence needs to be >= 23bp
    for i in range(len(seq) - 23 + 1):
        guide = seq[i:i+20]
        pam = seq[i+20:i+23]

        if pam[1:] == "GG":
            guides.append({
                "guide": guide,
                "pam": pam,
                "position": i,
                "strand": "+"
            })
    
    return guides

# Step 4: GC Content Calculation
def gc_content(seq):
    gc = (seq.count("G")) + (seq.count("C"))
    return gc/len(seq)

# Main
if __name__ == "__main__":
    dna = get_dna_sequence()
    # print("Sequence length: ", len(dna))
    # guides = find_grnas(dna)
    # print("Found guides:", len(guides))
    print(gc_content(dna))