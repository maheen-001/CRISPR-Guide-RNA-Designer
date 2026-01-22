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

# Main
if __name__ == "__main__":
    dna = get_dna_sequence()
    print("Sequence length: ", len(dna))