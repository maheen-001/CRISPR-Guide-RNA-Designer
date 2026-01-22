# This project was made by Maheen Abbasi on Jan. 22, 2026.

# Step 1: Accept DNA from user + ensure it is valid
def get_dna_sequence():
    seq = input("Enter DNA sequence: ").upper().strip()

    valid_bases = {"A", "T", "C", "G"}
    for base in seq:
        if base not in valid_bases:
            raise ValueError("Invalid DNA sequence: only A, T, C, G are allowed. Please check your sequence.")
    
    return seq

# Main
if __name__ == "__main__":
    dna = get_dna_sequence()
    print("Sequence length: ", len(dna))