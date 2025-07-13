import streamlit as st
import datetime

# Genetic code dictionary
genetic_code = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
    'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

def transcribe_dna_to_rna(dna_sequence):
    return dna_sequence.replace("T", "U")

def translate_full_sequence(rna_sequence):
    protein = []
    for i in range(0, len(rna_sequence) - 2, 3):
        codon = rna_sequence[i:i+3]
        protein.append(genetic_code.get(codon, 'X'))
    return ''.join(protein)

def find_all_proteins(rna_sequence):
    proteins = []
    i = 0
    while i < len(rna_sequence) - 2:
        codon = rna_sequence[i:i+3]
        if codon == "AUG":
            protein = []
            for j in range(i, len(rna_sequence) - 2, 3):
                codon_j = rna_sequence[j:j+3]
                aa = genetic_code.get(codon_j, "X")
                if aa == "*":
                    proteins.append("".join(protein))
                    break
                protein.append(aa)
            i += 3
        else:
            i += 1
    return proteins if proteins else ["No proteins found"]

def validate_dna_sequence(seq):
    invalid = set(seq) - {"A", "T", "C", "G"}
    if invalid:
        raise ValueError(f"Invalid DNA characters: {invalid}")

def validate_rna_sequence(seq):
    invalid = set(seq) - {"A", "U", "C", "G"}
    if invalid:
        raise ValueError(f"Invalid RNA characters: {invalid}")

def reverse_complement(seq):
    complement = str.maketrans("ATCG", "TAGC")
    return seq.translate(complement)[::-1]

def get_reading_frames(seq):
    frames = [seq[i:] for i in range(3)]
    rev = reverse_complement(seq)
    frames += [rev[i:] for i in range(3)]
    return [transcribe_dna_to_rna(f) for f in frames], rev

def get_protein_positions(rna_sequence):
    positions = []
    i = 0
    while i < len(rna_sequence) - 2:
        codon = rna_sequence[i:i+3]
        if codon == "AUG":
            start = i
            for j in range(i, len(rna_sequence) - 2, 3):
                codon_j = rna_sequence[j:j+3]
                if genetic_code.get(codon_j) == "*":
                    end = j + 3
                    positions.append((start, end))
                    break
            i += 3
        else:
            i += 1
    return positions

def process_sequence(seq, seq_type):
    output = ""
    if seq_type == 'DNA':
        validate_dna_sequence(seq)
        rna_seq = transcribe_dna_to_rna(seq)
        full_translation = translate_full_sequence(rna_seq)
        output += f"\n🔬 Full Translation (without start codon requirement):\n{full_translation}\n"

        rna_frames, rev_comp = get_reading_frames(seq)
        output += f"\n🔁 Reverse Complement:\n{rev_comp}\n"

        for i, rna in enumerate(rna_frames, 1):
            direction = "Forward" if i <= 3 else "Reverse"
            output += f"\n🧬 Reading Frame {i} ({direction}):\n{rna}\n"
            proteins = find_all_proteins(rna)
            for idx, prot in enumerate(proteins, 1):
                output += f"   🔹 Protein {idx}: {prot}\n"

    elif seq_type == 'RNA':
        validate_rna_sequence(seq)
        full_translation = translate_full_sequence(seq)
        output += f"\n🔬 Full Translation (without start codon requirement):\n{full_translation}\n"

        # Treat RNA like DNA to extract reading frames and reverse complement
        rna_frames = [seq[i:] for i in range(3)]
        rev = seq[::-1]  # simple reverse for RNA
        rna_frames += [rev[i:] for i in range(3)]

        output += f"\n🔁 Reverse (mirrored) RNA:\n{rev}\n"

        for i, rna in enumerate(rna_frames, 1):
            direction = "Forward" if i <= 3 else "Reverse"
            output += f"\n🧬 Reading Frame {i} ({direction}):\n{rna}\n"
            proteins = find_all_proteins(rna)
            for idx, prot in enumerate(proteins, 1):
                output += f"   🔹 Protein {idx}: {prot}\n"

    return output

# --- Streamlit Interface ---
st.set_page_config(page_title="TranscriptPro", layout="centered")
st.title("🧬 TranscriptPro")
st.subheader("DNA/RNA Sequence Translator & Analyzer")

seq_type = st.radio("Select Sequence Type:", ("DNA", "RNA"))
seq_input = st.text_area("Paste your DNA or RNA sequence below:", height=200)

if st.button("🧪 Translate Sequence"):
    try:
        seq = seq_input.strip().upper()
        output = process_sequence(seq, seq_type)

        st.success("✅ Translation Completed")
        st.text_area("🧾 Output:", output, height=500)

    except ValueError as ve:
        st.error(f"❌ Error: {ve}")
