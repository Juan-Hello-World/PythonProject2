from itertools import product
from time import time
import pandas as pd

# 1. BASE DE DATOS DE AMINOACIDOS Y PROPIEDADES BIOFISICAS
amino_acids = {
    'A': {'mass': 89, 'hydrophobicity': 1.8, 'polarity': 'nonpolar'},
    'C': {'mass': 121, 'hydrophobicity': 2.5, 'polarity': 'polar'},
    'D': {'mass': 133, 'hydrophobicity': -3.5, 'polarity': 'negative'},
    'E': {'mass': 147, 'hydrophobicity': -3.5, 'polarity': 'negative'},
    'F': {'mass': 165, 'hydrophobicity': 2.8, 'polarity': 'nonpolar'},
    'G': {'mass': 75, 'hydrophobicity': -0.4, 'polarity': 'nonpolar'},
    'H': {'mass': 155, 'hydrophobicity': -3.2, 'polarity': 'positive'},
    'I': {'mass': 131, 'hydrophobicity': 4.5, 'polarity': 'nonpolar'},
    'K': {'mass': 146, 'hydrophobicity': -3.9, 'polarity': 'positive'},
    'L': {'mass': 131, 'hydrophobicity': 3.8, 'polarity': 'nonpolar'},
    'M': {'mass': 149, 'hydrophobicity': 1.9, 'polarity': 'nonpolar'},
    'N': {'mass': 132, 'hydrophobicity': -3.5, 'polarity': 'polar'},
    'P': {'mass': 115, 'hydrophobicity': -1.6, 'polarity': 'nonpolar'},
    'Q': {'mass': 146, 'hydrophobicity': -3.5, 'polarity': 'polar'},
    'R': {'mass': 174, 'hydrophobicity': -4.5, 'polarity': 'positive'},
    'S': {'mass': 105, 'hydrophobicity': -0.8, 'polarity': 'polar'},
    'T': {'mass': 119, 'hydrophobicity': -0.7, 'polarity': 'polar'},
    'V': {'mass': 117, 'hydrophobicity': 4.2, 'polarity': 'nonpolar'},
    'W': {'mass': 204, 'hydrophobicity': -0.9, 'polarity': 'nonpolar'},
    'Y': {'mass': 181, 'hydrophobicity': -1.3, 'polarity': 'polar'}
}

# Ordenamiento lexicográfico
sorted_aas = sorted(amino_acids.keys())

# 2. VALIDACION DE PEPTIDOS
def is_valid_peptide(sequence):
    """Check if peptide meets all biological constraints"""
    # 1. No mas de 2 repeticiones consecutivas
    for i in range(len(sequence) - 2):
        if sequence[i] == sequence[i + 1] == sequence[i + 2]:
            return False
    # Calcular propiedades
    mass = sum(amino_acids[aa]['mass'] for aa in sequence)
    avg_hydro = sum(amino_acids[aa]['hydrophobicity'] for aa in sequence) / len(sequence)
    nonpolar_count = sum(1 for aa in sequence if amino_acids[aa]['polarity'] == 'nonpolar')
    # 2. Restricción de masa (500-2000 Da)
    if not 500 <= mass <= 2000:
        return False
    # 3. Hidrofobicidad (>1.0)
    if avg_hydro <= 1.0:
        return False
    # 4. Relacion no polar (≤50%)
    if nonpolar_count / len(sequence) > 0.5:
        return False
    return True


# 3. GENERADOR DE PÉPTIDOS CON ANÁLISIS
def generate_peptides(k, max_results=50, verbose=True):
    """Generate and validate peptides with performance tracking"""
    results = []
    start_time = time()
    if verbose:
        print(f"\nGenerating peptides of length {k}...")
        print("Constraints:")
        print("- Mass: 500-2000 Da")
        print("- Avg hydrophobicity > 1.0")
        print("- ≤2 consecutive repeats")
        print("- ≤50% nonpolar amino acids")
        print("\nProcessing...\n")
    for seq in product(sorted_aas, repeat=k):
        if is_valid_peptide(seq):
            seq_str = ''.join(seq)
            mass = sum(amino_acids[aa]['mass'] for aa in seq)
            hydro = sum(amino_acids[aa]['hydrophobicity'] for aa in seq) / k
            nonpolar = sum(1 for aa in seq if amino_acids[aa]['polarity'] == 'nonpolar') / k
            results.append({
                'sequence': seq_str,
                'mass': mass,
                'hydrophobicity': round(hydro, 2),
                'nonpolar_ratio': round(nonpolar, 2)
            })
            if len(results) >= max_results:
                break
    # Performance analysis
    total_time = round(time() - start_time, 2)
    total_sequences = len(sorted_aas) ** k
    valid_ratio = round(len(results) / max_results * 100, 2) if max_results else 0
    if verbose:
        print(f"\nCOMPLETED IN {total_time}s")
        print(f"Total possible sequences: {total_sequences:,}")
        print(f"Valid peptides found: {len(results)}")
        print(f"Validation ratio: {valid_ratio}%")
    return pd.DataFrame(results)


# 4. EJECUCION Y EXPORTACION DE RESULTADOS
if __name__ == "__main__":
    print("=== ANTIMICROBIAL PEPTIDE DESIGN TOOL ===")
    print(f"Available amino acids: {', '.join(sorted_aas)}\n")

    # Parameters (adjustable)
    PEPTIDE_LENGTH = 7
    MAX_RESULTS = 50

    # Run generator
    results_df = generate_peptides(
        k=PEPTIDE_LENGTH,
        max_results=MAX_RESULTS
    )
    # Display and save results
    print("\nTOP CANDIDATES:")
    print(results_df.head(10))
    # Export to CSV
    results_df.to_csv(f"peptide_candidates_k{PEPTIDE_LENGTH}.csv", index=False)
    print(f"\nResults saved to 'peptide_candidates_k{PEPTIDE_LENGTH}.csv'")