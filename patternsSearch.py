import itertools

import pandas as pd
import string

def readData():
    data = pd.read_excel("data_nervous_genes.xlsx")
    sequences = data["protein_sequence"].head()

    return sequences

def buscar_patrones(sequences):
    all_patterns = dict()
    pattern_ocurrence = 5
    for protein in sequences:
        protein_len = len(protein)
        patterns = set()
        for pattern in range(1, protein_len+1):
            for i in range(protein_len - pattern + 1):
                sub_seq = protein[i:i+pattern]
                if protein.count(sub_seq) >= pattern_ocurrence:
                    patterns.add(sub_seq)
                    all_patterns[protein] = patterns

    for key, value in all_patterns.items():
        print(f"Protein sequence: {key}, patterns: {value}")


if __name__ == "__main__":
    sequences = ['ABCABCAABABABABAB', 'BCABCBCCC', 'CABCA', 'ACBCA', 'ABCABCABCABCABCBACB']
    #sequences = readData()
    buscar_patrones(sequences)

