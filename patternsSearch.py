import itertools

import pandas as pd
import string


def readData():
    data = pd.read_excel("data_nervous_genes.xlsx")
    sequences = data["protein_sequence"].head(100)

    return sequences


def buscar_patrones(sequences):
    all_patterns = dict()
    pattern_ocurrence = 10
    for protein in sequences:
        protein_len = len(protein)
        patterns = set()
        for pattern_length in range(2, protein_len + 1):
            for pattern in range(protein_len - pattern_length + 1):
                sub_seq = protein[pattern:pattern + pattern_length]
                if protein.count(sub_seq) >= pattern_ocurrence:
                    patterns.add(sub_seq)


        list_ordered_patterns = sorted(patterns, key=len, reverse=True)
        all_patterns[protein] = list_ordered_patterns


    with open("100ProteinasOcurrence10.txt", "w") as archivo:
        for key, value in all_patterns.items():
            print(f"Protein sequence: {key}, patterns: {value}", file=archivo)


if __name__ == "__main__":
    #sequences = ['ABCABCAABABABABABABAB', 'BCABCBCCCCC', 'CABCA', 'ACBCAABCABCABCABCABCABCABCAC', 'ABCABCABCABCABCBACB', '']
    sequences = readData()
    buscar_patrones(sequences)
