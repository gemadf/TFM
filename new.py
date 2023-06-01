import pandas as pd
import string
import time
import difflib
from collections import defaultdict
import Levenshtein
import ast

def readData():
    data = pd.read_excel("data_nervous_genes.xlsx")
    sequences = data["protein_sequence"].head(5)

    """count_by_disease = data.groupby('disease_id').size().reset_index(name='count')
    sorted_by_count = count_by_disease.sort_values('count', ascending=False)

    print(sorted_by_count.head(5))

    data = data.loc[data["disease_id"] == "C0002395"]

    dataB = pd.read_excel("proteinas_en_comun_Alzheimer.xlsx")

    data = data[~((data["disease_id"] == "C0002395") &
                  (data["protein_id"].isin(dataB["protein_id"])) &
                  (data["gene_id"].isin(dataB["gene_id"])))]
    sequences = data["protein_sequence"]"""

    return sequences


def guardar_patrones_len1(sequences, pattern_freqMin):
    all_patterns = dict()
    longuitud_max = 0
    for protein in sequences:
        longuitud = len(protein)
        if longuitud > longuitud_max:
            longuitud_max = longuitud
        aux = set()
        auxPos = dict()
        all_patterns[protein] = []
        posicionPatterns = dict()  # Guarda los patrones que aparecen en la secuencia con sus posiciones asociadas
        patterns = set()  # Guarda solo los patrones que aparecen en la secuencia
        for index, letter in enumerate(protein):
            if letter not in aux:
                aux.add(letter)
            auxPos[letter] = auxPos.get(letter, []) + [index]

        for key, value in auxPos.items():
            patterns.add(key)
            posicionPatterns[key] = posicionPatterns.get(key, []) + value

        all_patterns[protein] = posicionPatterns

    # Each pattern associated to the proteins the pattern is in
    pattern_proteins = {}
    for protein, patterns in all_patterns.items():
        for pattern, positions in patterns.items():
            if pattern not in pattern_proteins:
                pattern_proteins[pattern] = {}
            if protein not in pattern_proteins[pattern]:
                pattern_proteins[pattern][protein] = []
            pattern_proteins[pattern][protein].extend(positions)


    for pattern, proteins in pattern_proteins.items():
        if len(proteins) >= min_ocurrence:
            pattern_freqMin[pattern] = proteins

    return pattern_freqMin, posicionPatterns, longuitud_max

def buscar_patrones_cada_proteina(sequences):
    pattern_freqMin = {}
    pattern_freqMin, posicionPatterns, longuitud_max = guardar_patrones_len1(sequences, pattern_freqMin)

    if bool(pattern_freqMin):
        for pattern_length in range(2, longuitud_max + 1):
            # Si se intenta acceder a una clave que no existe se creara una lista vacia
            auxPos = {}
            sub_seqs = []
            for pattern, proteins in pattern_freqMin.items():
                if len(pattern) == pattern_length - 1:
                    for prot, positions in proteins.items():
                        protein_len = len(prot)
                        if protein_len < pattern_length - 1:
                            continue
                        for position in positions:
                            if (protein_len < position + pattern_length):
                                continue
                            sub_seq = prot[position:position + pattern_length]

                            if sub_seq in pattern_freqMin:
                                continue
                            # Si la ultima letra que es la nueva del patron ya esta min_freq, el patron es posible
                            # min freq tb
                            ultima_letra = sub_seq[-1]
                            pos_ultima_letra = position + pattern_length - 1
                            if ultima_letra in pattern_freqMin and pos_ultima_letra in pattern_freqMin[ultima_letra][prot]:
                                if sub_seq not in auxPos:
                                    auxPos[sub_seq] = {}
                                if prot not in auxPos[sub_seq]:
                                    auxPos[sub_seq][prot] = []
                                auxPos[sub_seq][prot].append(position)
                                if sub_seq not in sub_seqs:
                                    sub_seqs.append(sub_seq)

                sub_seqs_copy = sub_seqs.copy()
                for p in sub_seqs_copy:
                    if len(auxPos[p]) < min_ocurrence:
                        del auxPos[p]
                        sub_seqs.remove(p)



            # Si no se encuentra ningun patron de longuitud pattern_length se sale del bucle. No hay mas patrones posible a encontrar
            if not bool(auxPos):
                break

            for pattern, proteins in auxPos.items():
                for prot, pos in proteins.items():
                    if pattern not in pattern_freqMin:
                        pattern_freqMin[pattern] = {}
                    if prot not in pattern_freqMin[pattern]:
                        pattern_freqMin[pattern][prot] = []
                    pattern_freqMin[pattern][prot].extend(pos)



        # Ordenar de mayor a menor tamaño. Las subcadenas del mismo tamaño se ordenan por orden alfabetico
        dict_ordered_patterns = dict(sorted(pattern_freqMin.items(), key=lambda x: (-len(x[0]), x[0])))

        df = pd.DataFrame(dict_ordered_patterns.items(), columns=['pattern', 'proteins'])
        df.to_csv('prueba.csv', index=False)

    return pattern_freqMin


def remplazar_sequence_for_ID():
    df_a = pd.read_csv('prueba.csv')
    df_b = pd.read_excel("data_nervous_genes.xlsx")
    proteinas_dict = dict(df_b[['protein_sequence', 'protein_id']].values)

    for i, row in df_a.iterrows():
        proteins_str = row['proteins']
        proteins_dict = ast.literal_eval(proteins_str)
        new_proteins_dict = {}
        for protein, positions in proteins_dict.items():
            if protein in proteinas_dict:
                new_protein = proteinas_dict[protein]
                new_proteins_dict[new_protein] = positions
            else:
                new_proteins_dict[protein] = positions
        df_a.at[i, 'proteins'] = str(new_proteins_dict)

    # Guardar el DataFrame actualizado en un archivo CSV
    df_a.to_csv('results.csv', index=False)


def patrones_similares(sequences):
    similar_patterns = {}  # Guarda los patrones similares relacionados con el patron similar del que parten
    pattern_freqMin = {}  # Guarda para cada patron similar las proteinas en las que se ha encontrado


    for pattern1 in sequences:
        for pattern2 in sequences:
            if pattern1 != pattern2:
                # Calcular distancia de Levenshtein entre patrones
                similarity = Levenshtein.distance(pattern1, pattern2) / max(len(pattern1), len(pattern2))
                # Para admitir una inserción o una delección el valor debe ser 1, y dividimos para normalizar y
                # adaptarlo a las distintas longuitudes
                umbral = 1 / max(len(pattern1), len(pattern2))
                #print("Patron 1: ", pattern1, " Patron 2: ", pattern2, " Similariad: ", similarity)
                #print(umbral)
                if similarity <= umbral:
                    if pattern1 not in similar_patterns:
                        similar_patterns[pattern1] = set()
                    if pattern2 not in similar_patterns:
                        similar_patterns[pattern2] = set()
                    similar_patterns[pattern1].add(pattern2)
                    similar_patterns[pattern2].add(pattern1)

    print(similar_patterns)



if __name__ == "__main__":
    inicio = time.time()
    pattern_freqMin = dict()
    #sequences = ['MTL', 'MTAL', 'MTOALX', 'MTJ', 'MTIO', 'MTK', 'AL', 'OIALCP', 'ALP', 'PALMT', 'ALX', 'ALZ', 'ALWQ']
    #sequences = ['MTL', 'MTAL', 'MTAL', 'MTJ', 'MTIO', 'MTK']
    #sequences = ['ABCABNABHABMABYAB']
    min_ocurrence = 5
    sequences = readData()
    pattern_freqMin = buscar_patrones_cada_proteina(sequences)
    #remplazar_sequence_for_ID()

    #lista_similares = ['ABCDE', 'ABCEE', 'ABCDF', 'A', 'ABCE', 'ABCDEF']
    #patrones_similares(sequences)
    fin = time.time()

    tiempo_total = fin - inicio
    print(tiempo_total, "segundos")

    """with open('prueba1.txt', 'r') as f1, open('prueba2.txt', 'r') as f2:
        contenido1 = f1.read()
        contenido2 = f2.read()

    diff = difflib.unified_diff(contenido1.splitlines(), contenido2.splitlines())

    if contenido1 == contenido2:
        print("Los archivos son iguales")
    else:
        print("Los archivos son diferentes")
        print('\n'.join(diff))"""