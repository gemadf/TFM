import pandas as pd
import string
import time
import difflib
from collections import defaultdict
import Levenshtein

def readData():
    data = pd.read_excel("data_nervous_genes.xlsx")
    sequences = data["protein_sequence"].head(100)

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
            #if len(value) >= min_frequence:
            patterns.add(key)
            posicionPatterns[key] = posicionPatterns.get(key, []) + value



        all_patterns[protein] = posicionPatterns

    # Each pattern associated to the proteins the pattern is in
    pattern_proteins = {}
    for protein, patterns in all_patterns.items():
        for pattern, positions in patterns.items():
            if pattern not in pattern_proteins:
                pattern_proteins[pattern] = {}
            #pattern_proteins[pattern].add(protein)
            if protein not in pattern_proteins[pattern]:
                pattern_proteins[pattern][protein] = []
            pattern_proteins[pattern][protein].extend(positions)


    for pattern, proteins in pattern_proteins.items():
        if len(proteins) >= min_ocurrence:
            pattern_freqMin[pattern] = proteins

    return pattern_freqMin, posicionPatterns, longuitud_max

"""if ultima_letra in pattern_freqMin:
    for proteinaa in pattern_freqMin[ultima_letra]:
        if pos_ultima_letra in pattern_freqMin[ultima_letra][proteinaa]:
            if sub_seq not in auxPos:
                auxPos[sub_seq] = {}
            if proteinaa not in auxPos[sub_seq]:
                auxPos[sub_seq][proteinaa] = []
            auxPos[sub_seq][proteinaa].append(position)
    if len(auxPos[sub_seq]) < min_ocurrence:
        del auxPos[sub_seq]"""

def buscar_patrones_cada_proteina(sequences):
    pattern_freqMin = {}
    pattern_freqMin, posicionPatterns, longuitud_max = guardar_patrones_len1(sequences, pattern_freqMin)

    #freq_aux = pattern_freqMin

    #for protein in sequences:
    #    protein_len = len(protein)
        #patterns = set()  # Guarda solo los patrones que aparecen en la secuencia
        #posicionPatterns = dict()  # Guarda los patrones que aparecen en la secuencia con sus posiciones asociadas

        # Comprueba si el diccionario con la posiciones de patrones NO esta vacío
    if bool(pattern_freqMin):
        for pattern_length in range(2, longuitud_max + 1):
        #for pattern_length in range(2, 6):
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
                            # Si la ultima letra que es la nueva del patron ya esta min_freq, el patron es min freq tb
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

        with open("prueba.txt", "w") as archivo:

            for key, value in dict_ordered_patterns.items():
                print(key, value, file=archivo)

    return pattern_freqMin


if __name__ == "__main__":
    inicio = time.time()
    pattern_freqMin = dict()
    #sequences = ['MTL', 'MTAL', 'MTOALX', 'MTJ', 'MTIO', 'MTK', 'AL', 'OIALCP', 'ALP', 'PALMT', 'ALX', 'ALZ', 'ALWQ']
    #sequences = ['MTL', 'MTAL', 'MTAL', 'MTJ', 'MTIO', 'MTK']
    #sequences = ['ABCABNABHABMABYAB']
    min_ocurrence = 5
    sequences = readData()
    pattern_freqMin = buscar_patrones_cada_proteina(sequences)


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