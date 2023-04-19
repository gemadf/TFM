import itertools

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


def buscar_patrones_cada_proteina(sequences):
    for protein in sequences:
        all_patterns[protein] = []
        protein_len = len(protein)
        patterns = set() #Guarda solo los patrones que aparecen en la secuencia
        posicionPatterns = dict() #Guarda los patrones que aparecen en la secuencia con sus posiciones asociadas

        #La primera pasada guarda los patrones de longuitud 1 con frecuencia minima, y las posiciones donde aparecen
        guardar_patrones_len1(protein, patterns, posicionPatterns)

        # Comprueba si el diccionario con la posiciones de patrones NO esta vacío
        if bool(posicionPatterns):
            #for pattern_length in range(2, protein_len + 1):
            for pattern_length in range(2,6):
                print(pattern_length)
                #Si se intenta acceder a una clave que no existe se creara una lista vacia
                auxPos = defaultdict(list)
                for key, value in posicionPatterns.items():
                    if len(key) == pattern_length - 1:
                        for position in value:
                            if (protein_len < position + pattern_length):
                                continue
                            sub_seq = protein[position:position + pattern_length]
                            #Si la ultima letra que es la nueva del patron ya esta min_freq, el patron es min freq tb
                            ultima_letra = sub_seq[-1]
                            pos_ultima_letra = position + pattern_length - 1
                            if ultima_letra in patterns and pos_ultima_letra in posicionPatterns[ultima_letra]:
                                auxPos[sub_seq].append(position)


                # Si no se encuentra ningun patron de longuitud pattern_length se sale del bucle. No hay mas patrones posible a encontrar
                if not bool(auxPos):
                    break

                for seq, pos in auxPos.items():
                    #if len(pos) >= min_frequence:
                    patterns.add(seq)
                    posicionPatterns[seq] = posicionPatterns.get(seq, []) + pos

        #Ordenar de mayor a menor tamaño. Las subcadenas del mismo tamaño se ordenan por orden alfabetico
        dict_ordered_patterns = dict(sorted(posicionPatterns.items(), key=lambda x: (-len(x[0]), x[0])))
        all_patterns[protein] = dict_ordered_patterns

    """with open("patternsAndAsocciatedPositionsEachProtein1000.txt", "w") as archivo:
       for key, value in all_patterns.items():
          print(value, key, file=archivo)"""

    return all_patterns

def guardar_patrones_len1(protein, patterns, posicionPatterns):
    aux = set()
    auxPos = dict()
    for index, letter in enumerate(protein):
        if letter not in aux:
            aux.add(letter)
        auxPos[letter] = auxPos.get(letter, []) + [index]

    for key, value in auxPos.items():
        #if len(value) >= min_frequence:
        patterns.add(key)
        posicionPatterns[key] = posicionPatterns.get(key, []) + value


def patrones_freq_min_identicos(all_patterns):
    # Each pattern associated to the proteins the pattern is in
    pattern_proteins = {}

    for protein, patterns in all_patterns.items():
        for pattern, positions in patterns.items():
            if pattern not in pattern_proteins:
                pattern_proteins[pattern] = set()
            pattern_proteins[pattern].add(protein)

    pattern_freqMin = {}
    for pattern, proteins in pattern_proteins.items():
        if len(proteins) > min_ocurrence:
            pattern_freqMin[pattern] = proteins

    # Ordenar de mayor a menor tamaño. Las subcadenas del mismo tamaño se ordenan por orden alfabetico
    dict_ordered_patterns = dict(sorted(pattern_freqMin.items(), key=lambda x: (-len(x[0]), x[0])))

    with open("1000patternsIdenticos_frequMin5.txt", "w") as archivo:
        for key, value in dict_ordered_patterns.items():
            print(key, value, file=archivo)

    df = pd.DataFrame(pattern_freqMin.items(), columns=['pattern', 'proteins'])
    df.to_csv('1000patternsIdenticos_frequMin5.csv', index=False)

def patrones_similares(all_patterns):
    pattern_proteins = {}  # Each pattern associated to the proteins the pattern is in
    similar_patterns = {}  # Guarda los patrones similares relacionados con el patron similar del que parten
    pattern_freqMin = {}  # Guarda para cada patron similar las proteinas en las que se ha encontrado


    for protein, patterns in all_patterns.items():
        for pattern, positions in patterns.items():
            if pattern not in pattern_proteins:
                pattern_proteins[pattern] = set()
            pattern_proteins[pattern].add(protein)

    for pattern1, proteins1 in pattern_proteins.items():
        for pattern2, proteins2 in pattern_proteins.items():
            if pattern1 != pattern2 and pattern2 not in similar_patterns:
                # Calcular distancia de Levenshtein entre patrones
                similarity = Levenshtein.distance(pattern1, pattern2) / max(len(pattern1), len(pattern2))
                # 0.1 indica que admite 1 insercion, deleccion o sustitucion
                if similarity >= 0.1:
                    # Verificar min_freq de las proteinas proteínas contengan ambos patrones
                    if len(proteins1 & proteins2) >= min_ocurrence:
                        if pattern1 not in similar_patterns:
                            similar_patterns[pattern1] = set()
                        if pattern2 not in similar_patterns:
                            similar_patterns[pattern2] = set()
                        similar_patterns[pattern1].add(pattern2)
                        similar_patterns[pattern2].add(pattern1)

    # Actualizar con las proteinas donde aparece su patron similar
    """for pattern1, patterns_similares in similar_patterns.items():
        for pattern2 in patterns_similares:
            if pattern1 in pattern_proteins:
                pattern_freqMin[pattern1].update(pattern_proteins[pattern1])
            if pattern2 in pattern_proteins:
                pattern_freqMin[pattern2].update(pattern_proteins[pattern2])"""

    # Patrones similares relacionados con el patron similar del que parten
    # Ordenar de mayor a menor tamaño. Las subcadenas del mismo tamaño se ordenan por orden alfabetico
    dict_ordered_patterns = dict(sorted(similar_patterns.items(), key=lambda x: (-len(x[0]), x[0])))

    with open("10patternsSimilaresYSusPatronesRelacionados_freqMin2.txt", "w") as archivo:
        for key, value in dict_ordered_patterns.items():
            print(key, value, file=archivo)

    df = pd.DataFrame(similar_patterns.items(), columns=['pattern', 'proteins'])
    df.to_csv('10patternsSimilaresYSusPatronesRelacionados_freqMin2.csv', index=False)


    # Para cada patron similar guarda las proteinas en las que se ha encontrado
    # Ordenar de mayor a menor tamaño. Las subcadenas del mismo tamaño se ordenan por orden alfabetico
    """dict_ordered_patterns = dict(sorted(pattern_freqMin.items(), key=lambda x: (-len(x[0]), x[0])))

    with open("10patternsSimilaresProteinas_freqMin2.txt", "w") as archivo:
        for key, value in dict_ordered_patterns.items():
            print(key, value, file=archivo)

    df = pd.DataFrame(pattern_freqMin.items(), columns=['pattern', 'proteins'])
    df.to_csv('10patternsSimilaresProteinas_freqMin2.csv', index=False)"""


def distance_levenshtein(pattern1, pattern2):
    d = dict()
    for i in range(len(pattern1)+1):
        d[i] = dict()
        d[i][0] = i
    for i in range(len(pattern2)+1):
        d[0][i] = i
    for i in range(1, len(pattern1)+1):
        for j in range(1, len(pattern2)+1):
            d[i][j] = min(d[i][j-1]+1,
                          d[i-1][j]+1,
                          d[i-1][j-1]+(not pattern1[i-1] == pattern2[j-1]))
    return d[len(pattern1)][len(pattern2)]

if __name__ == "__main__":
    inicio = time.time()
    all_patterns = dict()
    #sequences = ['MTILFLTMVISYFGCMKAAPMKEANIRGQGGLAYPGVRTHGTLESVNGPKAGSRGLTSLADTFEHVIEELLDEDQKVRPNEENNKDADLYTSRVMLSSQVPLEPPLLFLLEEYKNYLDAANMSMRVRRHSDPARRGELSVCDSISEWVTAADKKTAVDMSGGTVTVLEKVPVSKGQLKQYFYETKCNPMGYTKEGCRGIDKRHWNSQCRTTQSYVRALTMDSKKRIGWRFIRIDTSCVCTLTIKRGR']
    #sequences = ['ABCABNABHABMABYAB']
    min_ocurrence = 5
    sequences = readData()
    all_patterns = buscar_patrones_cada_proteina(sequences)

    patrones_freq_min_identicos(all_patterns)
    #patrones_similares(all_patterns)

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