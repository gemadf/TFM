import itertools

import pandas as pd
import string
import time
import difflib
from collections import defaultdict

def readData():
    data = pd.read_excel("data_nervous_genes.xlsx")
    sequences = data["protein_sequence"].head(1000)

    return sequences


def buscar_patrones(sequences):
    all_patterns = dict()
    for protein in sequences:
        all_patterns[protein] = []
        protein_len = len(protein)
        patterns = set()
        posicionPatterns = dict()

        #La primera pasada guarda los patrones de longuitud 1 con frecuencia minima, y las posiciones donde aparecen
        guardar_patrones(protein, patterns, posicionPatterns)

        # Comprueba si el diccionario con la posiciones de patrones NO esta vacío
        if bool(posicionPatterns):
            for pattern_length in range(2, protein_len + 1):
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
                    if len(pos) >= min_frequence:
                        patterns.add(seq)
                        posicionPatterns[seq] = posicionPatterns.get(seq, []) + pos


        #Ordenar de mayor a menor tamaño. Las subcadenas del mismo tamaño se ordenan por orden alfabetico
        list_ordered_patterns = sorted(patterns, key=lambda x: (-len(x), x))
        all_patterns[protein] = list_ordered_patterns

    with open("prueba.txt", "w") as archivo:
       for key, value in all_patterns.items():
          print(value, key, file=archivo)

def guardar_patrones(protein, patterns, posicionPatterns):
    aux = set()
    auxPos = dict()
    for index, letter in enumerate(protein):
        if letter not in aux:
            aux.add(letter)
        auxPos[letter] = auxPos.get(letter, []) + [index]

    for key, value in auxPos.items():
        if len(value) >= min_frequence:
            patterns.add(key)
            posicionPatterns[key] = posicionPatterns.get(key, []) + value



if __name__ == "__main__":
    inicio = time.time()
    #sequences = ['MTILFLTMVISYFGCMKAAPMKEANIRGQGGLAYPGVRTHGTLESVNGPKAGSRGLTSLADTFEHVIEELLDEDQKVRPNEENNKDADLYTSRVMLSSQVPLEPPLLFLLEEYKNYLDAANMSMRVRRHSDPARRGELSVCDSISEWVTAADKKTAVDMSGGTVTVLEKVPVSKGQLKQYFYETKCNPMGYTKEGCRGIDKRHWNSQCRTTQSYVRALTMDSKKRIGWRFIRIDTSCVCTLTIKRGR']
    #sequences = ['ABCABNABHABMABYAB']
    min_frequence = 5
    sequences = readData()
    buscar_patrones(sequences)

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