import time
import Levenshtein
import math

def patrones_similares(pattern_freqMin):
    similar_patterns = {}  # Guarda los patrones similares relacionados con el patron similar del que parten
    num_op = 3
    similar_patterns = {}  # Guarda los patrones similares relacionados con el patron similar del que parten

    patterns = list(pattern_freqMin.keys())
    num_patterns = len(patterns)

    for i in range(num_patterns):
        pattern1 = patterns[i]
        proteins1 = pattern_freqMin[pattern1]
        len_pattern1 = len(pattern1)

        for j in range(i + 1, num_patterns):
            pattern2 = patterns[j]
            proteins2 = pattern_freqMin[pattern2]
            len_pattern2 = len(pattern2)

            # Calcular distancia de Levenshtein entre patrones
            similarity = Levenshtein.distance(pattern1, pattern2) / max(len(pattern1), len(pattern2))
            # Para admitir una inserción, una delección o una sustitución el valor debe ser 1, y dividimos para normalizar y
            # adaptarlo a las distintas longitudes
            max_length = max(len_pattern1, len_pattern2)
            operaciones_max = math.ceil(0.1 * max_length)
            umbral = operaciones_max / max_length
            # print("Patron 1: ", pattern1, " Patron 2: ", pattern2, " Similariad: ", similarity)
            # print(umbral)

            if similarity <= umbral:
                if pattern1 not in similar_patterns:
                    similar_patterns[pattern1] = set()
                if pattern2 not in similar_patterns:
                    similar_patterns[pattern2] = set()

                if pattern1 not in pattern_freqMin:
                    pattern_freqMin[pattern2] = {}

                for proteina, posiciones in proteins1.items():
                    if proteina not in pattern_freqMin[pattern2]:
                        pattern_freqMin[pattern2][proteina] = []
                        if posiciones:
                            pattern_freqMin[pattern2][proteina].extend(posiciones)
                    else:
                        for posicion in posiciones:
                            if posicion not in pattern_freqMin[pattern2][proteina]:
                                pattern_freqMin[pattern2][proteina].append(posicion)
                            pattern_freqMin[pattern2][proteina].sort()

                if pattern2 not in pattern_freqMin:
                    pattern_freqMin[pattern1] = {}

                for proteina, posiciones in proteins2.items():
                    if proteina not in pattern_freqMin[pattern1]:
                        pattern_freqMin[pattern1][proteina] = []
                        if posiciones:
                            pattern_freqMin[pattern1][proteina].extend(posiciones)
                    else:
                        for posicion in posiciones:
                            if posicion not in pattern_freqMin[pattern1][proteina]:
                                pattern_freqMin[pattern1][proteina].append(posicion)
                            pattern_freqMin[pattern1][proteina].sort()

                similar_patterns[pattern1].add(pattern2)
                similar_patterns[pattern2].add(pattern1)
    print(pattern_freqMin)
    return pattern_freqMin


if __name__ == "__main__":
    inicio = time.time()

    lista_similares = {'Al': {"ALFAS": [0], 'AFNERALDAL': [5, 7]}, 'AF': {"JDLKAJSDL": [8, 9, 10], "JLADA": [3, 7], "ALFAS": [6]}}
    patrones_similares(lista_similares)
    fin = time.time()

    tiempo_total = fin - inicio
    print(tiempo_total, "segundos")