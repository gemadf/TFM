import pandas as pd
import time
import ast
import csv
import math
from interfazGrafica import interfaz
from descarteProteinas import ejecutar
import metricas
from graficas import grafica

def readData(archivoEntrada, enfermedad, archivoTarget):
    data = pd.read_excel(archivoEntrada)
    dataC = pd.read_csv("resultados/proteinasDescartadas.csv")

    #Descarte de proteinas
    data = data[~data['protein_id'].isin(dataC['ProteinasDescartadas'])]

    # "C0002395"
    if(enfermedad != ''):
        data = data.loc[data["disease_id"] == enfermedad]
        dataB = pd.read_excel("proteinas_en_comun_Alzheimer.xlsx")

    if(archivoTarget != ''):
        #Eliminar las proteinas target
        data = data[~((data["disease_id"] == enfermedad) &
                      (data["protein_id"].isin(dataB["protein_id"])))]

    sequences = data["protein_sequence"]

    num_filas = sequences.shape[0]

    return sequences, num_filas

def guardar_patrones_len1(sequences, pattern_freqMin):
    all_patterns = dict()
    longitud_max = 0
    # Each pattern associated to the proteins the pattern is in
    pattern_proteins = {}
    for protein in sequences:
        longitud = len(protein)
        if longitud > longitud_max:
            longitud_max = longitud

        all_patterns[protein] = []
        # En cada iteración guarda los patrones que aparecen en la secuencia con sus posiciones asociadas a la proteina
        posicionPatterns = dict()
        for index, letter in enumerate(protein):
            posicionPatterns[letter] = posicionPatterns.get(letter, []) + [index]

        all_patterns[protein] = posicionPatterns


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

    df = pd.DataFrame(pattern_freqMin.items(), columns=['pattern', 'proteins'])
    df.to_csv('prueba2.csv', index=False)
    return pattern_freqMin, posicionPatterns, longitud_max

def buscar_patrones_identicos(sequences):
    pattern_freqMin = {}
    pattern_freqMin, posicionPatterns, longitud_max = guardar_patrones_len1(sequences, pattern_freqMin)

    if bool(pattern_freqMin):
        for pattern_length in range(2, longitud_max + 1):
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

            # Si no se encuentra ningun patron de longitud pattern_length se sale del bucle. No hay mas patrones posible a encontrar
            if not bool(auxPos):
                break

            for pattern, proteins in auxPos.items():
                for prot, pos in proteins.items():
                    if pattern not in pattern_freqMin:
                        pattern_freqMin[pattern] = {}
                    if prot not in pattern_freqMin[pattern]:
                        pattern_freqMin[pattern][prot] = []
                    pattern_freqMin[pattern][prot].extend(pos)
                    if len(pattern) > 2:
                        if pattern[:-1] in pattern_freqMin:
                            del pattern_freqMin[pattern[:-1]]
                        if pattern[1:] in pattern_freqMin:
                            del pattern_freqMin[pattern[1:]]



        # Ordenar de mayor a menor tamaño. Las subcadenas del mismo tamaño se ordenan por orden alfabetico
        dict_ordered_patterns = dict(sorted(pattern_freqMin.items(), key=lambda x: (-len(x[0]), x[0])))

        df = pd.DataFrame(dict_ordered_patterns.items(), columns=['pattern', 'proteins'])

        df.to_csv('prueba.csv', index=False)

    return pattern_freqMin

def remplazar_sequence_for_ID(pattern_freqMin):
    df_b = pd.read_excel("data_nervous_genes.xlsx")

    output = []

    for key, value in pattern_freqMin.items():
        for proteina, posiciones in value.items():
            output.append([key, proteina, posiciones])

    output = [sublista for sublista in output if len(sublista[0]) != 1]

    # Ordenar de mayor a menor tamaño. Las subcadenas del mismo tamaño se ordenan por orden alfabetico
    output_ordered = sorted(output, key=lambda x: (-len(x[0]), x[0]))


    proteinas_dict = dict(df_b[['protein_sequence', 'protein_id']].values)

    for item in output_ordered:
        protein_sequence = item[1]
        if protein_sequence in proteinas_dict:
            item[1] = proteinas_dict[protein_sequence]

    df_a = pd.DataFrame(output_ordered, columns=['Patron', 'Proteina', 'Posiciones'])

    # Guardar el DataFrame actualizado en un archivo CSV
    df_a.to_csv('resultados/patronesIdenticos.csv', index=False)

if __name__ == "__main__":
    inicio = time.time()

    datosInterfaz = interfaz()
    print(datosInterfaz)

    archivoEntrada = datosInterfaz["NombreArchivoEntrada"]
    enfermedad = datosInterfaz["CodigoEnfermedad"]
    archivoTarget = datosInterfaz["NombreArchivoTarget"]
    similitud = float(datosInterfaz["Similitud"])

    ejecutar(archivoEntrada, enfermedad, similitud)
    pattern_freqMin = dict()

    sequences, num_filas = readData(archivoEntrada, enfermedad, archivoTarget)

    min_ocurrence = math.floor(num_filas * float(datosInterfaz["OcurrenciaMin"]))

    pattern_freqMin = buscar_patrones_identicos(sequences)
    remplazar_sequence_for_ID(pattern_freqMin)
    datosInterfaz = dict()

    metricas.metrica_distanciaProteinas()
    archivo = 'resultados/Metrica_distanciaProteinasMismoPatron.csv'
    nombreOutput = 'resultados/Figura_DistanciaProteinasMismoPatron'
    grafica(archivo, nombreOutput)


    metricas.patronesComun()
    archivo = 'resultados/Metrica_patronesComunes.csv'
    nombreOutput = 'resultados/Figura_distanciaProteinasPatronesComunes'
    grafica(archivo, nombreOutput)


    fin = time.time()

    tiempo_total = fin - inicio
    print(tiempo_total, "segundos")