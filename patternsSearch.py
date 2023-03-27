import itertools

import pandas as pd
import string
import time
import difflib


def readData():
    data = pd.read_excel("data_nervous_genes.xlsx")
    sequences = data["protein_sequence"].head(100)

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
                #print("1")
                #Se crean nuevos aux en cada iteracion. Asi solo contiene los patrones de longuitud pattern_length
                aux = set()
                auxPos = dict(posicionPatterns)
                for key, value in posicionPatterns.items():
                    #print("2")
                    if len(key) == pattern_length - 1:
                        for position in value:
                            #print("3")
                            if (protein_len < position + pattern_length):
                                continue
                            sub_seq = protein[position:position + pattern_length]
                            if sub_seq not in aux:
                                aux.add(sub_seq)
                            auxPos[sub_seq] = auxPos.get(sub_seq, []) + [position]

                for key, value in auxPos.items():
                    #print("4")
                    if len(value) >= min_frequence:
                        #print("5")
                        patterns.add(key)
                        posicionPatterns[key] = posicionPatterns.get(key, []) + [value]

                #Si no se encuentra ningun patron de longuitud pattern_length se sale del bucle. No hay mas patrones posible a encontrar
                if not bool(auxPos):
                    break
                #Se añaden los patrones encontrados a posicionPatterns
                posicionPatterns.update(auxPos)

        #Ordenar de mayor a menor tamaño. Las subcadenas del mismo tamaño se ordenan por orden alfabetico
        list_ordered_patterns = sorted(patterns, key=lambda x: (-len(x), x))
        all_patterns[protein] = list_ordered_patterns

    with open("prueba.txt", "w") as archivo:
       for key, value in all_patterns.items():
          print(f"Protein sequence: {key}, patterns: {value}", file=archivo)

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
    #sequences = ['MSLWQPLVLVLLVLGCCFAAPRQRQSTLVLFPGDLRTNLTDRQLAEEYLYRYGYTRVAEMRGESKSLGPALLLLQKQLSLPETGELDSATLKAMRTPRCGVPDLGRFQTFEGDLKWHHHNITYWIQNYSEDLPRAVIDDAFARAFALWSAVTPLTFTRVYSRDADIVIQFGVAEHGDGYPFDGKDGLLAHAFPPGPGIQGDAHFDDDELWSLGKGVVVPTRFGNADGAACHFPFIFEGRSYSACTTDGRSDGLPWCSTTANYDTDDRFGFCPSERLYTQDGNADGKPCQFPFIFQGQSYSACTTDGRSDGYRWCATTANYDRDKLFGFCPTRADSTVMGGNSAGELCVFPFTFLGKEYSTCTSEGRGDGRLWCATTSNFDSDKKWGFCPDQGYSLFLVAAHEFGHALGLDHSSVPEALMYPMYRFTEGPPLHKDDVNGIRHLYGPRPEPEPRPPTTTTPQPTAPPTVCPTGPPTVHPSERPTAGPTGPPSAGPTGPPTAGPSTATTVPLSPVDDACNVNIFDAIAEIGNQLYLFKDGKYWRFSEGRGSRPQGPFLIADKWPALPRKLDSVFEERLSKKLFFFSGRQVWVYTGASVLGPRRLDKLGLGADVAQVTGALRSGRGKMLLFSGRRLWRFDVKAQMVDPRSASEVDRMFPGVPLDTHDVFQYREKAYFCQDRFYWRVSSRSELNQVDQVGYVTYDILQCPED']
    #sequences = ['ELERKELEFDTNMDAVQMVITEAQKVDTRAKNAGVTIQDTLNTLDGLLHLMDQPLSVDEEGLVLLEQKLSRAKTQINSQLRPMMSELEERARQQR']
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