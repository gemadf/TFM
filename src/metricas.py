import time
import pandas as pd
import Levenshtein
import csv

def readData():
    data = pd.read_excel("data_nervous_genes.xlsx")
    """count_by_disease = data.groupby('disease_id').size().reset_index(name='count')
    sorted_by_count = count_by_disease.sort_values('count', ascending=False)

    # print(sorted_by_count.head(5))

    data = data.loc[data["disease_id"] == "C0002395"]

    dataB = pd.read_excel("proteinas_en_comun_Alzheimer.xlsx")

    data = data[~((data["disease_id"] == "C0002395") &
                  (data["protein_id"].isin(dataB["protein_id"])))]"""
    sequences = data["protein_sequence"]

    return sequences

def similitudProteinas(sequences):
    output = []

    for row1 in sequences:
        for row2 in sequences:
            if row1 != row2:
                #print(row1,row2)
                similarity = Levenshtein.distance(row1, row2) / max(len(row1), len(row2))
                output.append([row1, row2, similarity*100])
    return output

def metrica_distanciaProteinas():
    # Leer los archivos CSV
    data = pd.read_csv("C0002395_Disease_patronesIdenticos_ocurrence20%_SinDescarte.csv")
    df_b = pd.read_csv("C0002395Disease_%Similitud.csv")

    # Crear un diccionario de similaridades
    proteinas_dict = dict(zip(zip(df_b['Proteina1'], df_b['Proteina2']), df_b['Similaridad']))

    # Agrupar por el patrón de proteína
    grupos = data.groupby('Patron')

    # Crear una lista de tuplas con los índices únicos de las filas en cada grupo
    indices = [grupo.index for patron, grupo in grupos]

    # Generar todas las combinaciones únicas de índices
    index_combinations = [(i, j) for grp in indices for i in grp for j in grp if i != j]

    # Filtrar las combinaciones que no son duplicadas y tienen diferencias en las filas correspondientes
    filtered_combinations = [comb for comb in index_combinations if not data.loc[comb[0]].equals(data.loc[comb[1]])]

    # Filtrar las combinaciones que existen en el diccionario de similaridades
    output = [[data.loc[comb[0], 'Patron'], data.loc[comb[0], 'Proteina'], data.loc[comb[1], 'Proteina'],
               proteinas_dict.get((data.loc[comb[0], 'Proteina'], data.loc[comb[1], 'Proteina']), '')] for comb in
              filtered_combinations]

    # Crear un DataFrame a partir de la lista de resultados
    df = pd.DataFrame(output, columns=['Patron', 'Proteina1', 'Proteina2', 'Similitud'])

    # Guardar el DataFrame en un archivo CSV
    df.to_csv('Metrica_distanciaProteinasMismoPatron_C0002395_Disease_patronesIdenticos_ocurrence20%_SinDescarte.csv',
              index=False)

    """data = pd.read_csv("C0002395_Disease_patronesIdenticos_ocurrence20%_SinDescarte.csv").head(100)

    df_b = pd.read_csv("C0002395Disease_%Similitud.csv")

    proteinas_dict = df_b.set_index(['Proteina1', 'Proteina2'])['Similaridad'].to_dict()

    # Agrupar por el patrón de proteína
    output = []
    grupos = data.groupby('Patron')

    # Iterar sobre los grupos
    for patron, grupo in grupos:
        for index1, row1 in grupo.iterrows():
            for index2, row2 in grupo.iterrows():
                if index1 != index2 and not (row1 == row2).all():
                    if (row1['Proteina'], row2["Proteina"]) in proteinas_dict:
                        output.append([patron, row1["Proteina"], row2["Proteina"],
                                       proteinas_dict[(row1['Proteina'], row2['Proteina'])]])

    df = pd.DataFrame(output, columns=['Patron', 'Proteina1', 'Proteina2', 'Similitud'])
    df.to_csv('Metrica_distanciaProteinasMismoPatron_C0002395_Disease_patronesIdenticos_ocurrence20%_SinDescarte.csv',
              index=False)"""

def patronesComun():
    # Leer el archivo CSV y cargar los datos en una lista de diccionarios
    registros = []
    with open("C0002395_Disease_patronesIdenticos_ocurrence10%_SinDescarte.csv", 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            registros.append(row)

    # Diccionario para almacenar la cantidad de patrones únicos por proteína
    patrones_por_proteina = {}

    # Iterar sobre los registros y extraer los patrones únicos de cada proteína
    for registro in registros:
        proteina = registro['Proteina']
        patron = registro['Patron']
        if proteina not in patrones_por_proteina:
            patrones_por_proteina[proteina] = set()
        patrones_por_proteina[proteina].add(patron)

    # Diccionario para almacenar las proteinas que tienen en común cada par de proteinas
    proteinas_comunes = {}
    pares_proteinas_procesados = set()
    # Filtrar las proteínas que tienen al menos 10 patrones únicos en común
    for proteina1, patrones1 in patrones_por_proteina.items():
        for proteina2, patrones2 in patrones_por_proteina.items():
            if proteina1 != proteina2 and (proteina2, proteina1) not in pares_proteinas_procesados:
                patrones_comunes = patrones1.intersection(patrones2)
                if len(patrones_comunes) >= 96:
                    par_proteinas = (proteina1, proteina2)
                    proteinas_comunes[par_proteinas] = patrones_comunes
                    pares_proteinas_procesados.add(par_proteinas)
    output = []

    df_b = pd.read_csv("C0002395Disease_%Similitud.csv")

    proteinas_dict = df_b.set_index(['Proteina1', 'Proteina2'])['Similaridad'].to_dict()

    for par_proteinas, patrones_comunes in proteinas_comunes.items():
        proteina1, proteina2 = par_proteinas
        pattern_lengths = {}
        for pattern in patrones_comunes:
            length = len(pattern)
            key = f'Longitud {length}'
            if key in pattern_lengths:
                pattern_lengths[key] += 1
            else:
                pattern_lengths[key] = 1
        sorted_pattern_lengths = dict(sorted(pattern_lengths.items(), key=lambda x: int(x[0][9:]), reverse=True))
        if proteina1 != proteina2:
            if (proteina1, proteina2) in proteinas_dict:
                output.append([sorted_pattern_lengths, proteina1, proteina2,
                               proteinas_dict[(proteina1, proteina2)]])

    df = pd.DataFrame(output, columns=['Patrones', 'Proteina1', 'Proteina2', 'Similitud'])
    df.to_csv('Metrica_patronesComunes10%_C0002395_Disease_patronesIdenticos_ocurrence10%_SinDescarte.csv',
              index=False)

def remplazar_sequence_for_ID(output):
    df_b = pd.read_excel("data_nervous_genes.xlsx")

    # Ordenar de mayor a menor tamaño. Las subcadenas del mismo tamaño se ordenan por orden alfabetico
    output_ordered = sorted(output, key=lambda x: (-len(x[0]), x[0]))

    proteinas_dict = dict(df_b[['protein_sequence', 'protein_id']].values)

    for item in output_ordered:
        protein_sequence1 = item[0]
        protein_sequence2 = item[1]
        if protein_sequence1 in proteinas_dict and protein_sequence2 in proteinas_dict:
            item[0] = proteinas_dict[protein_sequence1]
            item[1] = proteinas_dict[protein_sequence2]



    df_a = pd.DataFrame(output_ordered, columns=['Proteina1', 'Proteina2', 'Similaridad'])

    # Guardar el DataFrame actualizado en un archivo CSV
    df_a.to_csv('AllProteins_%Similitud.csv', index=False)


if __name__ == "__main__":
    inicio = time.time()
    pattern_freqMin = dict()

    #sequences = readData()
    #similitud = similitudProteinas(sequences)
    #remplazar_sequence_for_ID(similitud)
    patronesComun()
    """output = patronesComun(data)
    remplazar_sequence_for_ID(output)"""

    fin = time.time()

    tiempo_total = fin - inicio
    print(tiempo_total, "segundos")