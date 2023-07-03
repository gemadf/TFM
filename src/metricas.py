import time
import pandas as pd
import Levenshtein
import csv

def readData():
    data2 = pd.read_csv("C0002395_Disease_patronesIdenticos_ocurrence20%_ConProteina.csv").head(3)

    data = pd.read_excel("data_nervous_genes.xlsx")
    count_by_disease = data.groupby('disease_id').size().reset_index(name='count')
    sorted_by_count = count_by_disease.sort_values('count', ascending=False)

    # print(sorted_by_count.head(5))

    data = data.loc[data["disease_id"] == "C0002395"]

    dataB = pd.read_excel("proteinas_en_comun_Alzheimer.xlsx")

    data = data[~((data["disease_id"] == "C0002395") &
                  (data["protein_id"].isin(dataB["protein_id"])) &
                  (data["gene_id"].isin(dataB["gene_id"])))]
    sequences = data["protein_sequence"]

    return data2, sequences

def similitudProteinas(sequences):
    output = []

    for row1 in sequences:
        for row2 in sequences:
            if row1 != row2:
                similarity = Levenshtein.distance(row1, row2) / max(len(row1), len(row2))
                output.append([row1, row2, similarity*100])
    return output

def distanciaProteinas(data):
    # Agrupar por el patrón de proteína
    output = []
    grupos = data.groupby('Patron')

    # Iterar sobre los grupos
    for patron, grupo in grupos:
        for index1, row1 in grupo.iterrows():
            for index2, row2 in grupo.iterrows():
                if index1 != index2 and not (row1 == row2).all():
                    similarity = Levenshtein.distance(row1['Proteina'], row2['Proteina']) / max(len(row1['Proteina']), len(row2['Proteina']))
                    output.append([patron, row1["Proteina"], row2["Proteina"], similarity*100])

    return output

def patronesComun(data):
    # Leer el archivo CSV y cargar los datos en una lista de diccionarios
    registros = []
    with open("C0002395_Disease_patronesIdenticos_ocurrence20%_ConProteina.csv", 'r') as file:
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
                if len(patrones_comunes) >= 10:
                    par_proteinas = (proteina1, proteina2)
                    proteinas_comunes[par_proteinas] = patrones_comunes
                    pares_proteinas_procesados.add(par_proteinas)
    output = []

    for par_proteinas, patrones_comunes in proteinas_comunes.items():
        proteina1, proteina2 = par_proteinas
        if proteina1 != proteina2:
            similarity = Levenshtein.distance(proteina1, proteina2) / max(len(proteina1), len(proteina2))
            output.append([patrones_comunes, proteina1, proteina2, similarity*100])

    return output

def remplazar_sequence_for_ID(output):
    df_b = pd.read_excel("data_nervous_genes.xlsx")

    # Ordenar de mayor a menor tamaño. Las subcadenas del mismo tamaño se ordenan por orden alfabetico
    output_ordered = sorted(output, key=lambda x: (-len(x[0]), x[0]))

    proteinas_dict = dict(df_b[['protein_sequence', 'protein_id']].values)

    for item in output_ordered:
        protein_sequence1 = item[1]
        protein_sequence2 = item[2]
        if protein_sequence1 in proteinas_dict and protein_sequence2 in proteinas_dict:
            item[1] = proteinas_dict[protein_sequence1]
            item[2] = proteinas_dict[protein_sequence2]



    df_a = pd.DataFrame(output_ordered, columns=['Patron', 'Proteina1', 'Proteina2', 'Similaridad'])

    # Guardar el DataFrame actualizado en un archivo CSV
    df_a.to_csv('aa_C0002395_20%.csv', index=False)


if __name__ == "__main__":
    inicio = time.time()
    pattern_freqMin = dict()

    data, sequences = readData()
    output = similitudProteinas(sequences)
    remplazar_sequence_for_ID(output)
    #output = distanciaProteinas(data)
    """output = patronesComun(data)
    remplazar_sequence_for_ID(output)"""

    fin = time.time()

    tiempo_total = fin - inicio
    print(tiempo_total, "segundos")