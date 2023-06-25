import time
import pandas as pd
import Levenshtein

def readData():
    data = pd.read_csv("C0002395_Disease_patronesIdenticos_ocurrence20%_ConProteina.csv").head(2000)

    return data


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

def proteinasSimilares(data):
    # Agrupar por el patrón de proteína
    output = []
    grupos = data.groupby('Patron')

    # Iterar sobre los grupos
    for patron, grupo in grupos:
        for index1, row1 in grupo.iterrows():
            for index2, row2 in grupo.iterrows():
                if index1 != index2 and not (row1 == row2).all():
                    lista_posiciones1 = eval(row1['Posiciones'])
                    lista_posiciones2 = eval(row2['Posiciones'])
                    if (len(lista_posiciones1) >= 0.02 * len(row1['Proteina']) and len(lista_posiciones2) >= 0.02 * len(row2['Proteina'])):
                        similarity = Levenshtein.distance(row1['Proteina'], row2['Proteina']) / max(len(row1['Proteina']),
                                                                                                 len(row2['Proteina']))
                        output.append([patron, row1["Proteina"], row2["Proteina"], similarity*100])

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
    df_a.to_csv('ProteinasSimilares10%_C0002395_20%.csv', index=False)


if __name__ == "__main__":
    inicio = time.time()
    pattern_freqMin = dict()

    data = readData()

    output = proteinasSimilares(data)
    remplazar_sequence_for_ID(output)

    fin = time.time()

    tiempo_total = fin - inicio
    print(tiempo_total, "segundos")