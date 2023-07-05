import pandas as pd
import Levenshtein
import time
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
import numpy as np

def readData():
    data = pd.read_excel("data_nervous_genes.xlsx")
    count_by_disease = data.groupby('disease_id').size().reset_index(name='count')
    sorted_by_count = count_by_disease.sort_values('count', ascending=False)

    # print(sorted_by_count.head(5))

    data = data.loc[data["disease_id"] == "C0002395"]

    dataB = pd.read_excel("proteinas_en_comun_Alzheimer.xlsx")

    data = data[~((data["disease_id"] == "C0002395") &
                  (data["protein_id"].isin(dataB["protein_id"])) &
                  (data["gene_id"].isin(dataB["gene_id"])))]
    sequences = data["protein_sequence"].head(200)

    return sequences

def similitudProteinas(sequences):
    output = []
    processed_combinations = set()

    for row1 in sequences:
        for row2 in sequences:
            if row1 != row2:
                combination = frozenset([row1, row2])
                if combination not in processed_combinations:
                    similarity = Levenshtein.distance(row1, row2) / max(len(row1), len(row2))
                    output.append([row1, row2, similarity*100])
                    processed_combinations.add(combination)
    return output

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
    df_a.to_csv('DistanciaProteinas_C0002395_20%.csv', index=False)

    return df_a


def descarte(data, threshold):
    filtered_data = data[data['Similaridad'] >= threshold]

    X = filtered_data['Similaridad'].values.reshape(-1, 1)

    # Normaliza los datos
    scaler = StandardScaler()
    X = scaler.fit_transform(X)

    # Realiza el clustering de densidad con DBSCAN
    dbscan = DBSCAN(eps=0.01, min_samples=2)  # Ajusta los parámetros según tus necesidades
    labels = dbscan.fit_predict(X)

    # Agrega las etiquetas de los clusters al DataFrame original
    filtered_data.loc[:, 'Cluster'] = labels

    filtered_data = filtered_data.reset_index(drop=True)

    central_proteins = []
    for cluster_label in set(labels):
        cluster_indices = np.where(labels == cluster_label)[0]
        central_index = cluster_indices[np.argmax(X[cluster_indices])]
        central_protein = filtered_data.loc[central_index, 'Proteina1']
        central_proteins.append(central_protein)

    # Filtrar el DataFrame original para obtener los elementos centrales de cada cluster
    central_data = data[data['Proteina1'].isin(central_proteins)]

    # Filtrar el DataFrame original para obtener las filas correspondientes al cluster -1
    cluster_minus1 = filtered_data[filtered_data['Cluster'] == -1]

    # Filtrar el DataFrame original para obtener las filas con similitud inferior al umbral
    below_threshold = data[data['Similaridad'] < threshold]

    # Concatenar las filas del cluster -1 y las filas con similitud inferior al umbral
    final_data = pd.concat([central_data, cluster_minus1, below_threshold])
    final_data = final_data.sort_index()
    final_data = final_data.drop_duplicates()

    # Imprimir los resultados
    print("Datos finales:")
    print(final_data)

    #central data tiene las proteinas centrales y los valoreles atípicos que no se encuandran en ningun cluster
    return central_data

def concatenarResultadosCluster_datosOriginales(central_data, df, threshold):
    df = df[df['Similaridad'] < threshold]


if __name__ == "__main__":
    inicio = time.time()
    pattern_freqMin = dict()
    threshold = 90

    data = readData()
    output = similitudProteinas(data)
    df = remplazar_sequence_for_ID(output)
    central_data = descarte(df, threshold)
    concatenarResultadosCluster_datosOriginales(central_data, df, threshold)

    fin = time.time()

    tiempo_total = fin - inicio
    print(tiempo_total, "segundos")