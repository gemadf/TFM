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
    sequences = data["protein_sequence"]

    return sequences



def similarity_percentage(pattern1, pattern2):
    return Levenshtein.distance(pattern1, pattern2) / max(len(pattern1), len(pattern2))

def descarte(data):

    data = data.tolist()
    # Crear una matriz de similitud utilizando el porcentaje de similitud
    num_samples = len(data)
    similarity_matrix = np.zeros((num_samples, num_samples))
    for i in range(num_samples):
        for j in range(num_samples):
            similarity_matrix[i, j] = similarity_percentage(data[i], data[j])

    # Definir los parámetros del algoritmo DBSCAN
    epsilon = 0.05  # Radio de vecindad (ajustar según tus necesidades)
    min_samples = 2  # Número mínimo de puntos para formar un cluster

    # Crear una instancia del algoritmo DBSCAN con la métrica de similitud
    dbscan = DBSCAN(eps=epsilon, min_samples=min_samples, metric='precomputed')

    # Ejecutar el algoritmo DBSCAN en la matriz de similitud
    labels = dbscan.fit_predict(similarity_matrix)

    data = remplazar_sequence_for_ID(data)

    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)

    print(f"Número de clusters: {n_clusters}")

    # Crear un diccionario para almacenar los clusters
    clusters = {}
    for i, label in enumerate(labels):
        if label not in clusters:
            clusters[label] = []
        clusters[label].append(i)

    # Imprimir los resultados del clustering
    for label, indices in clusters.items():
        print(f"Cluster: {label}")
        for i in indices:
            print(f"Datos: {data[i]}")

        # Filtrar las proteínas por similitud > 90%
        filtered_indices = [i for i in indices if similarity_percentage(data[i], data[i]) > 0.9]

        # Imprimir solo las proteínas que superan el umbral de similitud
        print("Proteínas con similitud > 90%:")
        for i in filtered_indices:
            print(f"Datos: {data[i]}")

def remplazar_sequence_for_ID(output):
    df_b = pd.read_excel("data_nervous_genes.xlsx")

    proteinas_dict = dict(df_b[['protein_sequence', 'protein_id']].values)

    for i in range(len(output)):
        protein_sequence = output[i]
        if protein_sequence in proteinas_dict:
            output[i] = proteinas_dict[protein_sequence]

    return output

if __name__ == "__main__":
    inicio = time.time()
    pattern_freqMin = dict()
    threshold = 90

    data = readData()
    #df = remplazar_sequence_for_ID(output)
    #central_data = descarte(df, threshold)
    central_data = descarte(data)

    fin = time.time()

    tiempo_total = fin - inicio
    print(tiempo_total, "segundos")