import pandas as pd
import Levenshtein
import time
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
import numpy as np
from scipy.spatial.distance import pdist, squareform
from pyclustering.cluster.dbscan import dbscan
from pyclustering.utils import timedcall
from Levenshtein import distance

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

def levenshtein_similarity(pattern1, pattern2):
    return Levenshtein.distance(pattern1, pattern2) / max(len(pattern1), len(pattern2))

def descarte(data):
    # Datos de ejemplo
    data = data.tolist()

    # Crear matriz de similitud
    num_points = len(data)
    similarity_matrix = [[0] * num_points for _ in range(num_points)]
    for i in range(num_points):
        for j in range(num_points):
            similarity_matrix[i][j] = levenshtein_similarity(data[i], data[j])

    # Parámetros del algoritmo DBSCAN
    eps = 0.02  # Umbral de similitud
    min_samples = 2  # Número mínimo de muestras para formar un cluster

    # Ejecutar el algoritmo DBSCAN
    dbscan_instance = dbscan(similarity_matrix, eps, min_samples, metric='precomputed')
    dbscan_instance.process()
    clusters = dbscan_instance.get_clusters()

    # Aplicar umbral de similitud del 90%
    threshold = 0.9
    filtered_clusters = []
    for cluster_id, cluster in enumerate(clusters):
        filtered_cluster = []
        for point_index in cluster:
            similarity_percentage = 1 - (similarity_matrix[point_index][point_index] / eps)
            if similarity_percentage >= threshold:
                filtered_cluster.append(point_index)
        if filtered_cluster:
            filtered_clusters.append(filtered_cluster)

    data = remplazar_sequence_for_ID(data)

    # Imprimir los resultados
    for cluster_id, cluster in enumerate(filtered_clusters):
        cluster_data = [data[i] for i in cluster]
        print(f'Cluster {cluster_id}: {", ".join(cluster_data)}')


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
    threshold = 0.9

    data = readData()
    descarte(data)

    fin = time.time()

    tiempo_total = fin - inicio
    print(tiempo_total, "segundos")