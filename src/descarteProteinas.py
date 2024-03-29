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

def readData(archivoEntrada, enfermedad):
    data = pd.read_excel(archivoEntrada)

    if (enfermedad != ''):
        data = data.loc[data["disease_id"] == enfermedad]

        dataB = pd.read_excel("proteinas_en_comun_Alzheimer.xlsx")

        data = data[~((data["disease_id"] == enfermedad) &
                      (data["protein_id"].isin(dataB["protein_id"])) &
                      (data["gene_id"].isin(dataB["gene_id"])))]

    sequences = data["protein_sequence"]

    return sequences

def levenshtein_similarity(pattern1, pattern2):
    return Levenshtein.distance(pattern1, pattern2) / max(len(pattern1), len(pattern2))

def descarte(data, threshold):
    # Datos de ejemplo
    data = data.tolist()

    # Crear matriz de similitud
    num_points = len(data)
    similarity_matrix = [[0] * num_points for _ in range(num_points)]
    for i in range(num_points):
        for j in range(num_points):
            similarity_matrix[i][j] = levenshtein_similarity(data[i], data[j])

    # Parámetros del algoritmo DBSCAN
    eps = 0.15  # Umbral de similitud
    min_samples = 2  # Número mínimo de muestras para formar un cluster

    # Ejecutar el algoritmo DBSCAN
    dbscan_instance = dbscan(similarity_matrix, eps, min_samples, metric='precomputed')
    dbscan_instance.process()
    clusters = dbscan_instance.get_clusters()

    filtered_clusters = []
    discarded_data = []

    for cluster_id, cluster in enumerate(clusters):
        filtered_cluster = []
        min_avg_distance = float('inf')
        central_point_index = None

        # Calcular la distancia promedio para cada punto del cluster
        for point_index in cluster:
            total_distance = 0
            for other_index in cluster:
                total_distance += similarity_matrix[point_index][other_index]
            avg_distance = total_distance / len(cluster)
            if avg_distance < min_avg_distance:
                min_avg_distance = avg_distance
                central_point_index = point_index

        # Verificar si el punto central supera el umbral
        similarity_percentage = 1 - (min_avg_distance / eps)
        if similarity_percentage >= threshold:
            filtered_cluster.append(central_point_index)
        else:
            discarded_data.extend([data[i] for i in cluster if i != central_point_index])

        if filtered_cluster:
            filtered_clusters.append(filtered_cluster)

    data = remplazar_sequence_for_ID(data)

    # Imprimir los resultados
    for cluster_id, cluster in enumerate(filtered_clusters):
        cluster_data = [data[i] for i in cluster]
        print(f'Cluster {cluster_id}: {", ".join(cluster_data)}')

    discarded_data = remplazar_sequence_for_ID(discarded_data)
    # Guardar los datos descartados en un archivo CSV utilizando Pandas
    if discarded_data:
        df = pd.DataFrame({'ProteinasDescartadas': discarded_data})
        df.to_csv('resultados/proteinasDescartadas.csv', index=False)


def remplazar_sequence_for_ID(output):
    df_b = pd.read_excel("data_nervous_genes.xlsx")

    proteinas_dict = dict(df_b[['protein_sequence', 'protein_id']].values)

    for i in range(len(output)):
        protein_sequence = output[i]
        if protein_sequence in proteinas_dict:
            output[i] = proteinas_dict[protein_sequence]

    return output

def ejecutar(archivoEntrada, enfermedad, similitud):
    data = readData(archivoEntrada, enfermedad)
    descarte(data, similitud)

