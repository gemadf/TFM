import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def grafica(archivo, nombreOutput):
    # Cargar los datos desde el archivo CSV
    data = pd.read_csv(archivo)

    plt.hist(data['Similitud'], bins=100)

    plt.xlabel('Valor Numérico')
    plt.ylabel('Frecuencia')
    plt.title('Distribución de Valores Numéricos')

    plt.savefig(nombreOutput)

