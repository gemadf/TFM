import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

if __name__ == "__main__":
    # Cargar los datos desde el archivo CSV
    data = pd.read_csv('Metrica_distanciaProteinasMismoPatron_C0002395_Disease_patronesIdenticos_ocurrence20%_SinDescarte.csv')

    plt.hist(data['Similitud'], bins=100)

    plt.xlabel('Valor Numérico')
    plt.ylabel('Frecuencia')
    plt.title('Distribución de Valores Numéricos')
    #plt.xlim(70, 100)

    plt.savefig("FiguraRangoCompleto_Metrica_distanciaProteinasMismoPatron_C0002395_Disease_patronesIdenticos_ocurrence20%_SinDescarte")
    plt.show()

