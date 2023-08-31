import pandas as pd
import Levenshtein

def readData(archivoEntrada):
    data = pd.read_excel(archivoEntrada).head(10)

    sequences = data["protein_sequence"]

    return sequences

def similitudProteinas(sequences):
    output = []

    for row1 in sequences:
        for row2 in sequences:
            if row1 != row2:
                similarity = Levenshtein.distance(row1, row2) / max(len(row1), len(row2))
                output.append([row1, row2, similarity*100])
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
    df_a.to_csv('AllProteins_%Similitud.csv', index=False)

if __name__ == "__main__":
    archivoEntrada = "data_nervous_genes.xlsx"
    sequences = readData(archivoEntrada)


    output = similitudProteinas(sequences)
    remplazar_sequence_for_ID(output)
