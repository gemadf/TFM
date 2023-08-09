import tkinter as tk


def interfaz():
    datos = dict()

    ventana = tk.Tk()
    ventana.geometry("800x500")

    ventana.grid_rowconfigure(0, weight=1)
    ventana.grid_rowconfigure(1, weight=1)
    ventana.grid_rowconfigure(2, weight=1)
    ventana.grid_rowconfigure(3, weight=1)
    ventana.grid_columnconfigure(0, weight=1)
    ventana.grid_columnconfigure(1, weight=1)

    titulo = tk.Label(ventana, text="DECUBRIMIENTO DE PATRONES SIMBOLICOS:")
    titulo.grid(row=0, column=0, columnspan=2, sticky="nsew", padx=10, pady=10)

    label1 = tk.Label(ventana, text="Nombre del archivo con los datos de entrada*:", anchor="w")
    label1.grid(row=1, column=0, sticky="nsew", padx=10, pady=10)

    input_frame1 = tk.Frame(ventana)
    input_frame1.grid(row=1, column=1, sticky="nsew", padx=10, pady=10)
    input1 = tk.Entry(input_frame1)
    input1.pack(fill="both", padx=5, pady=5)

    label2 = tk.Label(ventana, text="Id de la enfermedad para la cual se deseean obtener los patrones:", anchor="w")
    label2.grid(row=2, column=0, sticky="nsew", padx=10, pady=10)

    input_frame2 = tk.Frame(ventana)
    input_frame2.grid(row=2, column=1, sticky="nsew", padx=10, pady=10)
    input2 = tk.Entry(input_frame2)
    input2.pack(fill="both", padx=5, pady=5)

    label3 = tk.Label(ventana, text="Nombre del archivo con las proteínas target de la enfermedad:", anchor="w")
    label3.grid(row=3, column=0, sticky="nsew", padx=10, pady=10)

    input_frame3 = tk.Frame(ventana)
    input_frame3.grid(row=3, column=1, sticky="nsew", padx=10, pady=10)
    input3 = tk.Entry(input_frame3)
    input3.pack(fill="both", padx=5, pady=5)

    label4 = tk.Label(ventana, text="% de ocurrencia mínima (0-1)*:", anchor="w")
    label4.grid(row=4, column=0, sticky="nsew", padx=10, pady=10)

    input_frame4 = tk.Frame(ventana)
    input_frame4.grid(row=4, column=1, sticky="nsew", padx=10, pady=10)
    input4 = tk.Entry(input_frame4)
    input4.pack(fill="both", padx=5, pady=5)

    label5 = tk.Label(ventana, text="% de similitud entre dos patrones para que se realice descarte (0-1)*:", anchor="w")
    label5.grid(row=5, column=0, sticky="nsew", padx=10, pady=10)

    input_frame5 = tk.Frame(ventana)
    input_frame5.grid(row=5, column=1, sticky="nsew", padx=10, pady=10)
    input5 = tk.Entry(input_frame5)
    input5.pack(fill="both", padx=5, pady=5)

    boton = tk.Button(ventana, text="Obtener datos", command=lambda: obtener_datos(ventana, datos, input1, input2, input3, input4, input5))
    boton.grid(row=6, column=0, columnspan=2, padx=10, pady=10)

    ventana.mainloop()

    return datos

def obtener_datos(ventana, datos, input1, input2, input3, input4, input5):

    datos["NombreArchivoEntrada"] = input1.get()
    datos["CodigoEnfermedad"] = input2.get()
    datos["NombreArchivoTarget"] = input3.get()
    datos["OcurrenciaMin"] = input4.get()
    datos["Similitud"] = input5.get()

    ventana.destroy()

if __name__ == "__main__":
    datos = interfaz()

