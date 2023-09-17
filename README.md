# Búsca de patrones identicos para secuencias genéticas

## Instrucciones de ejecución

### Ejecutar el script "similitudAllProteins.py"

Debemos ejecutar el script para obtener el porcentaje de similitud entre todas las combinaciones de proteinas del conjunto. En el script, debemos modificar la ruta del archivo para el cual queremos obtener los resultados. 

Una vez tengamos los resultados, debemos mover el archivo a la misma ruta desde donde ejecutemos el script, ya que el script final necesita ese archivo para obtener los resultados de las métricas.

### Ejecutar el script patrones_identicos.py

Al ejecutar el script, nos aparecerá una interfaz en la cual deberemos introducir la información especificada. Los datos obligatorios se marcan con un *.

Como salida de la ejecución obtendremos:
    - Conjunto de proteínas descartadas.
    - Conjunto de patrones identicos encontrados con sus proteínas y posiciones asociadas.
    - Resultados de las métricas:
        - Patron, proteínas y similitud de las proteínas analizadas.
        - Figura con la distribución de las similitudes encontradas trás aplicar cada métrica