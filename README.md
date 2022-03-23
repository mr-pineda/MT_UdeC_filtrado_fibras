# Filtrado de Fibras

## Descripción:
El presente repositorio contiene el algoritmo utilizado al momento de la presentación de la memoria de título **"OPTIMIZACIÓN EN ALGORITMO DE SEGMENTACIÓN AUTOMÁTICA DE FASCÍCULOS DE MATERIA BLANCA SUPERFICIAL EN DATOS DE TRACTOGRAFÍAS PRBABILÍSTICAS"**, desarrollado por **Felipe Pineda S.** para optar al grado de ingeniero civíl electrónico en la Universidad de Concepción.

En este se desarrolla una implementación en C del algoritmo de filtrado desarrollado en Python por **Cristobal Mendoza** en la memomria de título "**...ANOTAR_NOMBRE_MEMORIA_CRISTOBAL...**". El objetivo principal de esta implementación es disminuir el tiempo de ejecución en el algoritmo de filtrado de fibras ruidosas usando la distancia SSPD como métrica comparativa. El objetivo secundario es desarrollar una librería que simplifique el trabajo con fibras en C para futuros trabajos en este ámbito.

## Contenido:
El principal contenido son las librearías desarrolladas en C (1), los archivos main.c para las implementaciones de ISOMAP y SSPD_3D y archivos desarrollados en MATLAB destinados a realizar pruebas de visualización de los fascículos.

Los otros elementos son generados por el IDE utilizado ([Code::Blocks](https://www.codeblocks.org/)) y están relacionados con la configuración del entorno de programación para el equipo utilizado. Para la correcta compilación de los scripts, se recomienda eliminar estos archivos y configurar el entorno de programación según el equipo utilizado. 

- **Librerias en C (1)**: Separadas en sus archivos de encabezado (.h) y definición (.c) ubicados en el directorio "headers_sources".
   - _**"directory.h"**: Librería dedicada a facilitar la navegación entre archivos y directorios (Requiere librería "dirent.h")._
   - _**"moreMath.h"**: Contiene funciones matemáticas adicionales para fascílitar la obtención métricas básicas (como la distancia euclídea en 2 y 3 dimensiones) y la operación entre matrices._
   - _**"fiberlib.h"**: Librería dedicada al trabajo con fibras. Contiene las estructuras 'fiber' y 'fasciculus' las cuales respresentan a las fibras individuales y a un conjunto de fibras respectivamente. Ademas contiene funciones para la lectura y escritura de archivos relacionados con fibras, métricas de fibras entre otras (Requiere el uso de la librería "mormeMath.h")._
   - _**"sspdlib.h"**: Contiene las funciones utilizadas para el cómputo de la distancia SSPD (Requiere el uso de las ibrerías "moreMath.h" y "fiberlib.h")._
   - _**"isomap.h"**: Contiene las funciones utilizadas para aplicar los algoritmos de ISOMAP y MDS cláisco **(2)**._

- **Archivos main.c**: Las implementaciones de los algoritmos de filtrado por ISOMAP y SSPD_3D estan ubicados en sus respectivos directorios. De forma adicional al filtrado se añadió al script una sección que mide el tiempo de filtrado por fascículo. En el caso del filtrado por ISOMAP, se muestra el tiempo empleado por la reduccion dimensional (ISOMAP o MDS) y el tiempo empleado total (reducción dimensional y filtrado) por fascículo.

- **Archivos de MATLAB (.m)**: Contiene una implmenetación simplificada del algoritmo de MDS y un visualizador de fibras en espacio tridimensional y bidimensional. Este algoritmo fue desarrollado con el fin de corroborar los resultados obtenidos mediante la implementación en C

## Observaciones:

 **(1)** Además de las librerías desarrolladas, se reaquiere añadir las librerías "strings.h", "limits.h", "math.h", "stdbool.h" y "time.h". En el desarrollo de este algoritmo se utilizó el compilador [MinGW-w64](https://www.mingw-w64.org/), el cual estaba incluído en el IDE, permitiendo el uso de las librearias "dirent.h" para el trabajo con directorios y "omp.h" para el uso de programación paralela.
    
 **(2)** En el algoritmo de MDS clásico se utiliza ala obtención de los valores propios de una matriz a travéz del método de Jacobi. Para esto se utilizó la implementación desarrollada por **John Burkardt** disponible en https://people.sc.fsu.edu/~jburkardt/c_src/jacobi_eigenvalue/jacobi_eigenvalue.html . Este algoritmo esta bajo la licencia [GNU LGPL](https://www.gnu.org/licenses/lgpl-3.0.en.html).
