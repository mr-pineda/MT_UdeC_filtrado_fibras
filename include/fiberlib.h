#ifndef FIBERS_H_INCLUDED
#define FIBERS_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <math.h>
#include <stdbool.h>
#include <unistd.h>
#include <limits.h>
#include <time.h>
#include <errno.h>

typedef struct {
    // Estructura que se usa para definir una fibra:
    //     numPts: Entero que indica la cantidad de pun tos que tiene la fibra.
    //     coordinates: Arreglo din�mico de floats, generalmente de tama�o numPts x 3 que almacena la coordenadas de la fibra.
    //                  Entonces coordinates[n] corresponde al n-�simo punto de la fibra
    int numPts;
    float** cordinates;
}fiber;

typedef struct {
    // Estructura dque se usa para definir un conjunto de fibras.
    //     name: string que almacena el nombre del conjunto (o nombre del bundle).
    //     numFibers: entero que almacena la cantidad de fibras en el conjuntos.
    //     fibers: Arreglo dinamico de fibras, de tama�o numFibers. Al trabajar con arreglos din�micos se puede hacer un bundle
    //             con una cantidad variable de fibras, que tengan una cantidad variable de puntos.
    char name[60];
    int numFibers;
    fiber* fibers;
}fasciculus;

fasciculus readBundle(char* fileName,char* path); //Lee los archivos .bundles y .bundlesdata y guarda su informaci�n en una variable fasciculus
void writeBundles(fasciculus bundle, char* fileName,char* path); //Escribe los datos de variable fasciculues en archivos .bundles y .bundlesdata
void clearBundle(fasciculus bundle); //Libera la memoria utilizada en un dato del tipo fasciculus
fasciculus selectFibers(fasciculus bundle, int* selected); //se seleccionas las fibras de un fasciculus a partir de un vector de indices
float fiberlen(fiber fib); //Se entrega el largo de la fibra ingresada
void fiberSwap(fiber* f1, fiber* f2, int dim);// Se intercambian los valores de 2 fibras

#endif // FIBERS_H_INCLUDED