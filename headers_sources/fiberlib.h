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

typedef unsigned int u_int;
typedef unsigned short u_short;
typedef unsigned long u_long;
typedef unsigned long long u_long2;
typedef long long long2;
typedef unsigned short byte;

typedef float f_data;

/**
 * @brief Representation of an individual fiber.
 *
 * An struct thar represent an individual fiber in M-dimetion (An array of M-dimentional points).
 * Three-dimensional points are usually used, but certain applications used in this project require
 * two-dimensional points. In case it is required to use for other not fibers related applications,
 * the possibility of using four-dimensional points or more is left open, but it has not been studied
 * in depth what possible issues can be found.
 */
typedef struct
{
    int numPts;          /**< Number of M-dimentional points.*/
    f_data **cordinates; /**< Array of consecutive M-dimentional points coordinates.*/
} fiber;

/**
 * @brief Representation of a bundle of fibers (as array)
 *
 * An struct that represent a bundle of fibers. Each fiber can contain a variable number of points.
 */
typedef struct
{
    char name[60]; /**< Name of the bundle.*/
    int numFibers; /**< Number of fibers in the bundle*/
    fiber *fibers; /**< Array of fibers.*/
} fasciculus;

/**
 * @brief Representation of a bundle of fibers (as doubly linked list)
 *
 * An struct that represent a bundle of fibers, but defined has doubly linked list. In some aplications
 * is more useful use a linked list than an array.
 */
typedef struct
{
    int numPts;          /**< Number of M-dimentional points.*/
    f_data **cordinates; /**< Array of consecutive M-dimentional points coordinates.*/
    fiber_list *next;    /**< Pointer to the next fiber.*/
    fiber_list *prev;    /**< Pointer to the previous fiber.*/
} fiber_list;

/**
 * @brief Extract the bundles data from files.
 * Open a pair of .bundles .bundlesdata files, extract the data and store it in a fasciculus structure.
 * @param fileName The name of the file. (Must have the .bundles extention).
 * @param path The path where both files are alocated (Relative path can be used).
 * @return fasciculus (Bundle of fiber as array).
 */
fasciculus readBundle(char *fileName, char *path);

/**
 * @brief Write the data of a fasciculus into files.
 * Write the fasciculus data (name, number of points and coordinates) into a pair of
 * .bundles .bundlesdata files.
 * @param fileName The name of the file. (Must have the .bundles extention).
 * @param path The path where both files are alocated (Relative path can be used).
 * @param path
 */
void writeBundles(fasciculus bundle, char *fileName, char *path);

/**
 * @brief Release memory use by fasciculus.
 * Release the memory allocation used by the fasciculus and all the fibers contaided
 * in it. The fibers pointer will get NULL adress and the numPts data will be zero;
 * @param bundle fasciculus data that will be deleted-
 */
void clearBundle(fasciculus bundle);


fasciculus selectFibers(fasciculus bundle, int *selected);
float fiberlen(fiber fib);
void fiberSwap(fiber *f1, fiber *f2, int dim);

#endif // FIBERS_H_INCLUDED