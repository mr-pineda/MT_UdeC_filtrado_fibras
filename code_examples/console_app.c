#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <omp.h>

#include "../include/directory.h"
#include "../include/fiberlib.h"
#include "../include/moreMath.h"
#include "../include/sspdlib.h"
#include "../include/isomap.h"

#define DEFAULT_THREADS 1;

int main()
{
    char currentDir[PATH_MAX]; // Directorio de trabajo (donde se ubica el .exe)
    char inputDir[PATH_MAX];   // Directorio donde se ubican archivos .bundle y .bundlesdata
    char outputDir[PATH_MAX];  // Directorio donde se guardarán los archivos creados

    /// CAMBIAR ENTRE MDS / ISOMAP / SSPD_3D SEGUN USO *****
    // Alterna entre imprimir ISOMAP o MDS segun sea necesario
    char *algorithm;
    int algorithmSel;
    const int MAX_THREADS = omp_get_max_threads();

    /// ******************************************

    /// Variables relacionadas con fibras
    fasciculus bundleIn;   // Guarda los datos leidos de los archivos
    fasciculus bundleOut;  // Guarda los datos que se escribirán en archivos
    float **matrix = NULL; // Almacena distancias ssdpd entre fibras
    int *idx_v = NULL;     // vector indices de fibras cuya distancia esta en un percentil especifico

    /// control de tiempo
    clock_t clockSspd, clockIsoMds, clockTotal;
    double timeSspd, timeIsoMds, timeTotal;

    printf("EJEMPLO DE FILTRADO DE CLUSTER EN CONSOLA\n\n");
    printf(" Modo de uso:\n"
           "   Algoritmo disenado para filtrar las fibras segun su distancia SSPD. Debe ubicar\n"
           "   los pares de archivos \".bundles\" y \".bundlesdata\" que desea filtrar en una\n"
           "   carpeta llamada \"InputBundles\", la cual debe estar ubicada en el mismo\n"
           "   directorio del \".exe\" correspondiente a esta consola. El resultado de los\n"
           "   faciculos filtrados se ubicaran en la carpeta \"Filtered_ALGORITHM\", donde \n"
           "   ALGORITHM corresponde a un indicador del algoritmo utilizado.\n"
           "   Este programa usa la libreria OpenMP para ser compatible con procesos multihebra.\n\n");
    printf(" Seleccione un algoritmo:\n"
           "   1 - SSPD 3D: Computo de distancias entre fibras tridimensionales.\n"
           "   2 - MDS + SSPD 2D: Se aplica previamente una reduccion dimensional usando distancia euclidea entre puntos.\n"
           "                      Luego se calcula la distancia entre fibras bidimensionales.\n"
           "   3 - ISOMAP + SSPD 2D: Se aplica previemente una reduccion dimensional usando distancia geodesica entre puntos. \n"
           "                         Luego se calcula la distancia entre fibras bidimensionales.\n\n"
           "   Seleccion: ");

    scanf("%d", &algorithmSel);
    while (algorithmSel < 1 || algorithmSel > 3)
    {
        printf("   Seleccion invalida. Ingrese nuevamente: ");
        scanf("%d", &algorithmSel);
    }

    switch (algorithmSel)
    {
    case 1:
        algorithm = "SSPD3D";
        break;

    case 2:
        algorithm = "MDS";
        break;

    default:
        algorithm = "ISOMAP";
        break;
    }

    get_current_dir(currentDir);
    strcpy(inputDir, currentDir);
    strcat(inputDir, "/InputBundles");

    strcpy(outputDir, currentDir);
    strcat(outputDir, "/Filtered_");
    strcat(outputDir, algorithm);
    mkdir(outputDir);

    if (algorithmSel == 1)
    {
        printf("   Se realizara filtrado de cluster por SSPD 3D\n\n");
    }
    else
    {
        printf("   Se realizara filtrado de cluster por %s + SSPD 2D\n\n", algorithm);
    }
    int numFiles;
    char **names = extFinder(inputDir, ".bundles", &numFiles);

    int threadsUsed;
    printf("   Ingrese el numero de hebras a utilizar (Máxima disponibles %d): ", MAX_THREADS);
    scanf("%d", &threadsUsed);
    if (threadsUsed >= MAX_THREADS)
    {
        printf("   Se utilizara la maxima cantidad de hebras.\n");
        threadsUsed = MAX_THREADS;
    }
    else if (threadsUsed <= 1)
    {
        printf("   Se utilizara 1 hebra.\n");
        threadsUsed = 1;
    }
    else
    {
        printf("   Se utilizaran %d hebras.\n", threadsUsed);
    }
    omp_set_num_threads(threadsUsed);

    float perc;
    printf("Ingrese percentil de fibras a conservar (entre 10%% y 90%%): ");
    scanf("%f", &perc);
    while (perc < 10.0 || perc > 90.0)
    {
        printf("Error ingrese valor de (entre 10%% y 90%%): ");
        scanf("%f", &perc);
    }
    printf("\n\n");

    clockTotal = clock();
    fasciculus bundleBuffer;
    if (algorithmSel == 1)
    {
        for (int i = 0; i < numFiles; i++)
        {
            bundleIn = readBundle(names[i], inputDir);
            int fibersNumber = bundleIn.numFibers;
            printf("bundle: %s  (Contiene %d fibras)\n", names[i], bundleIn.numFibers);

            clockSspd = clock();
            float suma[fibersNumber];
            matrix = ssdpd_matrix_py_op(bundleIn, 3);
            matrix_sum_py(fibersNumber, matrix, suma);
            idx_v = percentile(suma, bundleIn.numFibers, perc);
            bundleOut = selectFibers(bundleIn, idx_v);

            writeBundles(bundleOut, names[i], outputDir);
            timeSspd = (double)(clock() - clockSspd) / CLOCKS_PER_SEC;
            int minSspd = (int)(timeSspd / 60);
            int secSspd = ((int)timeSspd) % 60;
            int msecSspd = ((int)(timeSspd * 1000.0)) % 1000;

            printf("\tFasciculo filtrado:\n"
                   "\t Tiempo: %d min %2d seg %3d ms.\n"
                   "\t Completado %d de %d archivos \n\n",
                   minSspd, secSspd, msecSspd,
                   i + 1, numFiles);

            free(idx_v);
            idx_v = NULL;
            freeMatrix(matrix, bundleIn.numFibers);
            clearBundle(bundleIn);
            clearBundle(bundleOut);
        }
    }
    else
    {
        for (int i = 0; i < numFiles; i++)
        {
            bundleIn = readBundle(names[i], inputDir);
            int fibersNumber = bundleBuffer.numFibers;
            printf("bundle: %s  (Contiene %d fibras)\n", names[i], bundleIn.numFibers);

            clockIsoMds = clock();
            strcpy(bundleBuffer.name, bundleIn.name);
            bundleBuffer.numFibers = bundleIn.numFibers;
            bundleBuffer.fibers = (fiber *)malloc(sizeof(fiber) * bundleBuffer.numFibers);

            for (int j = 0; j < fibersNumber; j++)
            {
                int pointsNumber = bundleBuffer.fibers[j].numPts;
                bundleBuffer.fibers[j].numPts = pointsNumber;
                bundleBuffer.fibers[j].cordinates = (float **)malloc(sizeof(float *) * pointsNumber);

                float temp[2][pointsNumber];
                if (algorithmSel == 2)
                {
                    ISOMAP_2D(pointsNumber, bundleIn.fibers[j].cordinates, 6, 20, temp, dbl_centered_matrix_op);
                }
                else
                {
                    MDS_2D(pointsNumber, bundleIn.fibers[j].cordinates, 20, temp, dbl_centered_matrix_op);
                }

                for (int k = 0; k < pointsNumber; k++)
                {
                    bundleBuffer.fibers[j].cordinates[k] = (float *)malloc(sizeof(float) * 2);
                    bundleBuffer.fibers[j].cordinates[k][0] = temp[0][k];
                    bundleBuffer.fibers[j].cordinates[k][1] = temp[1][k];
                }
            }

            timeIsoMds = (double)(clock() - clockIsoMds) / CLOCKS_PER_SEC;
            int minIsoMds = (int)(timeIsoMds / 60);
            int secIsoMds = ((int)timeIsoMds) % 60;
            int msecIsoMds = ((int)(timeIsoMds * 1000.0)) % 1000;
            printf("\tTiempo %s: %d min %2d seg %3d ms.\n", algorithm, minIsoMds, secIsoMds, msecIsoMds);
            clockSspd = clock();

            float suma[fibersNumber];
            matrix = ssdpd_matrix_py_op(bundleBuffer, 2);
            matrix_sum_py(fibersNumber, matrix, suma);

            idx_v = percentile(suma, bundleBuffer.numFibers, perc);

            bundleOut = selectFibers(bundleIn, idx_v);

            writeBundles(bundleOut, names[i], outputDir);
            timeSspd = (double)(clock() - clockSspd) / CLOCKS_PER_SEC;
            int minSspd = (int)(timeSspd / 60);
            int secSspd = ((int)timeSspd) % 60;
            int msecSspd = ((int)(timeSspd * 1000.0)) % 1000;

            int timeFilter = timeIsoMds + timeSspd;
            int minFilter = (int)(timeFilter / 60);
            int secFilter = ((int)timeFilter) % 60;
            int msecFilter = ((int)(timeFilter * 1000.0)) % 1000;

            printf("\tTiempo SSPD 2D: %d min %2d seg %3d ms.\n", minSspd, secSspd, msecSspd);
            printf("\tFasciculo filtrado:\n"
                   "\t tiempo: %d min %2d seg %3d ms.\n"
                   "\t Completado %d de %d archivos \n\n",
                   minFilter, secFilter, msecFilter,
                   i + 1, numFiles);

            free(idx_v);
            idx_v = NULL;
            freeMatrix(matrix, bundleIn.numFibers);
            clearBundle(bundleIn);
            clearBundle(bundleBuffer);
        }
    }

    timeTotal = (double)(clock() - clockTotal) / CLOCKS_PER_SEC;
    int hrTotal = ((int)timeTotal) / 3600;
    int minTotal = ((int)(timeTotal / 60)) % 60;
    int secTotal = ((int)timeTotal) % 60;
    int msecTotal = ((int)(timeTotal * 1000.0)) % 1000;
    printf("\n\nTiempo total: %d hrs %2d min %2d seg %3d ms.\n\n", hrTotal, minTotal, secTotal, msecTotal);
    system("pause");

    for (int i = 0; i < numFiles; i++)
    {
        free(names[i]);
        names[i] = NULL;
    }
    names = NULL;

    return 0;
}
