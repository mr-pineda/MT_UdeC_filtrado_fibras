#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <math.h>
#include <dirent.h>
#include <stdbool.h>
#include <omp.h>

#include "../headers_sources/directory.h"
#include "../headers_sources/fiberlib.h"
#include "../headers_sources/moreMath.h"
#include "../headers_sources/sspdlib.h"

#define MAX_THREADS omp_get_max_threads()

int main(){

    ///Variable globales:
    /// Variables relacionadas con directorio
    char current_dir[PATH_MAX]; //Directorio de trabajo (donde se ubica el .exe)
    char input_dir[PATH_MAX];   //Directorio donde se ubican archivos .bundle y .bundlesdata
    char output_dir[PATH_MAX];  //Directorio donde se guardarán los archivos creados
    char newName[60];           //Prefijo que se le añadirá a los archivos de salida

    /// Variables relacionadas con fibras
    fasciculus bundleIn;    //Guarda los datos leidos de los archivos
    fasciculus bundleOut;   //Guarda los datos que se escribirán en archivos
    float** matrix = NULL;  //Almacena distancias ssdpd entre fibras
//    float* suma = NULL;     //vector de suma de distancias entre fibras.
    int* idx_v = NULL;      //vector indices de fibras cuya distancia esta en un percentil especifico


    ///control de tiempo
//    clock_t start;
    clock_t start_matrix;
    double t_prom = 0.0;
    double clock_time;

    get_current_dir(current_dir);

    strcpy(input_dir,current_dir); //Se crea una copia del directorio
    strcat(input_dir,"/FinalBundles"); //Se especifica el sub-directorio de archivos de entrada (Debe estar ubicado en el mismo directorio del .exe)

    strcpy(output_dir,current_dir);
    strcat(output_dir,"/SSDPD3D_C"); //Algoritmo análogo para el directorio de salida
    mkdir(output_dir);//Se crea un carpeta con el nombre del directorio de salida

    int numFiles;

    printf("Filtrado de cluster por SSPD 3D\n\n");

    char**names = extFinder(input_dir,".bundles",&numFiles);

    int num_threads = 3;
//    printf(" Ingrese el numero de hebras a utilizar [1 - %d]: ",MAX_THREADS);
//    scanf("%d",&num_threads);
//    while(num_threads < 1 || num_threads > MAX_THREADS){
//        printf("Error, ingrese valor entre 1 y %d: ",MAX_THREADS);
//        scanf("%d",&num_threads);
//    }
    omp_set_num_threads(num_threads);


    float perc = 70.0;
//    printf("Ingrese percentil de fibras a conservar (entre 10%% y 90%%): ");
//    scanf("%f",&perc);
//    while(perc < 10.0 || perc>90.0){
//        printf("Error ingrese valor de (entre 10%% y 90%%): ");
//        scanf("%f",&perc);
//    }
//    printf("\n\n");

//    start=clock();

    for(int i=0; i< numFiles ; i++){
        bundleIn = readBundle(names[i],input_dir);
        int nFib = bundleIn.numFibers;
        printf(" bundle: %s   (Contiene %d fibras)\n",names[i],nFib);
        start_matrix = clock();

        float suma[nFib];
        matrix = ssdpd_matrix_py_op(bundleIn,3);
        matrix_sum_py(nFib, matrix, suma);
        idx_v = percentile(suma,bundleIn.numFibers, perc);

        bundleOut = selectFibers(bundleIn, idx_v);

        strcpy(newName,"SSDPD_3D_");
        strcat(newName,names[i]);

        writeBundles(bundleOut,newName,output_dir);
        clock_time = (double)(clock()-start_matrix)/CLOCKS_PER_SEC;
        t_prom = (1.0/(1.0+(double)i))*(clock_time+(double)i*t_prom);

        int ms = ((int)(clock_time*1000.0))%1000;
        int m = ((int)(clock_time/60))%60;
        int s = ((int)clock_time)%60;

        int ms_p = ((int)(t_prom*1000.0))%1000;
        int m_p = ((int)(t_prom/60))%60;
        int s_p = ((int)t_prom)%60;
        printf("\tCompletado\n"
               "\t tiempo: %d min %d seg %d ms.\n"
               "\t tiempo prom: %d min %d seg %d ms.\n"
               "\t Completado %d de %d archivos \n\n",
               m,s,ms,
               m_p,s_p,ms_p,
               i+1,numFiles);

        freeMatrix(matrix,bundleIn.numFibers);
//        free(suma);
//        suma=NULL;
        free(idx_v);
        idx_v=NULL;
        clearBundle(bundleIn);
        clearBundle(bundleOut);
    }

//    clock_time=(double)(clock()-start)/CLOCKS_PER_SEC;
//    int horas,minutos, segundos, miliseg;
//    miliseg = ((int)(clock_time*1000.0))%1000;
//    segundos=((int)clock_time)%60;
//    minutos=((int)(clock_time/60))%60;
//    horas = ((int)clock_time)/3600;
//    printf("\n\nTiempo del algoritmo: %d hrs %d min %d seg %d ms.\n\n",horas,minutos,segundos,miliseg);
//    system("pause");

    for(int i=0;i<numFiles;i++){
        free(names[i]);
        names[i]=NULL;
    }
    names=NULL;

    return 0;
}
