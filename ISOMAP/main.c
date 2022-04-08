#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <omp.h>

#include "../headers_sources/moreMath.h"
#include "../headers_sources/directory.h"
#include "../headers_sources/fiberlib.h"
#include "../headers_sources/sspdlib.h"
#include "../headers_sources/isomap.h"

#define MAX_THREADS omp_get_max_threads()

int main(){
    ///Variable globales:
    /// Variables relacionadas con directorio
    char current_dir[PATH_MAX]; //Directorio de trabajo (donde se ubica el .exe)
    char input_dir[PATH_MAX];   //Directorio donde se ubican archivos .bundle y .bundlesdata
    char output_dir[PATH_MAX];  //Directorio donde se guardarán los archivos creados
    char newName[60];           //Prefijo que se le añadirá a los archivos de salida

    /// CAMBIAR ENTRE MDS / ISOMAP SEGUN USO *****
    // Alterna entre imprimir ISOMAP o MDS segun sea necesario
    char* algorithm = "ISOMAP";
    /// ******************************************

    /// Variables relacionadas con fibras
    fasciculus bundleIn;    //Guarda los datos leidos de los archivos
    fasciculus bundleOut;   //Guarda los datos que se escribirán en archivos
    float** matrix = NULL;  //Almacena distancias ssdpd entre fibras
    int* idx_v = NULL;      //vector indices de fibras cuya distancia esta en un percentil especifico

    ///control de tiempo
    clock_t start,start_matrix, start_isomap;
    double clock_time, clock_time_iso;
    double t_prom = 0.0;

    get_current_dir(current_dir);

    strcpy(input_dir,current_dir); //Se crea una copia del directorio
    strcat(input_dir,"/FinalBundles"); //Se especifica el sub-directorio de archivos de entrada (Debe estar ubicado en el mismo directorio del .exe)

    strcpy(output_dir,current_dir);
    strcat(output_dir,"/");
    strcat(output_dir,algorithm);
    strcat(output_dir,"_C");
    mkdir(output_dir);//Se crea un carpeta con el nombre del directorio de salida

    printf("Filtrado de cluster por %s y SSPD-2D\n\n", algorithm);

    int numFiles;
    char**names = extFinder(input_dir,".bundles",&numFiles);




/// Sección para elegir el número de núcleos a utilizar por consola
/// DESCOMENTAR PARA UTILIZAR
/*
    Al inicio, el programa determina el número máximo de hebras disponibles.
    Se le solicita al usuario elegir entre 1 y el máximo de hebras.
*/
    int num_threads = 1; //NÚMERO DE HEBRAS POR SCRIPT
//    printf(" Hebras disponibles del procesador: %d\n",MAX_THREADS);


//    printf(" Ingrese el número de hebras a utilizar [1 - %d]: ",MAX_THREADS);
//    scanf("%d",&num_threads);
//    while(num_threads < 1 || num_threads > MAX_THREADS){
//        printf("Error, ingrese valor entre 1 y %d: ",MAX_THREADS);
//        scanf("%d",&num_threads);
//    }
    omp_set_num_threads(num_threads);


/// Sección percentil de fibras descartadas por consola
/// DESCOMENTAR PARA UTILIZAR
    float perc =70.0;
//    printf("Ingrese percentil de fibras a conservar (entre 10%% y 90%%): ");
//    scanf("%f",&perc);
//    while(perc < 10.0 || perc>90.0){
//        printf("Error ingrese valor de (entre 10%% y 90%%): ");
//        scanf("%f",&perc);
//    }
//    printf("\n\n");

    start=clock();
    fasciculus bundle2D;
    for(int i=0; i< numFiles ; i++){
        bundleIn = readBundle(names[i],input_dir);
        printf("bundle: %s   (Contiene %d fibras)\n",names[i],bundleIn.numFibers);

        start_isomap=clock();
        strcpy(bundle2D.name,bundleIn.name);
        bundle2D.numFibers = bundleIn.numFibers;
        bundle2D.fibers = (fiber*)malloc(sizeof(fiber)*bundle2D.numFibers);

        int nFib=bundle2D.numFibers;
        for(int j=0;j<nFib;j++){

            bundle2D.fibers[j].numPts=bundleIn.fibers[j].numPts;
            bundle2D.fibers[j].cordinates = (float**)malloc(sizeof(float*)*bundle2D.fibers[j].numPts);
            {
                int num_points=bundle2D.fibers[j].numPts;

                float temp[2][num_points];
                /// CAMBIAR AQUI ALGORITMO DE REDUCCION DIMENSIONAL ISOMAP / MDS
                ISOMAP_2D(num_points,bundleIn.fibers[j].cordinates,6,20,temp,dbl_centered_matrix_op);

                for(int k=0 ; k<num_points ; k++){
                    bundle2D.fibers[j].cordinates[k] = (float*)malloc(sizeof(float)*2);
                    bundle2D.fibers[j].cordinates[k][0] = temp[0][k];
                    bundle2D.fibers[j].cordinates[k][1] = temp[1][k];
                }
            }

        }

        clock_time_iso = (double)(clock()-start_isomap)/CLOCKS_PER_SEC;

        int ms_iso = ((int)(clock_time_iso*1000.0))%1000;
        int m_iso = ((int)(clock_time_iso/60))%60;
        int s_iso = ((int)clock_time_iso)%60;
        printf("\tTiempo %s: %d min %d seg %d ms.\n",algorithm,m_iso,s_iso,ms_iso);
        start_matrix = clock();

        float suma[nFib];
        matrix = ssdpd_matrix_py(bundle2D,bundle2D,2);
        matrix_sum_py(nFib, matrix, suma);

        idx_v = percentile(suma,bundle2D.numFibers,perc);

        bundleOut = selectFibers(bundleIn,idx_v);

        strcpy(newName,algorithm);
        strcat(newName,"_");
        strcat(newName,names[i]);

        writeBundles(bundleOut,newName,output_dir);
        clock_time = (double)(clock()-start_matrix)/CLOCKS_PER_SEC;
        clock_time += clock_time_iso;
        t_prom = (1.0/(1.0+(double)i))*(clock_time+(double)i*t_prom);

        int ms_ss = ((int)(clock_time*1000.0))%1000;
        int m_ss = ((int)(clock_time/60))%60;
        int s_ss = ((int)clock_time)%60;

        int ms_p = ((int)(t_prom*1000.0))%1000;
        int m_p = ((int)(t_prom/60))%60;
        int s_p = ((int)t_prom)%60;

        printf("\tFILTRO COMPLETADO:\n"
               "\t tiempo: %d min %d seg %d ms.\n"
               "\t tiempo prom: %d min %d seg %d ms.\n"
               "\t Completado %d de %d archivos \n\n",
               m_ss,s_ss,ms_ss,
               m_p,s_p,ms_p,
               i+1,numFiles);

        free(idx_v);
        idx_v=NULL;
        freeMatrix(matrix,bundleIn.numFibers);
        clearBundle(bundleIn);
        clearBundle(bundle2D);
    }

    clock_time=(double)(clock()-start)/CLOCKS_PER_SEC;
    int horas,minutos, segundos, miliseg;
    miliseg = ((int)(clock_time*1000.0))%1000;
    segundos=((int)clock_time)%60;
    minutos=((int)(clock_time/60))%60;
    horas = ((int)clock_time)/3600;
    printf("\n\nTiempo del algoritmo: %d hrs %d min %d seg %d ms.\n\n",horas,minutos,segundos,miliseg);
    system("pause");

    for(int i=0;i<numFiles;i++){
        free(names[i]);
        names[i]=NULL;
    }
    names=NULL;

    return 0;
}
