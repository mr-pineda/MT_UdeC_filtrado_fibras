#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <math.h>
#include <stdbool.h>
#include <unistd.h>
#include <limits.h>
#include <time.h>

#include "directory.h"

bool endAnalyser(char* terminator, char* file_name){
    int term_len = strlen(terminator);
    int file_len = strlen(file_name);
    if(term_len > file_len)
        return false;
    return strcmp(terminator,file_name +(file_len-term_len)) == 0;
}

char** extFinder(char* cdir,char* ext,int* n){
    DIR *d;
    struct dirent *dir;
    d = opendir(cdir);
    unsigned int file_count = 0;
    char ** fileNames;
    if (d) {
        while ((dir = readdir(d)) != NULL) {
            if(endAnalyser(ext,(dir->d_name))){
                file_count++;
            }
        }
        if (file_count == 0){
            closedir(d);
            printf("\n  No se encontraron archivos %s\n\n",ext);
            return NULL;
        }
        else{
            printf("\n  Se encontraron %u archivos '%s'\n\n",file_count,ext);
            fileNames = (char**)malloc(sizeof(char*)*file_count);
            rewinddir(d);
            file_count = 0;
            while ((dir = readdir(d)) != NULL) {
                if(endAnalyser(ext,(dir->d_name))){
                    fileNames[file_count] = (char*)malloc(sizeof(char)*FILENAME_MAX);
                    sprintf(fileNames[file_count],"%s",dir->d_name);
                    file_count++;
                }
            }
            closedir(d);
        }
    }
    *n = file_count;
    return fileNames;
}

void replaceChar(char* strDest,int len, char orig, char newChar){

    for(int i=0; i<len; i++){
        if(strDest[i] == orig)
            strDest[i] = newChar;
    }
}

int* name2num(char** nameList,int numFiles){
    int* numList=(int*)malloc(sizeof(int) * numFiles);
    for(int i = 0; i< numFiles;i++){
        numList[i]=atoi(nameList[i]);
    }
    return numList;
}

char** num2name(int* numList, int numFiles){
    char** nameList = (char**)malloc(sizeof(char*)*numFiles);
    for(int i=0;i<numFiles;i++){
        nameList[i] = (char*)malloc(sizeof(char)*FILENAME_MAX);
        sprintf(nameList[i],"%d.bundles",numList[i]);
    }
    return nameList;
}

void cleanNames(char** nameList,int numFiles){
    for(int i =0;i<numFiles;i++){
        free(nameList[i]);
        nameList[i]=NULL;
    }
    free(nameList);
    nameList=NULL;
}

void get_current_dir(char* c_dir){
    if (getcwd(c_dir, PATH_MAX) != NULL) {//Obtenci�n del directorio actual del .exe
        replaceChar(c_dir,strlen(c_dir),'\\','/');//Se cambias los \ por / para evitar conflictos de string
    } else {
        perror("getcwd() error");//Error en caso de que no obtenga el directorio
        exit(1);
    }
}
