#ifndef DIRECTORY_H_INCLUDED
#define DIRECTORY_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <math.h>
#include <stdbool.h>
#include <unistd.h>
#include <limits.h>
#include <time.h>

bool endAnalyser(char* terminator, char* file_name); // Verifica si un archivo tiene una extensión específica
char** extFinder(char* cdir,char* ext,int* n); // Crea una lista de strings de archivos con una extensión específica
void replaceChar(char* strDest,int len, char orig, char newChar); //Reemplaza los caracteres de un string por otro especificado
int* name2num(char** nameList,int numFiles);
char** num2name(int* numList, int numFiles);
void cleanNames(char** nameList,int numFiles);
void get_current_dir(char c_dir[]);

#endif // DIRECTORY_H_INCLUDED
