#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <math.h>
#include <stdbool.h>
#include <unistd.h>
#include <limits.h>
#include <time.h>

#include "moreMath.h"

float eucDisPP(float* p1, float*p2){
    return sqrt((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1])+(p1[2]-p2[2])*(p1[2]-p2[2]));
};

float eucDisPP_2d(float* p1, float*p2){
    return sqrt((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1]));
};

float* matrix_sum(float** matrix,int rows, int cols, short axis){
    float* sum;
    if(axis){
        sum = (float*)malloc(sizeof(float)*rows);
        for (int i=0; i<rows; i++ ){
            sum[i]=0;
            for(int j=0;j<cols;j++){
                sum[i] += matrix[i][j];
            }
        }
        return sum;
    }
    else{
        sum = (float*)malloc(sizeof(float)*cols);
        for (int i=0; i<cols; i++ ){
            sum[i]=0;
            for(int j=0;j<rows;j++){
                sum[i] += matrix[j][i];
            }
        }
        return sum;
    }
};

float* matrix_sum_op(float** matrix,int len){
    float* sum;
    sum = (float*)malloc(sizeof(float)*len);
    for (int i=0; i<len; i++ ){
        sum[i]=0;
        for(int j=0;j<len;j++){
            sum[i] += matrix[i][j];
        }
    }
    return sum;
};

void matrix_sum_py(int len, float** matrix,float sum[len]){
    for (int i=0; i<len; i++ ){
        sum[i]=0;
        for(int j=0;j<len;j++){
            sum[i] += matrix[i][j];
        }
    }
};

void m_transpose(int row, int col, float matrix[row][col], float mat_tr[col][row]){
    for(int i=0 ; i<row ; i++){
        for(int j=0 ; j<col ; j++){
            mat_tr[j][i] = matrix[i][j];
        }
    }
}

void quicksort(float *arr, int N){
    int i, j=0, pivot;
    float temp;
    if (N < 2)
        return;
    pivot = arr[0];
    for (i = 1; i < N; i++) {
        if (arr[i] < pivot) {
            for (j = i - 1; j >= 0; j--) {
                if (arr[j] <= pivot)
                    break;
                temp = arr[j];
                arr[j] = arr[j + 1];
                arr[j + 1] = temp;
            }
        }
    }
    quicksort(arr, j);
    quicksort(arr + j + 1, N - j - 1);
}

int* percentile(float arr[], int arr_size, float perc){
    int i, j, count,num_ele=0;
    int* idx_vector;
    float percent;
    for (i = 0; i < arr_size; i++) {
        count = 0;
        for (j = 0; j < arr_size; j++) {
            if (arr[i] > arr[j]) {
                count++;
            }
        }
        percent = (count * 100) / (arr_size - 1);
        if(percent <= perc){
            num_ele++;
        }
    }

    idx_vector = (int*)malloc(sizeof(int)*(num_ele+1));
    idx_vector[0]=num_ele;
    num_ele  = 1;
    for (i = 0; i < arr_size; i++) {
        count = 0;
        for (j = 0; j < arr_size; j++) {
            if (arr[i] > arr[j]) {
                count++;
            }
        }
        percent = (count * 100) / (arr_size - 1);
        if(percent <= perc){
            idx_vector[num_ele++] = i;
        }
    }
    return idx_vector;
}


void freeMatrix(float** matrix,int rows){
    for(int i=0; i<rows; i++){
        free(matrix[i]);
    }
    free(matrix);
}
