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

#include "fiberlib.h"
#include "moreMath.h"

fasciculus readBundle(char *fileName, char *path)
{
    fasciculus bundle;
    bundle.fibers = NULL;

    char namePath[PATH_MAX];
    strcpy(namePath, path);
    strcat(namePath, "/");
    strcat(namePath, fileName);

    FILE *bundlesFile;
    bundlesFile = fopen(namePath, "r");

    if (bundlesFile == NULL)
    {
        printf("Error opening the file %s\n", fileName);
        exit(1);
    }

    char tempbuff[USHRT_MAX];
    char bundleName[60];
    for (int i = 0; i < 2; i++)
        fgets(tempbuff, USHRT_MAX, bundlesFile);
    fgets(tempbuff, 18, bundlesFile);
    fscanf(bundlesFile, "%s", bundleName);
    strcpy(bundle.name, strtok(bundleName, "'"));
    for (int i = 0; i < 2; i++)
        fgets(tempbuff, USHRT_MAX, bundlesFile);
    fgets(tempbuff, 21, bundlesFile);
    fscanf(bundlesFile, "%u", &bundle.numFibers);
    fclose(bundlesFile);

    FILE *bundlesdataFile;
    strcat(namePath, "data");
    bundlesdataFile = fopen(namePath, "rb");
    if (bundlesdataFile == NULL)
    {
        printf("Error opening the file %s\nFull path: %s", fileName, namePath);
        exit(2);
    }

    bundle.fibers = (fiber *)malloc(bundle.numFibers * sizeof(fiber));
    for (int i = 0; i < bundle.numFibers; i++)
    {
        fread(&bundle.fibers[i].numPts, sizeof(int), 1, bundlesdataFile);
        bundle.fibers[i].cordinates = (float **)malloc(bundle.fibers[i].numPts * sizeof(float *));
        for (int j = 0; j < bundle.fibers[i].numPts; j++)
        {
            bundle.fibers[i].cordinates[j] = (float *)malloc(3 * sizeof(float));

            for (int k = 0; k < 3; k++)
            {
                fread(&bundle.fibers[i].cordinates[j][k], sizeof(float), 1, bundlesdataFile);
            }
        }
    }
    fclose(bundlesdataFile);
    return bundle;
}

void writeBundles(fasciculus bundle, char *fileName, char *path)
{
    FILE *bundlesFile;
    FILE *bundlesdataFile;
    char name[PATH_MAX];
    strcpy(name, path);
    strcat(name, "/");
    strcat(name, fileName);
    bundlesFile = fopen(name, "w");
    char fileBundle[100];
    sprintf(fileBundle, "['%s',0]", bundle.name);

    if (bundlesFile == NULL)
    {
        printf("Error opening the file %s\n", fileName);
        exit(1);
    }
    fprintf(bundlesFile, "attributes = {\n    'binary' : 1,\n    'bundles' : %s,\n    'byte_order' : 'DCBA',\n    'curves_count' : %d,\n    'data_file_name' : '*.bundlesdata',\n    'format' : 'bundles_1.0',\n    'space_dimension' : 3\n  }", fileBundle, bundle.numFibers);
    fclose(bundlesFile);

    strcat(name, "data");
    bundlesdataFile = fopen(name, "wb");
    if (bundlesdataFile == NULL)
    {
        printf("Error opening the file %s\n", name);
        exit(2);
    }

    for (int i = 0; i < bundle.numFibers; i++)
    {
        fwrite(&bundle.fibers[i].numPts, sizeof(int), 1, bundlesdataFile);
        for (int j = 0; j < bundle.fibers[i].numPts; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                fwrite(&bundle.fibers[i].cordinates[j][k], sizeof(float), 1, bundlesdataFile);
            }
        }
    }

    fclose(bundlesdataFile);
};

void clearBundle(fasciculus bundle)
{
    for (int i = 0; i < bundle.numFibers; i++)
    {
        for (int j = 0; j < bundle.fibers[i].numPts; j++)
        {
            free(bundle.fibers[i].cordinates[j]);
            bundle.fibers[i].cordinates[j] = NULL;
        }
        free(bundle.fibers[i].cordinates);
        bundle.fibers[i].cordinates = NULL;
    }
    free(bundle.fibers);
    bundle.fibers = NULL;
}

fasciculus selectFibers(fasciculus bundleIn, int *selected)
{
    fasciculus bundleOut;
    bundleOut.numFibers = selected[0];
    strcpy(bundleOut.name, bundleIn.name);
    bundleOut.fibers = (fiber *)malloc(sizeof(fiber) * bundleOut.numFibers);
    for (int i = 1; i <= bundleOut.numFibers; i++)
    {
        bundleOut.fibers[i - 1].numPts = bundleIn.fibers[selected[i]].numPts;
        bundleOut.fibers[i - 1].cordinates = (float **)malloc(sizeof(float *) * bundleOut.fibers[i - 1].numPts);
        for (int j = 0; j < bundleOut.fibers[i - 1].numPts; j++)
        {
            bundleOut.fibers[i - 1].cordinates[j] = (float *)malloc(sizeof(float) * 3);
            for (int k = 0; k < 3; k++)
            {
                bundleOut.fibers[i - 1].cordinates[j][k] = bundleIn.fibers[selected[i]].cordinates[j][k];
            }
        }
    }
    return bundleOut;
}

float fiberlen(fiber fib)
{
    float len = 0;
    for (int i = 1; i < fib.numPts; i++)
    {
        len += eucDisPP(fib.cordinates[i - 1], fib.cordinates[i]);
    }
    return len;
}

void fiberSwap(fiber *f1, fiber *f2, int dim)
{
    int tempNum = f1->numPts;
    float **tempCord = (float **)malloc(sizeof(float *) * f1->numPts);
    for (int i = 0; i < f1->numPts; i++)
    {
        tempCord[i] = (float *)malloc(sizeof(float) * dim);
        for (int j = 0; j < 3; j++)
        {
            tempCord[i][j] = f1->cordinates[i][j];
        }
        free(f1->cordinates[i]);
    }
    free(f1->cordinates);

    f1->numPts = f2->numPts;
    f1->cordinates = (float **)malloc(sizeof(float *) * f2->numPts);
    for (int i = 0; i < f2->numPts; i++)
    {
        f1->cordinates[i] = (float *)malloc(sizeof(float) * dim);
        for (int j = 0; j < 3; j++)
        {
            f1->cordinates[i][j] = f2->cordinates[i][j];
        }
        free(f2->cordinates[i]);
    }
    free(f2->cordinates);

    f2->numPts = tempNum;
    f2->cordinates = (float **)malloc(sizeof(float *) * tempNum);
    for (int i = 0; i < tempNum; i++)
    {
        f2->cordinates[i] = (float *)malloc(sizeof(float) * dim);
        for (int j = 0; j < 3; j++)
        {
            f2->cordinates[i][j] = tempCord[i][j];
        }
        free(tempCord[i]);
    }
    free(tempCord);
}
