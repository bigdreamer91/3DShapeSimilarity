//
//  main.cpp
//  3DShapeRetrieval
//
//  Created by Geethanjali Jeevanatham on 5/17/16.
//  Copyright © 2016 Geethanjali Jeevanatham. All rights reserved.
//

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "jacobi_eigenvalue.hpp"

/* Usage
 * ./3DShapeRetrieval <filename>.off <filename1>.txt
 * files.txt is used which contains a list of all the filenames in the database
 *
 */


// Function Declaration

int ParseArgs(int argc, char **argv);
double variance(int i, int j, int nverts, float **matrixArray);
int readFile(char *filename, int ***array, float ***midPointX, float ***midPointY, float ***midPointZ,int *nverts,int *it_max);
int computeDemand(int *array,int W, int U, int it_max);
float computeFlowVal(char *str, char *filename);
void tRimLinkedLst();
void appendLinkedList(float val, char *filename);

// Global variables
static char *filename = 0;
static char *filename1 = 0;

// Global declaration of struct node for linked list
struct node{
    float val;
    char *filename;
    struct node* next;
}*head;

int main(int argc, char **argv) {
    
    if (!ParseArgs(argc, argv)) exit(1);
    
    head = NULL;
    struct node *ptr = NULL;
    
    FILE *fp1;
    float flowSum = 0.0;
    int val=3;
    
    // Read files listed in files.txt one by one
    if(!(fp1 = fopen(filename1,"r"))){
        fprintf(stderr, "Unable to open file %s\n",filename1);
        return 0;
    }
    
    int linecount = 0,txtCount=0,txtnCount=0;
    char buffer[1024];
    char str [40];
    while(fgets(buffer, 1023, fp1)){
        linecount++;
        char *bufferp = buffer;
        while(isspace(*bufferp)) bufferp++;
        
        if(*bufferp == '#') continue;
        if(*bufferp == '\0') continue;
        
        if(txtCount==0){
            if((sscanf(bufferp, "%d",&txtCount)!=1)||(txtCount==0)){
                fprintf(stderr, "Syntax error reading header on line %d in file %s\n", linecount,filename1);
                fclose(fp1);
                return NULL;
            }
        }
        else if(txtnCount < txtCount){
            
            if((sscanf(bufferp, "%s", str)!=1)){
                fprintf(stderr, "Syntax error with vertex coordinates on line %d in file %s\n",linecount,filename1);
                fclose(fp1);
                return NULL;
            }
            
            if(strcmp(str, filename)!=0){
                
                // PTD is calculated for each file in the files.txt
                flowSum = computeFlowVal(str,filename);
                //printf("flowSum - %f\n",flowSum);
                
                // each PTD value calculated is appended to a linked list to compare and find the top 5 results
                appendLinkedList(flowSum, str);
            }
            
            txtnCount++;
        }
    }
    
    // Printing results. printing the PTD value and its corresponding filename
    printf("\n");
    ptr = head;
    while(ptr!=NULL){
        printf("val - %f %s    \n",ptr->val,ptr->filename);
        ptr = ptr->next;
    }
    
    printf("\n");
    
    
    return 0;
}

// ParseArgs will parse the command line argument
int ParseArgs(int argc, char **argv)
{
    int print_usage = 0;

    argc--; argv++;
    while (argc > 0) {
        if ((*argv)[0] == '-') {
            if (!strcmp(*argv, "-help")) { print_usage = 1; }
            else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
            argv++; argc--;
        }
        else {
            if (!filename) filename = *argv;
            else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
            argv++; argc--;
            
            if (!filename1) filename1 = *argv;
            else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
            argv++; argc--;
        }
    }
    
    // Check filename
    if (!filename || print_usage) {
        printf("Usage: pca <filename>\n");
        return 0;
    }
    
    // Return OK status
    return 1;
}

// calulates the variance of the data
double variance(int i, int j, int nverts, float **matrixArray){
    int x;
    double mult = 0.0;
    for(x=0;x<nverts;x++){
        mult += (matrixArray[x][i] * matrixArray[x][j]);
    }
    mult = (mult / (nverts - 1));
    return mult;
}

// readfile reads the vertices from input file and creates the PCA applied data matrix, weighted point set array and the midpoints of each 25*25*25 cell
int readFile(char *filename, int ***array, float ***midPointX, float ***midPointY, float ***midPointZ,int *nverts,int *it_max){
    float **matrixArray = nullptr;
    float **transposeData = nullptr;
    float **matrixC = nullptr;
    float **matrixFinal = nullptr;
    float coVar[3][3];
    FILE *fp;
    double sumX = 0.0, sumY = 0.0, sumZ = 0.0;
    int i,ncount=0,x,y,z,nfaces,nedges;
    
    
    
    if(!(fp = fopen(filename,"r"))){
        fprintf(stderr, "Unable to open file %s\n",filename);
        return 0;
    }
    
    int linecount = 0;
    char buffer[1024];
    while(fgets(buffer, 1023, fp)){
        linecount++;
        char *bufferp = buffer;
        while(isspace(*bufferp)) bufferp++;
        
        if(*bufferp == '#') continue;
        if(*bufferp == '\0') continue;
        
        if(*nverts==0){
            if (!strstr(bufferp, "OFF")) {
                if((sscanf(bufferp, "%d%d%d", nverts,&nfaces,&nedges)!=3)||(*nverts==0)){
                    fprintf(stderr, "Syntax error reading header on line %d in file %s\n", linecount,filename);
                    fclose(fp);
                    return NULL;
                }
                
                matrixArray = new float*[*nverts];
                for(int i = 0; i < *nverts; ++i){
                    matrixArray[i] = new float[3];
                    matrixArray[i][0] = 0.0;
                    matrixArray[i][1] = 0.0;
                    matrixArray[i][2] = 0.0;
                }
            }
            
        }
        else if(ncount < *nverts){
            
            // saves the three vertices read in each line in matrixArray against the line number
            if((sscanf(bufferp, "%f%f%f", &matrixArray[ncount][0],&matrixArray[ncount][1],&matrixArray[ncount][2])!=3)){
                fprintf(stderr, "Syntax error with vertex coordinates on line %d in file %s\n",linecount,filename);
                fclose(fp);
                return NULL;
            }
            
            // summation of all the vertices in each column
            // code commented since program works without performing (X_{i} - \bar{X}) computation
           /* sumX += matrixArray[ncount][0];
            sumY += matrixArray[ncount][1];
            sumZ += matrixArray[ncount][2]; */
            ncount++;
        }
    }
    
    // mean of the vertices of each column
    // code commented since program works without performing (X_{i} - \bar{X}) computation
   /* sumX = sumX / *nverts;
    sumY = sumY / *nverts;
    sumZ = sumZ / *nverts; */
    
    // find 3 dimensional co-variance matrix
    coVar[0][0] = variance(0, 0, *nverts,matrixArray);
    coVar[0][1] = coVar[1][0] = variance(0, 1, *nverts,matrixArray);
    coVar[0][2] = coVar[2][0] = variance(0, 2, *nverts,matrixArray);
    coVar[1][1] = variance(1, 1, *nverts,matrixArray);
    coVar[1][2] = coVar[2][1] = variance(1, 2, *nverts,matrixArray);
    coVar[2][2] = variance(2, 2, *nverts,matrixArray);
    
    double a[3*3];
    int it_num;
    int rot_num;
    double v[3*3];
    double d[3];
    
    //calculating eigen vector and eign values calculation - making use of jacobi_eigenvalue.hpp open source code
    jacobi_eigenvalue(3, a, 100, v, d, it_num, rot_num);
    
    double eigVec[3][3];
    
    // eigen vectors calculated are copied into eigVec 2d matrix in decreasing order of eigen values
    eigVec[0][0] = v[6];
    eigVec[0][1] = v[7];
    eigVec[0][2] = v[8];
    eigVec[1][0] = v[3];
    eigVec[1][1] = v[4];
    eigVec[1][2] = v[5];
    eigVec[2][0] = v[0];
    eigVec[2][1] = v[1];
    eigVec[2][2] = v[2];
    
    transposeData = new float*[3];
    matrixC = new float*[3];
    
    for(i=0;i<3;i++){
        transposeData[i] = new float[*nverts];
        matrixC[i] = new float[*nverts];
    }
    
    // transpose of data matrix
    for(i=0;i<3;i++){
        for(int y = 0; y < *nverts;y++){
            transposeData[i][y] = matrixArray[y][i];
        }
    }
    
    // data matrix deleted to free memory
    for(i=0;i<*nverts;i++){
        delete [] matrixArray[i];
    }
    delete [] matrixArray;
    
    // eigVec * trasposeData is calculated and saved in matrixC. This will be the PCA applied final data
    for(int e=0;e<3;e++){
        for(int f=0; f<*nverts;f++){
            matrixC[e][f] = 0;
            for(int g=0; g<3;g++){
                matrixC[e][f] = matrixC[e][f] + (eigVec[e][g] * transposeData[g][f]);
            }
        }
    }
    
    // deleting trasnposeData to free memory
    for(i=0;i<3;i++){
        delete [] transposeData[i];
    }
    delete [] transposeData;
    
    matrixFinal = new float*[*nverts];
    for(i=0;i<*nverts;i++){
        matrixFinal[i] = new float[3];
    }
    
    // transpose of matrixC to get nverts X 3 dimension matrix
    for(i=0;i<*nverts;i++){
        for(int y = 0; y < 3;y++){
            matrixFinal[i][y] = matrixC[y][i];
        }
    }
    
    // deleting matrixC
    for(i=0;i<3;i++){
        delete [] matrixC[i];
    }
    delete [] matrixC;
    
    float leastX = 0.0, leastY = 0.0, leastZ = 0.0;
    float maxX = 0.0, maxY = 0.0, maxZ = 0.0;
    
    // finding the least and max values in each column
    for(i=0;i<*nverts;i++){
        if(matrixFinal[i][0] < leastX){
            leastX = matrixFinal[i][0];
        }
        if(matrixFinal[i][1] < leastY){
            leastY = matrixFinal[i][1];
        }
        if(matrixFinal[i][2] < leastZ){
            leastZ = matrixFinal[i][2];
        }
    }
    
    leastX = 0.0 - leastX;
    leastY = 0.0 - leastY;
    leastZ = 0.0 - leastZ;
    
    // each column values are adjusted so no vertice value is < 0 in the negative axis
    for(i=0;i<*nverts;i++){
        
        if(leastX!=0.0){
            matrixFinal[i][0] += leastX;
        }
        if(leastY!=0.0){
            matrixFinal[i][1] += leastY;
        }
        if(leastZ!=0.0){
            matrixFinal[i][2] += leastZ;
        }
        if(matrixFinal[i][0] > maxX){
            maxX = matrixFinal[i][0];
        }
        if(matrixFinal[i][1] > maxY){
            maxY = matrixFinal[i][1];
        }
        if(matrixFinal[i][2] > maxZ){
            maxZ = matrixFinal[i][2];
        }
    }
    
    // finding the max value in each column and finding the value to split which will act as a bin width.
    if(maxX > 1.0){
        maxX = maxX/25;
    }
    else if(maxX < 1.0){
        maxX = 1.0 / 25;
    }
    
    if(maxY > 1.0){
        maxY = maxY/25;
    }
    else if(maxY < 1.0){
        maxY = 1.0 / 25;
    }
    
    if(maxZ > 1.0){
        maxZ = maxZ/25;
    }
    else if(maxZ < 1.0){
        maxZ = 1.0 / 25;
    }
    
    for(x=0;x<25;x++){
        for(y=0;y<25;y++){
            for(z=0;z<25;z++){
                array[x][y][z] = 0;
            }
        }
    }
    
    // counting the number of vertices that fall in a 25*25*25 cell space. deleting matrixFinal to free memory
    for(i=0;i<*nverts;i++){
        x = matrixFinal[i][0]/maxX;
        y = matrixFinal[i][1]/maxY;
        z = matrixFinal[i][2]/maxZ;
        array[x][y][z] = array[x][y][z] + 1;
        
        // summation of all the vertices in each column for each cell to be used in midpoint calculation
        midPointX[x][y][z] = midPointX[x][y][z] + matrixFinal[i][0];
        midPointY[x][y][z] = midPointY[x][y][z] + matrixFinal[i][1];
        midPointZ[x][y][z] = midPointZ[x][y][z] + matrixFinal[i][2];
        delete [] matrixFinal[i];
    }
    delete [] matrixFinal;
    
    // for each cell in the 25*25*25 cell space finding the midpoint of all the vertices falling in the cell
    for(x=0;x<25;x++){
        for(y=0;y<25;y++){
            for(z=0;z<25;z++){
                if(array[x][y][z]!=0){
                    // midpoint is calculated as the summation of the vertices in each cell divided by the number of vertices in each cell
                    midPointX[x][y][z] = (float)(midPointX[x][y][z] / array[x][y][z]);
                    midPointY[x][y][z] = (float)(midPointY[x][y][z] / array[x][y][z]);
                    midPointZ[x][y][z] = (float)(midPointZ[x][y][z] / array[x][y][z]);
                    *it_max = *it_max+1;
                }
            }
        }
    }
    
    return 1;
}

// for each cell in the 25*25*25 cell space the weighted point set is adjusted to follow the constraint that \sum_{i=1}^{m}(f_{ij}) = (u_{j}W)/U)
int computeDemand(int *array,int W, int U, int it_max){
    int count=0;
            for(int z=0;z<it_max;z++){
                array[z] = ((array[z] * W) / U);
                if(array[z]!=0){
                    count++;
                }
            }
    
    return count;
}

// computes the min flow between supply and demand follwing all the constraints mentioned in PTD
float computeFlowVal(char *str, char *filename){
    int nverts = 0, nverts1 = 0,x,y,z,success,it_max=0,it_max1=0;
    int count = 0;
    int ***supplyArray = new int**[25];
    
    float ***midPointX = new float**[25];
    float ***midPointY = new float**[25];
    float ***midPointZ = new float**[25];
    
    
    int *supplyArrayFinal;
    float *midPointXFinal;
    float *midPointYFinal;
    float *midPointZFinal;
    
    
    float flowSum = 0.0;
    
    // initialization of data structures for the supply data file
    for(x=0;x<25;x++){
        supplyArray[x] = new int*[25];
        midPointX[x] = new float*[25];
        midPointY[x] = new float*[25];
        midPointZ[x] = new float*[25];
        for(y=0;y<25;y++){
            supplyArray[x][y] = new int[25];
            midPointX[x][y] = new float[25];
            midPointY[x][y] = new float[25];
            midPointZ[x][y] = new float[25];
            for(z=0;z<25;z++){
                supplyArray[x][y][z] = 0;
                midPointX[x][y][z] = 0.0;
                midPointY[x][y][z] = 0.0;
                midPointZ[x][y][z] = 0.0;
            }
        }
    }
    
    // save the PCA applied data, weighted point set and midpoints of supply data file
    success = readFile(filename, supplyArray, midPointX, midPointY, midPointZ, &nverts, &it_max);
    
    supplyArrayFinal = new int[it_max];
    midPointXFinal = new float[it_max];
    midPointYFinal = new float[it_max];
    midPointZFinal = new float[it_max];
    
    // processing step to remove zero values and keep only meaningful data
    for(x=0;x<25;x++){
        for(y=0;y<25;y++){
            for(z=0;z<25;z++){
                if(supplyArray[x][y][z]!=0){
                    supplyArrayFinal[count] = supplyArray[x][y][z];
                    midPointXFinal[count] = midPointX[x][y][z];
                    midPointYFinal[count] = midPointY[x][y][z];
                    midPointZFinal[count] = midPointZ[x][y][z];
                    count++;
                }
            }
            delete [] supplyArray[x][y];
            delete [] midPointX[x][y];
            delete [] midPointY[x][y];
            delete [] midPointZ[x][y];
        }
        delete [] supplyArray[x];
        delete [] midPointX[x];
        delete [] midPointY[x];
        delete [] midPointZ[x];
    }
    
    delete [] supplyArray;
    delete [] midPointX;
    delete [] midPointY;
    delete [] midPointZ;
    
    int ***demandArray = new int**[25];
    float ***midPointXD = new float**[25];
    float ***midPointYD = new float**[25];
    float ***midPointZD = new float**[25];
    int *demandArrayFinal = nullptr,*demandArrayFinal1 = nullptr;
    float *midPointXDFinal = nullptr,*midPointXDFinal1 = nullptr;
    float *midPointYDFinal = nullptr,*midPointYDFinal1 = nullptr;
    float *midPointZDFinal = nullptr,*midPointZDFinal1 = nullptr;
    
    // initialization of the data structures for the demand data file
    for(x=0;x<25;x++){
        demandArray[x] = new int*[25];
        midPointXD[x] = new float*[25];
        midPointYD[x] = new float*[25];
        midPointZD[x] = new float*[25];
        for(y=0;y<25;y++){
            demandArray[x][y] = new int[25];
            midPointXD[x][y] = new float[25];
            midPointYD[x][y] = new float[25];
            midPointZD[x][y] = new float[25];
            for(z=0;z<25;z++){
                demandArray[x][y][z] = 0;
                midPointYD[x][y][z] = 0.0;
                midPointXD[x][y][z] = 0.0;
                midPointZD[x][y][z] = 0.0;
            }
        }
    }

    // read the demand data file and save the PCA applied data array, weighted point sets and mid points
    success = readFile(str, demandArray, midPointXD, midPointYD, midPointZD, &nverts1, &it_max1);
    
    demandArrayFinal = new int[it_max1];
    midPointXDFinal = new float[it_max1];
    midPointYDFinal = new float[it_max1];
    midPointZDFinal = new float[it_max1];
    
    count = 0;
    
    // processing step to keep only useful data
    for(x=0;x<25;x++){
        for(y=0;y<25;y++){
            for(z=0;z<25;z++){
                if(demandArray[x][y][z]!=0){
                    demandArrayFinal[count] = demandArray[x][y][z];
                    midPointXDFinal[count] = midPointXD[x][y][z];
                    midPointYDFinal[count] = midPointYD[x][y][z];
                    midPointZDFinal[count] = midPointZD[x][y][z];
                    count++;
                }
            }
            delete [] demandArray[x][y];
            delete [] midPointXD[x][y];
            delete [] midPointYD[x][y];
            delete [] midPointZD[x][y];
        }
        delete [] demandArray[x];
        delete [] midPointXD[x];
        delete [] midPointYD[x];
        delete [] midPointZD[x];
    }
    
    delete [] demandArray;
    delete [] midPointXD;
    delete [] midPointYD;
    delete [] midPointZD;
    
    // changing the demand weighted point set to follow the constraint that total flow between all the supply points and a demand point = (u_{j}W)/U
    success = computeDemand(demandArrayFinal, nverts, nverts1, it_max1);
    
    demandArrayFinal1 = new int[success];
    midPointXDFinal1 = new float[success];
    midPointYDFinal1 = new float[success];
    midPointZDFinal1 = new float[success];
    
    count = 0;
    for(x = 0; x<it_max1;x++){
        if(demandArrayFinal[x]!=0){
            demandArrayFinal1[count] = demandArrayFinal[x];
            midPointXDFinal1[count] = midPointXDFinal[x];
            midPointYDFinal1[count] = midPointYDFinal[x];
            midPointZDFinal1[count] = midPointZDFinal[x];
            count++;
        }
    }
    delete [] demandArrayFinal;
    delete [] midPointXDFinal;
    delete [] midPointYDFinal;
    delete [] midPointZDFinal;
    
    float temp_distance = 0.0;
    int flag = 0, pos;
    float mindist;
    
    // calculates min flow. Finds the min dist between corresponding midpoints and allows flow of weights maintaining the constraints mentioned in PTD
    for(x=0;x<success;x++){
        while(demandArrayFinal1[x]>0){
            flag = 0;
            for(y=0;y<it_max;y++){
                if(supplyArrayFinal[y]>0){
                    temp_distance = sqrt(((midPointXFinal[y] - midPointXDFinal1[x])*(midPointXFinal[y] - midPointXDFinal1[x]))+((midPointYFinal[y] - midPointYDFinal1[x])*(midPointYFinal[y] - midPointYDFinal1[x]))+((midPointZFinal[y] - midPointZDFinal1[x])*(midPointZFinal[y] - midPointZDFinal1[x])));
                    if(flag==0){
                        mindist = temp_distance;
                        flag = 1;
                        pos = y;
                    }
                    else if(temp_distance < mindist){
                        mindist = temp_distance;
                        pos = y;
                    }
                }
            }
            while(supplyArrayFinal[pos]>0 && demandArrayFinal1[x]>0){
                supplyArrayFinal[pos] = supplyArrayFinal[pos]-1;
                demandArrayFinal1[x] = demandArrayFinal1[x]-1;
                // summation of all the flow. Only one flow is allowed per iteration hence only the mindist is summed for each iteration
                flowSum = flowSum + mindist;
            }
        }
    }
    
    // flowSum is dived by nverts as per the PTD objective function
    flowSum = (flowSum / nverts);
    
    return flowSum;
    
}

// maintains the linked list size to 5
void tRimLinkedLst(){
    int len = 0,count=0;
    struct node *current = head;
    struct node *temp = nullptr;
    while(current!=NULL){
        len++;
        current = current->next;
    }
    
    current = head;
    
    if(len > 5){
        while(count<5){
            temp = current;
            current = current->next;
            count++;
        }
        temp->next = NULL;
        
        while(current!=NULL){
            temp = current;
            current = current->next;
            temp = NULL;
            free(temp);
        }

    }

}

// appends each PTD value in a linked list and sorts them
void appendLinkedList(float val, char *filename){
    
    int len = 0, len1 = 0;
    std::string s(filename);
    char *c = new char[s.length()+1];
    strcpy(c,s.c_str());
    struct node *current = head;
    struct node *temp = nullptr;
    temp = (struct node*)malloc(sizeof(struct node));
    temp->val = val;
    temp->filename = c;
    temp->next = NULL;
    
    while(current!=NULL){
        len1++;
        current = current->next;
    }
    
    if(head==NULL){
        head = temp;
    }
    else{
        struct node *current = head;
        struct node *prev = NULL;
        while(current!=NULL){
            if(val < current->val){
                temp->next = current;
                if(prev==NULL){
                    head = temp;
                    break;
                    
                }
                else{
                    prev->next = temp;
                    break;
                }
            }
            prev = current;
            current = current->next;
            len++;
        }
        
        if(len == len1 && len < 6){
            prev->next = temp;
            temp->next = NULL;
        }
    }
    
    // maintains the linked list size to 5
    tRimLinkedLst();
    
}

