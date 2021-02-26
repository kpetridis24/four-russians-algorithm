
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <sys/time.h>

struct timeval t0;
struct timeval tic();
float toc(struct timeval begin);
void  testBMM(int **A, int **B, int **D, int size);
int   bin2dec(int *binvec, int len);
int** boolMatrix(int row, int col);
int** matrixOR(int **mat1, int **mat2, int size);
int** matrixAND(int **mat1, int **mat2, int size);
int*  vectOR(int *v1, int *v2, int len);
void  vectCpy(int *v1, int *v2, int start, int end);
int*  RowFromBottom(int **B, int partB, int t, int pos);
void  printArray(int **A, int *B, int len, int dim);
int** calc_rowSums(int **B, int n, int t, int partB);
int** BoolMatrixMult(int **A, int **B, int n, int t, int **Filt, bool filtered);


void main(){

    float dur;
    int n = 15,
        t = 3;
    
    int **A = boolMatrix(n, n);
    int **B = boolMatrix(n, n);
    int **Filt = boolMatrix(n, n);
    printf("A:\n"); printArray(A, NULL, n, 2); printf("\n");
    printf("B:\n"); printArray(B, NULL, n, 2); printf("\n");

    t0 = tic();
    int **C = BoolMatrixMult(A, B, n, t, Filt, false);
    dur = toc(t0);
    
    printf("C = A x B:\n"); printArray(C, NULL, n, 2); 
    printf("~ Duration: %f\n", dur);
    testBMM(A, B, C, n);
    
    free(A);
    free(B);
    free(C);
}


/* 4-Russian BMM */
int** BoolMatrixMult(int **A, int **B, int n, int t, int **Filt, bool filtered){
  
    int **Cpar  = (int **)malloc(n * sizeof(int *)), **RS, indx;
    int **C     = (int **)malloc(n * sizeof(int *));
    int *binvec = (int  *)malloc(t * sizeof(int));
    for(int g = 0; g < n; g++) {
        Cpar[g] = (int *)malloc(n * sizeof(int));
        C   [g] = (int *)malloc(n * sizeof(int));
    }
    

    for(int i = 0; i < n/t; i++){

        RS = calc_rowSums(B, n, t, i);          /* Get row sums of Bi */

        for(int j = 0; j < n; j++){

            vectCpy(binvec, A[j]+i*t, 0, t);    /* Get vector from Ai */
            indx = bin2dec(binvec, t);          /* Convert binary vector to integer */
            vectCpy(Cpar[j], RS[indx], 0, n);   /* Calc Ci using 4-Russ-Alg */
        }
        matrixOR(C, Cpar, n);
    }
    
    if( filtered ) matrixAND(C, Filt, n);
    return C;
}


/* All possible row sums of matrix */
int** calc_rowSums(int **B, int n, int t, int partB){

    int indx, bp, k;
    int **RS  = (int **)malloc(pow(2, t) * sizeof(int *));
    int *temp = (int  *)malloc(n * sizeof(int  ));
    for(int g = 0; g < pow(2, t); g++) RS[g] = (int *)malloc(n * sizeof(int));

    bp = 1;
    k  = 0;

    for(int j = 1; j < pow(2, t); j++){

        indx  = j - pow(2, k); 
        temp  = RowFromBottom(B, partB, t, k+1);
        RS[j] = vectOR( RS[indx], temp, n );
        
        if(bp == 1){
            bp = j + 1;
            k ++;
        }
        else bp --;
    }

    return RS;
}


/* Accesses array from bottom */
int* RowFromBottom(int **B, int partB, int t, int pos){
    return B[ t*partB+t-pos ];
}


/* copies vector */
void vectCpy(int *v1, int *v2, int start, int end){

    for(int i = start; i < end; i++)
        v1[i] = v2[i];
}


/* ORs two vectors */
int* vectOR(int *v1, int *v2, int len){
    
    int *v3 = (int *)malloc(len * sizeof(int));
    for(int i = 0; i < len; i++){
        v3[i] = ( v1[i] | v2[i] );
    }
    return v3;
}


/* ORs 2-D matrices */
int** matrixOR(int **mat1, int **mat2, int size){

    for(int i = 0; i < size; i++)
        for(int j = 0; j < size; j++)
            mat1[i][j] = ( mat1[i][j] | mat2[i][j] );

    return mat1;
}


/* ANDs 2-D matrices */
int** matrixAND(int **mat1, int **mat2, int size){

    for(int i = 0; i < size; i++)
        for(int j = 0; j < size; j++)
            mat1[i][j] = ( mat1[i][j] & mat2[i][j] );

    return mat1;
}


/* Binary vector to int */
int bin2dec(int *binvec, int len){

    int decval = 0;

    for(int i=0, i2=len-1; i < len; i++, i2--)
        decval += binvec[i2] * pow(2, i);
    
    return decval;
}


/* Random boolean matrix */
int** boolMatrix(int row, int col){

    srand(time( NULL ));
    int** mat = (int **)malloc(row * sizeof(int *));
    for(int i=0; i<row; i++) mat[i] = (int *)malloc(col * sizeof(int));

    for(int i=0; i<row; i++)
        for(int j=0; j<col; j++)
            mat[i][j] = rand() % 2;

    return mat;
}


void printArray(int **A, int *B, int len, int dim){

    if(dim == 2){
        for(int i = 0; i < len; i++){
            for(int j = 0; j < len; j++){
                printf("%d ", A[i][j]);
            }
            printf("\n");
        }
    }
    else{
        for(int i = 0; i < len; i++) printf("%d ", B[i]);
    }
}


/* BMM tester */
void testBMM(int **A, int **B, int **D, int size){

    int **C = (int **)malloc(size * sizeof(int *)), corr = 0;
    for(int i=0; i<size; i++) C[i] = (int *)malloc(size * sizeof(int));

    for(int r=0; r<size; r++)
        for(int i=0; i< size; i++)
            for(int j=0; j<size; j++)
                C[r][i] = C[r][i] | ( A[r][j] & B[j][i] );
            
    for(int i=0; i< size; i++)
        for(int j=0; j<size; j++)
            if( D[i][j] != C[i][j] )
            {
                printf("ERROR!\nIncorrect element calculation!\n");
                exit(0);
            }

    printf("CORRECT!\n");
}


struct timeval tic(){
    
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv;
}


float toc(struct timeval begin){
    
    struct timeval end;
    gettimeofday(&end, NULL);
    float stime = ((double)(end.tv_sec-begin.tv_sec)*1000)+
                            ((double)(end.tv_usec-begin.tv_usec)/1000);
    stime /= 1000;
    return stime;
}