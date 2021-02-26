
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <sys/time.h>

struct timeval t0;
struct timeval tic();
float toc(struct timeval begin);
int** one2twodim(int* A, int len);
int*  two2onedim(int** A, int row, int col);
void  testBMM(int **A, int **B, int **D, int size);
int   bin2dec(int *binvec, int len);
int** boolMatrix(int row, int col);
int** matrixOR(int **mat1, int **mat2, int size);
int** matrixAND(int **mat1, int **mat2, int size);
int   vectOR(int v1, int v2);
void  vectCpy(int *v1, int *v2, int start, int end);
int*  RowFromBottom(int **B, int partB, int t, int pos, int n);
void  printArray(int **A, int *B, int len, int dim);
int** calc_rowSums(int **B, int n, int t, int partB);
int* BoolMatrixMult(int *A, int **B, int n, int t, int **Filt, bool filtered);


__device__ void vectOR2(int *v1, int *v2, int n, int *v3, bool newMat){
    
    int id = threadIdx.x + blockDim.x * blockIdx.x;

    if( newMat ){
            v3[id] = ( v1[id] | v2[id] );
    }
    else{
        for(int j = 0; j < n; j++)
            v1[id*n+j] = ( v1[id*n+j] | v2[id*n+j] );
    }
    
}


__global__ void fourRus(int *A, int *RS, int n, int t, int partA, int *Cpar, int *dC){

    int id = threadIdx.x + blockDim.x * blockIdx.x;
    int bin, dec = 0;

    for(int k=0, k2=t-1; k < t; k++, k2--){
        bin =  A[n*id+t*partA+k2];
        dec += bin * pow(2, k);
    }

    for(int j = 0; j < n; j++)
        Cpar[n*id+j] = RS[n*dec+j];
    
    vectOR2(dC, Cpar, n, NULL, false);
}


__global__ void help(int *v1, int *v2, int n, int *v3, bool newMat){
    int id = threadIdx.x + blockDim.x * blockIdx.x;
    v3[id] = ( v1[id] | v2[id] );
}





int main(void){

    int n = 120,
        t = 3;
    
    int **A1 = boolMatrix(n, n);
    int *A = two2onedim(A1, n, n);
    int **B = boolMatrix(n, n);

    int *C = BoolMatrixMult(A, B, n, t, NULL, false);
    int **C2=one2twodim(C, n);
    testBMM(A1, B, C2, n);
    
    free(A);
    free(A1);
    free(B);
    free(C);
    free(C2);
}


int* BoolMatrixMult(int *A, int **B, int n, int t, int **Filt, bool filtered){
  
    int **Cp  = (int **)malloc(n * sizeof(int *));
    int *C   = (int *)malloc(n * n * sizeof(int));
    int *help = (int  *)malloc(n * n * sizeof(int));
    int *RS1 = (int  *)malloc(n * pow(2,t) * sizeof(int));
    for(int g = 0; g < n; g++) {
        Cp[g] = (int *)malloc(n * sizeof(int));
    }

    
    int **RS  = (int **)malloc(pow(2,t) * sizeof(int *));
    for(int g = 0; g < pow(2,t); g++) RS[g]=(int *)malloc(n * sizeof(int));

    int *Cpar, *dA, *dRS, size  = n*n*sizeof(int), *dC,
                          size2 = n*pow(2, t)*sizeof(int);

    cudaMalloc((void **)&dA  , size);
    cudaMalloc((void **)&dRS , size2);
    cudaMalloc((void **)&Cpar, size);
    cudaMalloc((void **)&dC  , size);

    cudaMemcpy(dA, A, size, cudaMemcpyHostToDevice);
    t0 = tic();

    for(int i = 0; i < n/t; i++){

        RS = calc_rowSums(B, n, t, i);   
        RS1 = two2onedim(RS, pow(2, t), n);     
        cudaMemcpy(dRS, RS1, size2, cudaMemcpyHostToDevice);

        fourRus<<<20, n/20>>>(dA, dRS, n, t, i, Cpar, dC);
        cudaDeviceSynchronize();
      
    }
    
    cudaMemcpy(C, dC, size, cudaMemcpyDeviceToHost);
    float dur = toc(t0);
    printf("~ Duration: %f\n", dur);
    
    cudaFree(dA);
    cudaFree(dRS);
    cudaFree(Cpar);
    cudaFree(dC);
    return C;
}



/* All possible row sums of matrix */
int** calc_rowSums(int **B, int n, int t, int partB){

    int indx, bp, k;
    int **RS  = (int **)calloc(pow(2, t) , sizeof(int *));
    int *temp = (int  *)calloc(n , sizeof(int  ));
    for(int g = 0; g < pow(2, t); g++) RS[g] = (int *)calloc(n , sizeof(int));

    bp = 1;
    k  = 0;

    int *dRS, *dtemp, *dres, size = n*sizeof(int);
    cudaMalloc((void **)&dRS, size);
    cudaMalloc((void **)&dtemp, size);
    cudaMalloc((void **)&dres, size);
    
    for(int j = 1; j < pow(2, t); j++){
        
        indx  = j - pow(2, k); 
        temp  = RowFromBottom(B, partB, t, k+1, n);

        cudaMemcpy(dRS, RS[indx], size, cudaMemcpyHostToDevice);
        cudaMemcpy(dtemp, temp, size, cudaMemcpyHostToDevice);

        help<<<20, n/20>>>(dRS, dtemp, n, dres, true);
        cudaMemcpy(RS[j], dres, size, cudaMemcpyDeviceToHost);
        
        if(bp == 1){
            bp = j + 1;
            k ++;
        }
        else bp --;
    }

    return RS;
}


/* Accesses array from bottom */
int* RowFromBottom(int **B, int partB, int t, int pos, int n){
    int *res=(int *)malloc(n * sizeof(int));
    for(int i=0; i<n; i++) res[i]=B[t*partB+t-pos][i];
    return res;
}


/* copies vector */
void vectCpy(int *v1, int *v2, int start, int end){

    for(int i = start; i < end; i++)
        v1[i] = v2[i];
}


/* ORs two vectors */
int vectOR(int v1, int v2){
    
    int v3 = ( v1 | v2 );
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


int** one2twodim(int* A, int len){

    int** D = (int **)calloc(len , sizeof(int *));
    for(int f = 0; f < len; f++) D[f] = (int *)calloc(len , sizeof(int));
    int cnt = 0;

    for(int i = 0; i < len; i++)
        for(int j = 0; j < len; j++)
            D[i][j] = A[cnt++];

    return D;
}


/* 2D to 1D */
int* two2onedim(int** A, int row, int col){

    int* C = (int *)calloc(row * col , sizeof(int));
    int cnt = 0;

    for(int i = 0; i < row; i++)
        for(int j = 0; j < col; j++)
            C[cnt++] = A[i][j];

    return C;
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
    printf("\n");
}


/* BMM tester */
void testBMM(int **A, int **B, int **D, int size){

    int **C = (int **)calloc(size , sizeof(int *));
    for(int i=0; i<size; i++) C[i] = (int *)calloc(size , sizeof(int));

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