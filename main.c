#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#define MAX_POINTS 1000
#define MAX_DIMS  10
#define LINE_LENGTH 1024

//bool type (false and true)
typedef enum bool{
    false,
    true
} bool;

// Data structure to deal with given points from file
typedef struct Data {
    double  array [MAX_POINTS][MAX_DIMS];
    int n;
    int m;

}Data;

//Prints all points in data matrix in format (separated by commas and each point in new line)
void print_data(Data *mat) {
    int n = mat->n;
    int m = mat->m;
    int i, j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            printf("%f", mat->array[i][j]);
            if (j != m - 1) {
                printf(",");
            }
        }
        printf("\n");
    }
}

// Given a path to file, return a Data structure that contains all points
Data load_data(char *path) {

    //Initialize some variables needed later
    char ch;
    char *token;
    double value;
    int i = 0;
    int j = 0;
    Data data;
    char line[LINE_LENGTH];

    //Open file and check if legit
    FILE *fptr;
    fptr = fopen(path, "r");

    //Deal with case that file is not legit
    assert(fptr != NULL);

    // Go line by line (Assuming LINE_Length is enough) and insert it to Data
    while (fgets(line, sizeof(line), fptr)) {
        j = 0;
        token = strtok(line, ",");
        while (token != NULL) {
            //First convert it to double and insert it
            value = strtod(token, NULL);
            data.array[i][j] = value;
            //Advance to next iteration
            token = strtok(NULL, ",");
            j++;
        }
        i++;
    }

    //Finish creating Data matrix and update relevant parameters
    data.n = i;
    data.m = j;
    fclose(fptr);
    return data;
}


//*********************************************************************************//
//************************   Matrix Related Functions  ****************************//
//*********************************************************************************//


// Data structure to deal with mathematical matrices
typedef struct Matrix {
    double **array;
    int n;
    int m;
} Matrix;

//Given integers n and m, return a zeros matrix shape (n,m)
Matrix zeros(int n, int m) {

    //Initialize variables needed later
    Matrix mat;
    mat.array= calloc(n,sizeof(double ));
    double *row;
    int i, j;

    //Deal with errors in allocating memory
    assert (mat.array!=NULL);

    //update parameters of Matrix M
    mat.n = n;
    mat.m = m;


    // Create array needed
    for (i = 0; i < n; i++) {
        row = calloc(m, sizeof(double));

        //Deal with errors in allocating memory
        assert (row!=NULL);

        mat.array[i] = row;
    }

    return mat;
}

//Frees Matrix (inner Free)
void free_mat(Matrix *mat) {
    int n = mat->n;
    int i;

    //free each row
    for (i = 0; i < n; i++) {
        free(mat->array[i]);
    }
    //free columns array
    free(mat->array);
}

//Prints Matrix in format (separated by commas and each point in a new line)
void print_mat(Matrix* mat) {
    int i, j;
    int n= mat->n;
    int m= mat->m;

    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            printf("%f", mat->array[i][j]);
            if (j != m - 1) {
                printf(",");
            }
        }
        printf("\n");
    }
}

//Given a Matrix A, returns a Matrix A.T s.t for all i and j A[i,j] = A.T[j,i]
Matrix transpose(Matrix* A){

    //Initialize Paramaters
    int i,j;
    int n = A->n;
    int m = A->m;

    Matrix A_t = zeros(m,n);
    for (i=0; i<n;i++){
        for (j=0; j<n; j++){
            A_t.array[i][j] = A->array[j][i];
        }
    }
    return A_t;
}

//returns a new matrix B which is exact copy of A. for all i,j A[i,j] = B[i,j]
Matrix mat_copy(Matrix* A){
    int i=0,j=0;
    Matrix B = zeros(A->n,A->m);
    for (i = 0; i < A->n; i++) {
        for (j = 0; j < A->m; j++) {
            B.array [i][j] = A->array[i][j];
        }
    }
    return B;
}

//Given a Matrix A, performs power operation of all diagonal
Matrix mat_pow_diagonal(Matrix* A, double pow){
    //Initialize Variables
    int n = A->n;
    int m = A->m;
    int i;

    //Assert Matrix is squared
    assert(n==m);

    //Create a copy and then power each element on diagonal
    Matrix B = mat_copy(A);
    for (i=0; i<n; i++){
        B.array[i][i] = powf(A->array[i][i],pow);
    }
    return B;
}

//returns true if and only if for all i,j A[i,j]=B[i,j]
bool mat_is_equal(Matrix* A,Matrix* B){
    int i,j;
    int n1 = A->n;
    int n2 = B->n;
    int m1 = A->m;
    int m2 = B->m;
    if ((n1!=n2) || (m1!=m2)){
        return false;
    }
    for (i=0; i<n1;i++){
        for (j=0; j<m1;j++){
            if (A->array[i][j]!= B->array[i][j]){
                return false;
            }
        }
    }
    return true;
}

//Performs Matrix multiplication
Matrix mat_mul (Matrix* A,Matrix* B){
    //Initialize parameters need later
    int i,j,k,sum;
    int n1 = A->n;
    int n2 = B->n;
    int m1 = A->m;
    int m2 = B->m;

    //Check if it is ok to multiply
    assert (m1==n2);

    //Initialize Matrix C=A*B
    Matrix C = zeros(n1,m2);

    //Perform Matrix Multiplication
    for (i=0; i<n1;i++){
        for (j=0;j<m2;j++){
            sum =0;
            for (k=0;k<m1;k++){
                sum += A->array[i][k] * B->array[k][j];
            }
            C.array[i][j] = sum;
        }
    }
    return C;
}

//Multiply a matrix by a scalar (entry by entry)
Matrix mat_scalar_mul(Matrix* A, double scalar){
    //Initialize Parameters needed later
    int i,j;
    int n = A->n;
    int m = A->m;

    Matrix C = zeros(n,m);

    //do the multiplicity entry by entry
    for (i=0;i<n;i++){
        for(j=0;j<m;j++){
            C.array[i][j] = A->array[i][j]*scalar;
        }
    }
    return C;

}

//Returns the identity Matrix I[i,j]=0 if i!=j and I[i,i] = 1
Matrix mat_identity (int n){
    int i=0;
    Matrix I = zeros(n,n);
    for (i=0;i<n;i++){
        I.array[i][i]=1;
    }
    return I;
}

// Returns a matrix C s.t C[i,j]=A[i,j]+B[i,j] for all i,j. Meaning C = A+B
Matrix mat_add (Matrix* A, Matrix*B){
    // Initialize some variables need later
    int i,j;
    int n1 = A->n;
    int n2 = B->n;
    int m1 = A->m;
    int m2 = B->m;

    //Make sure operation is legit
    assert(n1==n2 && m1==m2);

    Matrix C =zeros(n1,m1);
    for (i=0; i<n1; i++){
        for (j=0; j<m1; j++){
            C.array[i][j] = A->array[i][j] +B->array[i][j];
        }
    }
    return C;
}

// Returns a matrix C s.t C[i,j]=A[i,j]-B[i,j] for all i,j. Meaning C = A-B
Matrix mat_sub (Matrix* A, Matrix* B){
    // Initialize some variables need later
    int i,j;
    int n1 = A->n;
    int n2 = B->n;
    int m1 = A->m;
    int m2 = B->m;

    //Make sure operation is legit
    assert(n1==n2 && m1==m2);

    Matrix C =zeros(n1,m1);
    for (i=0; i<n1; i++){
        for (j=0; j<m1; j++){
            C.array[i][j] = A->array[i][j] - B->array[i][j];
        }
    }
    return C;
}

//returns a Matrix B s.t B[i,j] = A[i,j]^2 for all i,j.
Matrix mat_square(Matrix* A){
    Matrix B = zeros(A->n,A->m);
    int i,j;
    for (i=0;i<A->n;i++){
        for (j=0;j<A->m;j++){
            B.array [i][j] = powf(A->array[i][j],2);
        }
    }
    return B;
}

//returns a double which is the sum of all entries in the row_index row in the matrix
double mat_sum_by_row(Matrix* A, int row_index){
    double sum =0;
    int j;
    for (j=0; j<A->m;j++){
        sum += A->array[row_index][j];
    }
    return sum;
}

//returns a double which is the sum of all entries in the row_index row in the matrix
double mat_total_sum (Matrix* A){
    double sum =0;
    int i;
    for (i=0; i<A->n;i++){
        sum += mat_sum_by_row(A,i);
    }
    return sum;
}


//**************** Tests *******************
void test_transpose_and_copy(){
    printf("--------------------------------------------------\n");
    printf("---Tests for Transpose and copy Functions----\n");
    printf("--------------------------------------------------\n");
    Matrix A = zeros(3,3);
    A.array[0][0] = 1;
    A.array[0][1] = 2;
    A.array[0][2] = 3;
    A.array[1][0] = 4;
    A.array[1][1] = 5;
    A.array[1][2] = 6;
    A.array[2][0] = 7;
    A.array[2][1] = 8;
    A.array[2][2] = 9;


    printf("This is A:\n");
    print_mat(&A);
    printf("\n");
    printf("This is Trasnpose(A):\n");

    Matrix A_t = transpose(&A);
    print_mat(&A_t);
    printf("\n");

    Matrix B = mat_copy(&A);
    printf("This is B= mat_copy(A):\n");
    print_mat(&B);
    printf("\n");

    Matrix B_t = mat_copy(&A_t);
    printf("This is B.T= mat_copy(A.T):\n");
    print_mat(&B);
    printf("\n");

    free_mat(&A);
    free_mat(&A_t);
    free_mat(&B);
    free_mat(&B_t);

}
void test_mat_mul(){
    printf("--------------------------------------------------\n");
    printf("---Tests for Matrix Multiplication----\n");
    printf("--------------------------------------------------\n");
    Matrix A = zeros(2,4);
    A.array[0][0] = 3;
    A.array[0][1] = 2;
    A.array[0][2] = 1;
    A.array[0][3] = 5;

    A.array[1][0] = 9;
    A.array[1][1] = 1;
    A.array[1][2] = 3;
    A.array[1][3] = 0;


    Matrix B = zeros(4,3);
    B.array[0][0] = 2;
    B.array[0][1] = 9;
    B.array[0][2] = 0;

    B.array[1][0] = 1;
    B.array[1][1] = 3;
    B.array[1][2] = 5;

    B.array[2][0] = 2;
    B.array[2][1] = 4;
    B.array[2][2] = 7;

    B.array[3][0] = 8;
    B.array[3][1] = 1;
    B.array[3][2] = 5;

    printf("This is C= A*B (Matrix multiplication):\n");
    Matrix C = mat_mul(&A,&B);
    print_mat(&C);
    printf("\n");
}
void test_mat_scalar_mul(){
    printf("--------------------------------------------------\n");
    printf("---Tests for Scalar and matrix multiplication----\n");
    printf("--------------------------------------------------\n");
    Matrix B = zeros(4,3);
    B.array[0][0] = 2;
    B.array[0][1] = 9;
    B.array[0][2] = 0;

    B.array[1][0] = 1;
    B.array[1][1] = 3;
    B.array[1][2] = 5;

    B.array[2][0] = 2;
    B.array[2][1] = 4;
    B.array[2][2] = 7;

    B.array[3][0] = 8;
    B.array[3][1] = 1;
    B.array[3][2] = 5;

    printf("This is the Matrix B:\n");
    print_mat(&B);
    printf("\n");

    printf("This is C= 4*B (Matrix multiplication):\n");
    Matrix A = mat_scalar_mul(&B,4);
    print_mat(&A);
    printf("\n");
}
void test_power_diag(){
    printf("--------------------------------------------------\n");
    printf("---Tests for Power Diagonal----\n");
    printf("--------------------------------------------------\n");
    Matrix A = zeros(3,3);
    A.array[0][0] = 1;
    A.array[0][1] = 2;
    A.array[0][2] = 3;
    A.array[1][0] = 4;
    A.array[1][1] = 5;
    A.array[1][2] = 6;
    A.array[2][0] = 7;
    A.array[2][1] = 8;
    A.array[2][2] = 9;


    printf("This is the Matrix A:\n");
    print_mat(&A);
    printf("\n");


    printf("This is the Matrix B = sqrt(A) :\n");
    Matrix B= mat_pow_diagonal(&A,0.5);
    print_mat(&B);
    printf("\n");


}
void test_mat_iden(){
    printf("--------------------------------------------------\n");
    printf("-----Tests for Identity----\n");
    printf("--------------------------------------------------\n");
    int i;
    Matrix I;
    for (i=1; i<=5;i++){
        I = mat_identity(i);
        printf("This is the identity(%d) Matrix:\n",i);
        print_mat(&I);
        printf("\n");
        //free_mat(&I);
    }

}
void test_sums(){
    Matrix A = zeros(3,3);
    A.array[0][0] = 1;
    A.array[0][1] = 2;
    A.array[0][2] = 3;
    A.array[1][0] = 4;
    A.array[1][1] = 5;
    A.array[1][2] = 6;
    A.array[2][0] = 7;
    A.array[2][1] = 8;
    A.array[2][2] = 9;

    Matrix squared = mat_square(&A);
    Matrix B = mat_copy(&A);
    Matrix C = mat_sub(&squared,&B);
    double result =  mat_total_sum(&C);
    printf("The sum of all elements is (should be 240): %f \n",result);
    Matrix D = mat_add(&C,&A);
    result = mat_total_sum(&D);
    printf("The sum of all elements is (should be 285): %f \n",result);
    free_mat(&A);
    free_mat(&B);
    free_mat(&C);
    free_mat(&D);
    free_mat(&squared);
}
void all_mat_tests(){
    test_transpose_and_copy();
    test_mat_mul();
    test_mat_scalar_mul();
    test_power_diag();
    test_mat_iden();
    test_sums();
}
//********************************************


int main(){
    Data data = load_data("C:\\Users\\Tomer\\CLionProjects\\SpectralKMeans\\tests\\input.txt");
    all_mat_tests();


}
