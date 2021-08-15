#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#define MAX_POINTS 1000
#define MAX_DIMS  10
#define LINE_LENGTH 1024
#define EPSILON 0.001

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
    int i,j,k;
    double sum;
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

//given index i returns the i'th row in the Matrix A
double* mat_get_row (Matrix* A, int i){
    int n = A->n;
    int m = A->m;
    assert(i<n && i>=0);
    int j;
    double* row = calloc(m ,sizeof(double ));
    for (j=0; j<m; j++){
        row[j] = A->array[i][j];
    }
    return row;
}

//given index j returns the j'th row in the Matrix A
double* mat_get_col (Matrix* A, int j){
    int n = A->n;
    int m = A->m;
    assert(j<m && j>=0);
    int i;
    double* col = calloc(m , sizeof(double ));
    for (i=0; i<m; i++){
        col[i] = A->array[i][j];
    }
    return col;
}

//given index i returns the i'th row in the Data type data
double* data_get_row (Data* data, int i){
    int n = data->n;
    int m = data->m;
    assert(i<n && i>=0);
    int j;
    double* row = calloc(m ,sizeof(double ));
    for (j=0; j<m; j++){
        row[j] = data->array[i][j];
    }
    return row;
}

//given index j returns the j'th row in the Data type data
double* data_get_col (Data* data, int j){
    int n = data->n;
    int m = data->m;
    assert(j<m && j>=0);
    int i;
    double* col = calloc(m , sizeof(double ));
    for (i=0; i<m; i++){
        col[i] = data->array[i][j];
    }
    return col;
}

//*********************************************************************************//
//************************   Algorithm Related Functions  ****************************//
//*********************************************************************************//

//returns the l2 norm of 2 points (point1, point2)
double l2_norm (double* point1, double* point2, int n){
    int i;
    double sum =0;
    double value;

    for (i=0; i<n ; i++){
        value = point1[i] - point2[i];
        value = pow(value,2);
        sum += value;
    }
    return sqrt(sum);
}

//returns 1 if num>=0 else returns -1
int sign (double num){
    if (num<0){
        return -1;
    }
    else{
        return 1;
    }
}

//Build the weighted adjacency matrix W[i,j] = exp(-l2(p1,p2)/2)
Matrix build_W (Data* data){
    //Initialize Variables
    int i,j;
    int n = data->n;
    int d = data->m;
    Matrix W = zeros(n,n);
    double* point1;
    double* point2;
    double l2;
    double value;

    //Loop and calculate W
    for (i=0 ; i<n ; i++){
        for (j=i+1; j<n; j++){
            point1 = data_get_row(data, i);
            point2 = data_get_row(data, j);
             l2 = l2_norm(point1, point2, d);
            value = exp(-l2/2);
            W.array[i][j] = value;
            W.array[j][i] = value;
            free(point1);
            free(point2);
        }
    }
    return W;
}

//Build Diagonal Degree Matrix D, given the weighted Matrix (must be squared) W
Matrix build_D (Matrix* W){

    int n = W->n;
    int i;
    Matrix D = zeros(n,n);
    for (i=0; i<n ; i++){
        D.array[i][i] = mat_sum_by_row(W,i);
    }
    return D;
}

//Build D_half given the matrix D , given the Diagonal degree matrix D
Matrix build_D_half (Matrix* D){
    return mat_pow_diagonal(D, -0.5);
}

//Build l_norm Matrix (The laplacian) of the graph
Matrix laplacian (Matrix* D_half, Matrix* W){
    int n = D_half->n;

    //Calculate everything necessary
    Matrix temp_val1 = mat_mul(D_half,W); //D_half @ W
    Matrix temp_val2 = mat_mul(&temp_val1,D_half); // D_half @ W @ D_half
    Matrix I = mat_identity(n);
    Matrix ret_mat = mat_sub(&I,&temp_val2); //I - D_half @ W @ D_half

    //free everything that is needed
    free_mat(&temp_val1);
    free_mat(&temp_val2);
    free_mat(&I);

    return ret_mat;
}


//Struct to Store entry of the matrix (i,j)
typedef struct Index{
    int x;
    int y;
} Index;

//Returns a struct Index element  index such that A[index.x, index.y] is the maximum element in A which is not on the diagonal
Index off_diagonal_index (Matrix* A){
    //Initialize variables need later
    Index idx;
    idx.x = 0;
    idx.y = 1;
    int n = A->n;
    int i;
    int j;
    double num;
    double max_val = A->array[idx.x][idx.y];

    for (i=0; i<n; i++){
        for (j=i+1; j<n; j++){
            num = A->array[i][j];
            if (fabs(num)>fabs(max_val)){
                max_val = num;
                idx.x = i;
                idx.y = j;
            }
        }
    }
    return idx;
}

//returns a Matrix P built as described in Project's Description
Matrix build_P (Matrix* A){
    int n = A->n;
    Matrix P = mat_identity(n);
    double** arr = A->array;

    //find the index of the largest off diagonal element
    Index idx = off_diagonal_index(A);
    int i = idx.x;
    int j = idx.y;


    //Calculate Values of theta, c, s, t
    double theta = (arr[j][j] - arr[i][i])/(2.0 *arr[i][j]);
    double t = sign(theta)/ (fabs(theta) + sqrt(theta*theta +1));
    double c = 1 / sqrt(t*t+1);
    double s = t*c;

    //Update Values in the correct places in P
    P.array[i][i] = c;
    P.array[j][j] = c;
    P.array[i][j] = s;
    P.array[j][i] = -s;

    return P;

}

//Return the sum of squares of all non diagonal entries
double off (Matrix* A){
    //Initialize Parameters needed later
    double sum = 0;
    int i,j;
    int n = A->n;

    //Sum all squares of non diagonal entries
    for (i=0; i<n; i++){
        for (j=0; j<n; j++){
            if (i!=j){
                sum += pow(A->array[i][j],2);
            }
        }
    }
    return sum;
}

//Struct to store eigenvalues and eigenvectors output from Jacobi Algorithm
typedef struct Eigen{
    double* eigvals;
    Matrix eigvects;
} Eigen;

//A function to depp free the Eigen Struct
void free_eigen(Eigen* eigen){
    free_mat(&eigen->eigvects);
    free(eigen->eigvals);
}

//Calculate A_prime according to the definition in project based on indices (not matrix multiplication)
Matrix calc_A_prime(Matrix* A){
    //Initialize variables
    double** arr = A->array;
    int n = A->n;
    Matrix A_prime = mat_copy(A);
    int r;

    //find the index of the largest off diagonal element
    Index idx = off_diagonal_index(A);
    int i = idx.x;
    int j = idx.y;


    //Calculate Values of theta, c, s, t
    double theta = (arr[j][j] - arr[i][i])/(2.0 *arr[i][j]);
    double t = sign(theta)/ (fabs(theta) + sqrt(theta*theta +1));
    double c = 1 / sqrt(t*t+1);
    double s = t*c;

    for (r =0; r<n; r++){
        if (r!=i && r!=j){
            A_prime.array[r][i] = c * arr[r][i] - s * arr[r][j];
            A_prime.array[r][j] = c * arr[r][j] + s* arr [r][i];
        }
    }
    A_prime.array[i][i] = c*c*arr[i][i] + s*s*arr[j][j] - 2*s*c*arr[i][j];
    A_prime.array[j][j] = s*s*arr[i][i] + c*c*arr[j][j] + 2*s*c*arr[i][j];
    A_prime.array[i][j] = 0;
    A_prime.array[j][i] = 0;

    return A_prime;
}


Eigen jacobi_algorithm(Matrix* mat){

    int n = mat->n;
    int i;
    //Calculate A, A', P, V (Initialized)
    Matrix A = mat_copy(mat); //Create a Copy so won't destroy original
    Matrix P = build_P(&A);
    Matrix A_prime = calc_A_prime(&A);
    Matrix V = mat_copy(&P); //in the future V =P1 @P2 @ P3.....
    Matrix V_temp;
    double diff = off(&A) - off(&A_prime);
    int c =0;


    while (diff> EPSILON && c<100){


        //printf("Jacobi Iteration number %d.| diff is %f|off(A') =  %f\n",c,diff,off(&A_prime));
        //Free What needs to be Free
        free_mat(&A);
        free_mat(&P);

        //Build A_prime and P
        A = A_prime;
        P = build_P(&A);
        A_prime = calc_A_prime(&A);

        //Update Eigenvectors Matrix
        V_temp = mat_mul(&V,&P); //A temporary variable so V can be freed later
        free_mat(&V);
        V = V_temp;

        diff = off(&A) - off(&A_prime);
        c++;
    }

    //return Value as Struct - eigvals (array of doubles) and eigvecs (Matrix)
    Eigen eigen;
    eigen.eigvals = calloc(n,sizeof (double));
    for (i=0;i<n;i++){
        eigen.eigvals[i] = A_prime.array[i][i];
    }
    eigen.eigvects = V;

    return eigen;

}

//*********************************************************************************//
//******************************* **  Tests *************************************//
//*********************************************************************************//

//Tests for Matrix Operations
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

//Tests for algorithm (need to compare with python implemantation

void mat_to_file(Matrix* mat, char* file_name){
    FILE* fptr;
    char* pre = "C:\\Users\\Tomer\\CLionProjects\\spectral_k_means\\tests\\";
    char name [1024] = {0};
    strcat(name,pre);
    strcat(name,file_name);

    fptr = fopen(name,"w");
    int i, j;
    int n= mat->n;
    int m= mat->m;

    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            fprintf(fptr,"%f", mat->array[i][j]);
            if (j != m - 1) {
                fprintf(fptr,",");
            }
        }
        fprintf(fptr , "\n");
    }
    fclose(fptr);
}
void array_to_file(double * mat,int n ,char* file_name){
    FILE* fptr;
    char* pre = "C:\\Users\\Tomer\\CLionProjects\\spectral_k_means\\tests\\";
    char name [1024] = {0};
    strcat(name,pre);
    strcat(name,file_name);

    fptr = fopen(name,"w");
    int i, j;

    for (i = 0; i < n; i++) {
            fprintf(fptr,"%f",mat[i]);
            if (i != n - 1) {
                fprintf(fptr,",");
        }
    }
    fclose(fptr);
}
void test_WDH(Data* data){
    Matrix W  = build_W(data);
    mat_to_file(&W,"out_w.txt");
    Matrix D  = build_D(&W);
    mat_to_file(&D,"out_D.txt");
    Matrix D_half  = build_D_half(&D);
    mat_to_file(&D_half,"out_D_half.txt");
    Matrix l_norm = laplacian(&D_half,&W);
    mat_to_file(&l_norm,"out_lap.txt");
    Matrix P = build_P(&l_norm);
    mat_to_file(&P,"out_p.txt");
    Eigen eigen = jacobi_algorithm(&l_norm);
    mat_to_file(&eigen.eigvects,"out_eigvects.txt");
    array_to_file(eigen.eigvals,eigen.eigvects.n,"out_eigvals.txt");




    /*free_mat(&W);
    free_mat(&D);
    free_mat(&D_half);
    free_mat(&l_norm);
    free_mat(&P);
    free_eigen(&eigen);*/
}
void test_off_diag (){
    Matrix mat = zeros(5,5);
    int i;
    int j;
    int t=1;
    int c =1;
    int n = mat.n;

    for (i=0; i<n; i++){
        for (j=i; j<n; j++){
            mat.array [i][j] = t*c*sqrt(5);
            mat.array [j][i] = t*c*sqrt(5);
            t = t*-1;
            c++;
        }
    }
    printf("This is Mat: \n");
    print_mat(&mat);
    Index idx = off_diagonal_index(&mat);
    printf("largest value is in index (%d,%d): \n", idx.x, idx.y);
}





//********************************************


int main(){
    Data data = load_data("C:\\Users\\Tomer\\CLionProjects\\spectral_k_means\\tests\\input.txt");
    test_WDH(&data);
    //test_off_diag();


}
