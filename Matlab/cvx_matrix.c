#include <stdlib.h>
#include <stdio.h>

#include "cvx_matrix.h"

int cvxmat_alloc(cvx_mat *mat, int rows, int cols)
{
    mat->vals = malloc(rows * cols * sizeof(double));
    if (mat->vals == NULL) {
        fprintf(stderr, "*** ERROR: Out-of-memory rows = %d  cols = %d ***\n", rows, cols);
        return 0;
    }

    mat->rows = rows;
    mat->cols = cols;
    mat->N = rows*cols;

    // Slows us slightly, but much safer
    for (int i = 0; i < mat->N; i++) {
        mat->vals[i] = 0.0;
    }

    return 1;
}

double cvxmat_get(cvx_mat *mat, int row, int col) 
{
    return mat->vals[mat->cols * row + col];
}

void cvxmat_set(cvx_mat *mat, int row, int col, double val) 
{
    mat->vals[mat->cols * row + col] = val;
}

void cvxmat_EWinvert(cvx_mat *mat) 
{
    for (int i = 0; i < mat->N; i++) {
        mat->vals[i] = 1.0/mat->vals[i];
    }
}

void cvxmat_EWmultIP(cvx_mat *mat0, cvx_mat *mat1)
{
    for (int i = 0; i < mat0->N; i++) {
        mat0->vals[i] *= mat1->vals[i];
    }
}

void cvxmat_subractMat(cvx_mat *mat0, cvx_mat *mat1, cvx_mat *mat2)
{
    for (int i = 0; i < mat0->N; i++) {
        mat0->vals[i] = mat1->vals[i] - mat2->vals[i];
    }
}

void cvxmat_multMatMat(cvx_mat *C, cvx_mat *A, cvx_mat *B)
{
    double sum;
    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < B->cols; j++) {
            sum = 0.0;
            for (int k = 0; k < A->cols; k++) {
                //sum += A[i,k] * B[k,j]
                sum += A->vals[k + A->cols * i] * B->vals[j + B->cols * k];
            }
            C->vals[j + C->cols * i] = sum;
        }
    }
}

void cvxmat_inplace_transpose(cvx_mat *A)
{
    double *new_vals;
    double *temp_free = A->vals;
    new_vals = malloc(A->N * sizeof(double));
    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->cols; j++) {
            new_vals[A->cols * i + j] = A->vals[A->cols * j + i];
        }
    }
    A->vals = new_vals;
    free(temp_free);

}

void cvxmat_multAtA(cvx_mat *C, cvx_mat *A)
{
    cvxmat_inplace_transpose(A);
    double sum;
    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->cols; j++) {
            sum = 0.0;
            for (int k = 0; k < A->cols; k++) {
                //sum += A[i,k] * B[k,j]
                sum += A->vals[A->cols * i + k] * A->vals[A->cols * j + k];
                // sum += A->vals[A->cols * k + i] * A->vals[A->cols * k + j];
            }
            C->vals[j + C->cols * i] = sum;
        }
    }
}



void cvxmat_multAx(cvx_mat *b, cvx_mat *A, cvx_mat *x)
{
    double sum;
    for (int i = 0; i < A->rows; i++) {
        sum = 0.0;
        for (int j = 0; j < A->cols; j++) {
            sum += A->vals[A->cols * i + j] * x->vals[j];
        }
        b->vals[i] = sum;
    }
}

void cvxmat_multAx2(cvx_mat *b, cvx_mat *A, cvx_mat *x)
{
    double sum;
    double *Apos = &(A->vals[0]);
    double *bpos = &(b->vals[0]);
    double *xpos;
    int cols = A->cols;
    for (int i = 0; i < A->rows; i++) {
        xpos = &(x->vals[0]);
        sum = 0.0;
        for (int j = 0; j < cols; j++) {
            sum += (*Apos++) * (*xpos++);
        }
        *bpos = sum;
        *bpos++;
    }
}


void cvxmat_multAx22(cvx_mat *b, cvx_mat *A, cvx_mat *x)
{
	int size = A->rows;

    double *Apos1 = &(A->vals[0]);
    double *Apos2 = &(A->vals[size]);

    double *ypos = &(b->vals[0]);
    
    for(int i=0; i<size/2; i++)
    {

        double ytemp1 = 0;
        double ytemp2 = 0;
        double *xpos = &(x->vals[0]);

        for(int j=0; j<size; j++)
        {

            ytemp1 += (*Apos1++) * (*xpos);
            ytemp2 += (*Apos2++) * (*xpos);

            xpos++;
        }
        *ypos = ytemp1;
        ypos++;
        *ypos = ytemp2;

        ypos++;
        
        // skip next row
        Apos1 += size;
        Apos2 += size;
    }

}


void cvxmat_multAx23(cvx_mat *b, cvx_mat *A, cvx_mat *x)
{
	int size = A->rows;

    double *Apos1 = &(A->vals[0]);
    double *Apos2 = &(A->vals[size]);

    double *ypos = &(b->vals[0]);
    
    for(int i=0; i<size/2; i++)
    {

        register double ytemp1 = 0;

        register double ytemp2 = 0;
        register double *xpos = &(x->vals[0]);

        int blocksize = size / 8;

        double x0 = xpos[0];
        double x1 = xpos[1];
        
        for(int j=0; j<blocksize; j++)
        {

            ytemp1 += x0 * (Apos1[0]);
            ytemp2 += x0 * (Apos2[0]);

            x0 = xpos[2];
            ytemp1 += x1 * (Apos1[1]);
            ytemp2 += x1 * (Apos2[1]);
            x1 = xpos[3];

            ytemp1 += x0 * (Apos1[2]);
            ytemp2 += x0 * (Apos2[2]);

            x0 = xpos[4];
            ytemp1 += x1 * (Apos1[3]);
            ytemp2 += x1 * (Apos2[3]);
            x1 = xpos[5];

            ytemp1 += x0 * (Apos1[4]);
            ytemp2 += x0 * (Apos2[4]);

            x0 = xpos[6];
            ytemp1 += x1 * (Apos1[5]);
            ytemp2 += x1 * (Apos2[5]);
            x1 = xpos[7];

            xpos+=8;
            ytemp1 += x0 * (Apos1[6]);
            ytemp2 += x0 * (Apos2[6]);
            x0 = xpos[0];

            ytemp1 += x1 * (Apos1[7]);
            Apos1+=8;

            ytemp2 += x1 * (Apos2[7]);
            x1 = xpos[1];

            Apos2+=8;

            
        }
        
        ypos[0] = ytemp1;
        ypos[1] = ytemp2;

        ypos+=2;
        
        // skip next row
        Apos1 += size;
        Apos2 += size;
    }

}







void cvxmat_multAx4(cvx_mat *b, cvx_mat *A, cvx_mat *x)
{
	int size = A->rows;

    register double ytemp1 = 0;

    register double ytemp2 = 0;
    register double *xpos = &(x->vals[0]);

    int blocksize = size / 8;

    double x0 = xpos[0];

    double x1 = xpos[1];

    double *ypos = &(b->vals[0]);

    double *Apos1 = &(A->vals[0]);
    double *Apos2 = &(A->vals[size]);
    
    for(int j=0; j<blocksize; j++)
    {

        ytemp1 += x0 * (Apos1[0]);
        ytemp2 += x0 * (Apos2[0]);

        x0 = xpos[2];
        ytemp1 += x1 * (Apos1[1]);
        ytemp2 += x1 * (Apos2[1]);
        x1 = xpos[3];

        ytemp1 += x0 * (Apos1[2]);
        ytemp2 += x0 * (Apos2[2]);

        x0 = xpos[4];
        ytemp1 += x1 * (Apos1[3]);
        ytemp2 += x1 * (Apos2[3]);
        x1 = xpos[5];

        ytemp1 += x0 * (Apos1[4]);
        ytemp2 += x0 * (Apos2[4]);

        x0 = xpos[6];
        ytemp1 += x1 * (Apos1[5]);
        ytemp2 += x1 * (Apos2[5]);
        x1 = xpos[7];

        xpos+=8;
        ytemp1 += x0 * (Apos1[6]);
        ytemp2 += x0 * (Apos2[6]);
        x0 = xpos[0];

        ytemp1 += x1 * (Apos1[7]);
        Apos1+=8;

        ytemp2 += x1 * (Apos2[7]);
        x1 = xpos[1];

        Apos2+=8;

        
    }
    
    ypos[0] = ytemp1;
    ypos[1] = ytemp2;

    ypos+=2;
    
    // skip next row
    Apos1 += size;
    Apos2 += size;

}







			




void cvxmat_multAtx(cvx_mat *b, cvx_mat *A, cvx_mat *x)
{
    double sum;
    for (int i = 0; i < A->cols; i++) {
        sum = 0.0;
        for (int j = 0; j < A->rows; j++) {
            sum += A->vals[A->rows * j + i] * x->vals[j];
        }
        b->vals[i] = sum;
    }
}



void cvxmat_subractMatMult1(cvx_mat *mat0, double v0, cvx_mat *mat1, cvx_mat *mat2)
{
    for (int i = 0; i < mat0->N; i++) {
        mat0->vals[i] = v0 * mat1->vals[i] - mat2->vals[i];
    }
}

void cvxmat_updateG(cvx_mat *G, double relax, cvx_mat *xbar)
{
    // G=p*xbar+(1-p)*G;
    for (int i = 0; i < G->N; i++) {
        G->vals[i] = relax * xbar->vals[i] + (1.0-relax) * G->vals[i];
    }
}

void cvxmat_setvals(cvx_mat *mat, double val)
{   
    for (int i = 0; i < mat->N; i++) {
        mat->vals[i] = val;
    }
}

int copyNewMatrix(cvx_mat *mat_in, cvx_mat *mat_out)
{
    cvxmat_alloc(mat_out, mat_in->rows, mat_in->cols);

    int N = mat_out->rows * mat_out->cols;
    for (int i = 0; i < N; i++) {
        mat_out->vals[i] = mat_in->vals[i];
    }

    return 1;
}


/*
int main (void)
{
    printf ("In matrix.c main function\n");
    
    cvx_mat mat;
    allocMatrix(&mat, 100, 100);
    printf ("Matrix rows = %d, cols = %d\n", mat.rows, mat.cols);

    cvx_mat mat2;
    allocMatrix(&mat2, 1000000000, 1000000000);
    printf ("Matrix rows = %d, cols = %d\n", mat2.rows, mat2.cols);

    return 0;
}
*/