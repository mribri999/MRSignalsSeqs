#ifndef CVX_MATRIX_H
#define CVX_MATRIX_H

typedef struct {
    int rows;
    int cols;
    int N;
    double *vals;
} cvx_mat;

int cvxmat_alloc(cvx_mat *mat, int rows, int cols);
void cvxmat_setvals(cvx_mat *mat, double val);


double cvxmat_get(cvx_mat *mat, int i, int j);
void cvxmat_set(cvx_mat *mat, int i, int j, double val);

void cvxmat_updateG(cvx_mat *G, double relax, cvx_mat *xbar);
void cvxmat_multMatMat(cvx_mat *C, cvx_mat *A, cvx_mat *B);
void cvxmat_multAtA(cvx_mat *C, cvx_mat *A);
void cvxmat_multAx(cvx_mat *b, cvx_mat *A, cvx_mat *x);
void cvxmat_multAx2(cvx_mat *b, cvx_mat *A, cvx_mat *x);
void cvxmat_multAtx(cvx_mat *b, cvx_mat *A, cvx_mat *x);

int copyNewMatrix(cvx_mat *mat_in, cvx_mat *mat_out);


void cvxmat_EWinvert(cvx_mat *mat);
void cvxmat_EWmultIP(cvx_mat *mat0, cvx_mat *mat1);
void cvxmat_subractMat(cvx_mat *mat0, cvx_mat *mat1, cvx_mat *mat2);
void cvxmat_subractMatMult1(cvx_mat *mat0, double v0, cvx_mat *mat1, cvx_mat *mat2);

#endif /* CVX_MATRIX_H */