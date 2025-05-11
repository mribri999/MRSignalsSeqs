#ifndef CVX_OPMAXWELL_H
#define CVX_OPMAXWELL_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "cvx_matrix.h"

typedef struct {
    int active;
    int verbose;
    
    int N;
    int ind_inv;
    double dt;
    double tol;

    double sign;
    double sigX;
    double Xx;
    double zXbuff;
    double zXbar;
    double zX;

    cvx_mat x_store;

} cvxop_maxwell;

void cvxop_maxwell_init(cvxop_maxwell *opX, int N, double dt, int ind_inv, double tol, int verbose);
void cvxop_maxwell_update(cvxop_maxwell *opX, cvx_mat *txmx, double rr);
void cvxop_maxwell_add2taumx(cvxop_maxwell *opX, cvx_mat *taumx);
void cvxop_maxwell_add2tau(cvxop_maxwell *opX, cvx_mat *tau_mat);
int cvxop_maxwell_check(cvxop_maxwell *opX, cvx_mat *G);

#endif /* CVX_OPMAXWELL_H */