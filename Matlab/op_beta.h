#ifndef CVX_OPBETA_H
#define CVX_OPBETA_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "cvx_matrix.h"

typedef struct {
    int active;
    int verbose;

    int N;
    double dt;
    double weight;

    cvx_mat C;
} cvxop_beta;

void cvxop_beta_init(cvxop_beta *opC, int N, double dt, double weight, int verbose);
void cvxop_beta_reweight(cvxop_beta *opC, double weight_mod);
void cvxop_beta_add2taumx(cvxop_beta *opC, cvx_mat *taumx);
void cvxop_beta_destroy(cvxop_beta *opC);
void cvxop_beta_add2tau(cvxop_beta *opC, cvx_mat *tau_mat);


#endif /* CVX_OPBETA_H */