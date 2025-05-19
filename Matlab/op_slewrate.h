#ifndef CVX_OPSLEWRATE_H
#define CVX_OPSLEWRATE_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "cvx_matrix.h"

typedef struct {
    int active;
    int verbose;

    int N;
    int Naxis;
    int Ntotal;
    double dt;

    double regularize;

    double smax;
    double base_smax;
    double weight;
    double sigD;

    cvx_mat zD;
    cvx_mat zDbuff;
    cvx_mat zDbar;
    cvx_mat Dx;

} cvxop_slewrate;

void cvxop_slewrate_init(cvxop_slewrate *opD, int N, int Naxis, double dt, double smax, 
                  double init_weight, double regularize, int verbose);
void cvxop_slewrate_reweight(cvxop_slewrate *opD, double weight_mod);
void cvxop_slewrate_add2tau(cvxop_slewrate *opD, cvx_mat *tau_mat);
void cvxop_slewrate_add2taumx(cvxop_slewrate *opD, cvx_mat *taumx);
void cvxop_slewrate_update(cvxop_slewrate *opD, cvx_mat *txmx, double relax);
int cvxop_slewrate_check(cvxop_slewrate *opD, cvx_mat *G);
void cvxop_slewrate_destroy(cvxop_slewrate *opD);


#endif /* CVX_OPSLEWRATE_H */