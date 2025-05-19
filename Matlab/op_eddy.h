#ifndef CVX_OPEDDY_H
#define CVX_OPEDDY_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "cvx_matrix.h"

typedef struct {
    int active;
    int verbose;

    int N;
    int Nrows;

    int ind_inv;
    double dt;

    double weight;

    cvx_mat E0;
    cvx_mat E;

    cvx_mat norms;
    cvx_mat weights;
    cvx_mat checks;
    cvx_mat tolerances;
    cvx_mat goals;


    cvx_mat sigE;
    cvx_mat zE;
    cvx_mat zEbuff;
    cvx_mat zEbar;
    cvx_mat Ex;


} cvxop_eddy;

#define MAXROWS 16

void cvxop_eddy_init(cvxop_eddy *opE, int N, int ind_inv, double dt,
                     double init_weight, int verbose);
void cvxop_eddy_addrow(cvxop_eddy *opE, double lambda, double goal, double tol, double offset);
void cvxop_eddy_finishinit(cvxop_eddy *opE);

void cvxop_eddy_add2tau(cvxop_eddy *opE, cvx_mat *tau_mat);
void cvxop_eddy_add2taumx(cvxop_eddy *opE, cvx_mat *taumx);
void cvxop_eddy_update(cvxop_eddy *opE, cvx_mat *txmx, double relax);
int cvxop_eddy_check(cvxop_eddy *opE, cvx_mat *G);
void cvxop_eddy_reweight(cvxop_eddy *opE, double weight_mod);
void cvxop_eddy_destroy(cvxop_eddy *opE);

#endif /* CVX_OPEDDY_H */