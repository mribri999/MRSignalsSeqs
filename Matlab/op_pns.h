#ifndef CVX_OPPNS_H
#define CVX_OPPNS_H

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
    double thresh;
    double weight;

    int n_conv;

    cvx_mat coeff; 
    cvx_mat Px;

    cvx_mat zP;
    cvx_mat zPbar;
    cvx_mat zPbuff;
    cvx_mat Ptau;


} cvxop_pns;

void cvxop_pns_init(cvxop_pns *opP, int N, double dt, int ind_inv, double thresh, double init_weight, int verbose);
int cvxop_pns_check(cvxop_pns *opP, cvx_mat *G);
void cvxop_pns_add2tau(cvxop_pns *opP, cvx_mat *tau_mat);
void cvxop_pns_add2taumx(cvxop_pns *opP, cvx_mat *taumx);
void cvxop_pns_update(cvxop_pns *opP, cvx_mat *txmx, double rr);
int cvxop_pns_check(cvxop_pns *opP, cvx_mat *G);
void cvxop_pns_reweight(cvxop_pns *opP, double weight_mod);




#endif /* CVX_OPPNS_H */