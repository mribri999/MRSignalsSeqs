#ifndef CVX_OPGRADIENT_H
#define CVX_OPGRADIENT_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "cvx_matrix.h"

typedef struct {
    int active;
    int verbose;

    int N;
    double dt;

    double gmax;
    int ind_inv;

    cvx_mat Gfix;

} cvxop_gradient;

void cvxop_gradient_init(cvxop_gradient *opG, int N, double dt, double gmax, int ind_inv, int verbose);
void cvxop_gradient_limiter(cvxop_gradient *opG, cvx_mat *xbar);
void cvxop_gradient_setFixRange(cvxop_gradient *opG, int start, int end, double val);
void cvxop_init_G(cvxop_gradient *opG, cvx_mat *G);
int cvxop_gradient_check(cvxop_gradient *opG, cvx_mat *G);
double cvxop_gradient_getbval(cvxop_gradient *opG, cvx_mat *G);
int cvxop_gradient_overflowcheck(cvxop_gradient *opG, cvx_mat *G);
void cvxop_gradient_setGfix(cvxop_gradient *opG, int N_gfix, double *gfix) ;
void cvxop_gradient_destroy(cvxop_gradient *opG);



#endif /* CVX_OPGRADIENT_H */