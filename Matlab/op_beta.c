#include "op_beta.h"

/**
 * Initialize the opC struct
 * This is the operator that maximizes beta-value for diffusion purposes
 * At this point it is mostly superceded by bval.c and this won't be tested much anymore
 */
void cvxop_beta_init(cvxop_beta *opC, int N, double dt, double weight, int verbose) {
    opC->N = N;
    opC->dt = dt;
    opC->verbose = verbose;
    opC->weight = weight;

    cvxmat_alloc(&opC->C, N, 1);

    double tt;
    for (int i = 0; i < N; i++) {
        tt = N-i;
        opC->C.vals[i] = tt*(tt+1)/2.0;
    }

    double norm = 0.0;

    for (int i = 0; i < N; i++) {
        norm += opC->C.vals[i] * opC->C.vals[i];
    }

    norm = sqrt(norm);

    for (int i = 0; i < N; i++) {
        opC->C.vals[i] *= opC->weight / norm;
    }
}

/**
 * Reweight the C matrix (weight_mod * C)
 */
void cvxop_beta_reweight(cvxop_beta *opC, double weight_mod)
{
    opC->weight *= weight_mod;

    for (int i = 0; i < opC->C.N; i++) {
        opC->C.vals[i] *= opC->weight;
    }
}

/**
 * Add absolute value of columns to the tau matrix 
 */
void cvxop_beta_add2tau(cvxop_beta *opC, cvx_mat *tau_mat)
{
    if (opC->active > 0) {
        for (int i = 0; i < opC->N; i++) {
            tau_mat->vals[i] += fabs(opC->C.vals[i]);
        }
    }
}


/**
 * Step the gradient waveform (taumx)
 */
void cvxop_beta_add2taumx(cvxop_beta *opC, cvx_mat *taumx)
{
    if (opC->active > 0) {
        // MATH: taumx -= C
        for (int i = 1; i < taumx->N; i++) {
            taumx->vals[i] -= opC->C.vals[i];
        }

    }
}

void cvxop_beta_destroy(cvxop_beta *opC)
{
    free(opC->C.vals);
}