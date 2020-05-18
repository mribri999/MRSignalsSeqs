#include "op_slewrate.h"

/**
 * Initialize the opD struct
 * This is the operator that enforces the global slew rate limit
 */
void cvxop_slewrate_init(cvxop_slewrate *opD, int N, double dt, double smax, double init_weight, double regularize, int verbose)
{
    opD->active = 1;
    opD->N = N;
    opD->dt = dt;

    opD->weight = init_weight;
    opD->base_smax = smax*dt;
    opD->smax = smax*dt*opD->weight;
    opD->sigD = 1.0/(2.0*opD->weight);

    opD->regularize = regularize;

    opD->verbose = verbose;

    if (opD->verbose>0) {   
        printf("opD->smax = %.2e\n", opD->smax);
    }

    cvxmat_alloc(&opD->zD, N-1, 1);
    cvxmat_alloc(&opD->zDbuff, N-1, 1);
    cvxmat_alloc(&opD->zDbar, N-1, 1);
    cvxmat_alloc(&opD->Dx, N-1, 1);
}

/**
 * Reweight the constraint and update all the subsequent weightings, and also the current descent direction zB
 */
void cvxop_slewrate_reweight(cvxop_slewrate *opD, double weight_mod)
{
    if (opD->weight < 1.0e64) { // prevent overflow
        opD->weight *= weight_mod;
        opD->smax = opD->base_smax*opD->weight;
        opD->sigD = 1.0/(2.0*opD->weight);

        // zD is usually reset anyways, but this is needed to maintain the existing relaxation
        for (int i = 0; i < opD->zD.N; i++) {
            opD->zD.vals[i] *= weight_mod;
        }
    }
}


/**
 * Add absolute volume of columns to tau
 */
void cvxop_slewrate_add2tau(cvxop_slewrate *opD, cvx_mat *tau_mat)
{
    int N = tau_mat->rows;
    for (int i = 0; i < N; i++) {
        if ((i == 0) || (i == N-1)) {
            tau_mat->vals[i] += opD->weight;
        } else {
            tau_mat->vals[i] += 2.0*opD->weight;
        }
    }
}


/**
 * Step the gradient waveform (taumx)
 */
void cvxop_slewrate_add2taumx(cvxop_slewrate *opD, cvx_mat *taumx)
{   
    
    // MATH: taumx += D*zD
    taumx->vals[0] += -opD->weight*opD->zD.vals[0];
    for (int i = 1; i < opD->zD.N; i++) {
        taumx->vals[i] += opD->weight*(opD->zD.vals[i-1] - opD->zD.vals[i]);
    }
    taumx->vals[taumx->N-1] += opD->weight*opD->zD.vals[opD->zD.N-1];

}

/**
 * Primal dual update
 */
void cvxop_slewrate_update(cvxop_slewrate *opD, cvx_mat *txmx, double rr)
{
    //MATH: Dx = sigD * D * txmx
    cvxmat_setvals(&(opD->Dx), 0.0);
    for (int i = 0; i < opD->Dx.N; i++) {
        opD->Dx.vals[i] += opD->weight*opD->sigD*(txmx->vals[i+1] - txmx->vals[i]);
    }

    // MATH: zDbuff  = zD + Dx = zD + sigD.*(D*txmx) 
    for (int i = 0; i < opD->zDbuff.N; i++) {
        opD->zDbuff.vals[i] = opD->zD.vals[i] + opD->Dx.vals[i];
    }

    
    // MATH: zDbar = clip( zDbuff/sigD , [-smax, smax])
    for (int i = 0; i < opD->zDbar.N; i++) {
        double temp = opD->zDbuff.vals[i] / opD->sigD;
        if (temp > 0.99*opD->smax) {
            opD->zDbar.vals[i] = 0.99*opD->smax;
        } else if (temp < -0.99*opD->smax) {
            opD->zDbar.vals[i] = -0.99*opD->smax;
        } else {
            opD->zDbar.vals[i] = opD->regularize * temp; // This adds a minor smoothing to the descent
        }
    }

    // MATH: zDbar = zDbuff - sigD*zDbar
    for (int i = 0; i < opD->zDbar.N; i++) {
        opD->zDbar.vals[i] = opD->zDbuff.vals[i] - opD->sigD*opD->zDbar.vals[i];
    }

    // MATH: zD = rr*zDbar + (1-rr)*zD
    for (int i = 0; i < opD->zD.N; i++) {
        opD->zD.vals[i] = rr * opD->zDbar.vals[i] + (1.0 - rr) * opD->zD.vals[i];
    }
}

/**
 * Check if slew rates are under smax
 */
int cvxop_slewrate_check(cvxop_slewrate *opD, cvx_mat *G)
{
    cvxmat_setvals(&(opD->Dx), 0.0);
    for (int i = 0; i < opD->Dx.N; i++) {
        opD->Dx.vals[i] += opD->weight*opD->sigD*(G->vals[i+1] - G->vals[i]);
    }

    int slew_too_high = 0;
    int slew_bad = 0;

    for (int i = 0; i < opD->Dx.N; i++) {
        if (fabs(opD->Dx.vals[i]) > opD->smax*opD->sigD) {
            slew_too_high++;
            slew_bad = 1;
        }
    }

    if (opD->verbose>0) {  
        printf("    Slew check:     (%d)  [%.2e]  %d \n", slew_bad, opD->weight, slew_too_high);
    }

    return slew_bad;
}

/*
 * Free arrays
 */
void cvxop_slewrate_destroy(cvxop_slewrate *opD)
{
    free(opD->zD.vals);
    free(opD->zDbuff.vals);
    free(opD->zDbar.vals);
    free(opD->Dx.vals);
}