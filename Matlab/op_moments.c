#include "op_moments.h"

/**
 * Initialize the opB struct
 * This is the operator that sets
 */
void cvxop_moments_init(cvxop_moments *opQ, int N, int ind_inv, double dt,
                        double *moment_tols_in, double init_weight, int verbose)
{
    opQ->active = 1;
    opQ->N = N;
    opQ->ind_inv = ind_inv;
    opQ->dt = dt;
    opQ->verbose = verbose;
    opQ->weight = init_weight;

    opQ->Nm = 3; // Hard code maximum number of moments for now

    cvxmat_alloc(&opQ->moment_tol0, opQ->Nm, 1);
    cvxmat_alloc(&opQ->Q0, opQ->Nm, N);

    cvxmat_alloc(&opQ->norms, opQ->Nm, 1);

    cvxmat_alloc(&opQ->Q, opQ->Nm, N);
    cvxmat_alloc(&opQ->moment_tol, opQ->Nm, 1);

    cvxmat_alloc(&opQ->sigQ, opQ->Nm, 1);


    
    cvxmat_alloc(&opQ->zQ, opQ->Nm, 1);
    cvxmat_alloc(&opQ->zQbuff, opQ->Nm, 1);
    cvxmat_alloc(&opQ->zQbar, opQ->Nm, 1);
    cvxmat_alloc(&opQ->Qx, opQ->Nm, 1);



    // Copy initial moment tolerances into cvx_mat array
    for (int j = 0; j < opQ->Nm; j++) {
        opQ->moment_tol0.vals[j] = moment_tols_in[j];
    }

    // Set Q0 array rows to dt^momentnum
    for (int j = 0; j < opQ->Nm; j++) {
        for (int i = 0; i < N; i++) {
            double temp = pow((dt*i), (double)j);
            cvxmat_set(&(opQ->Q0), j, i, temp);
        }
    }
    
    // Scale so that Q0 returns true moments (in T*s/m)
    for (int j = 0; j < opQ->Nm; j++) {
        for (int i = 0; i < N; i++) {
            double temp = cvxmat_get(&(opQ->Q0), j, i) * dt;
            cvxmat_set(&(opQ->Q0), j, i, temp);
        }
    }
    
    // Gradient values need to be reversed after the 180
    for (int j = 0; j < opQ->Nm; j++) {
        for (int i = ind_inv; i < N; i++) {
            double temp = cvxmat_get(&(opQ->Q0), j, i) * -1.0;
            cvxmat_set(&(opQ->Q0), j, i, temp);
        }
    }


    // Calculate the row norms of the moment array and store
    for (int j = 0; j < opQ->Nm; j++) {
        for (int i = 0; i < N; i++) {
            double temp = cvxmat_get(&(opQ->Q0), j, i);
            opQ->norms.vals[j] += temp*temp;
        }
        opQ->norms.vals[j] = sqrt(opQ->norms.vals[j]);
    }

    // Scale the Q0 array and copy into Q
    // Q = Q0 * weight/norms
    // also scale moment_tol as scaled moment_tol0
    for (int j = 0; j < opQ->Nm; j++) {
        for (int i = 0; i < N; i++) {
            double temp = cvxmat_get(&(opQ->Q0), j, i);
            cvxmat_set(&(opQ->Q), j, i, opQ->weight * temp / opQ->norms.vals[j]);
        }
        opQ->moment_tol.vals[j] = opQ->weight *  opQ->moment_tol0.vals[j] / opQ->norms.vals[j];
    }


    if (opQ->verbose>0) {   
        printf("Q norms = %.2e  %.2e  %.2e\n", opQ->norms.vals[0], opQ->norms.vals[1], opQ->norms.vals[2]);
    }

    // Calculate sigQ as inverse of sum(abs(row of Q))
    for (int j = 0; j < opQ->Nm; j++) {
        for (int i = 0; i < N; i++) {
            double temp = cvxmat_get(&(opQ->Q), j, i);
            opQ->sigQ.vals[j] += fabs(temp);
        }
        opQ->sigQ.vals[j] = 1.0 / opQ->sigQ.vals[j];
    }
    
    if (opQ->verbose > 0) {   
        printf("sigQ = %.2e  %.2e  %.2e\n", opQ->sigQ.vals[0], opQ->sigQ.vals[1], opQ->sigQ.vals[2]);
    }
}

/*
 * Reweight the constraint and update all the subsequent weightings, and also the current descent direction zQ
 * basically weight_mod * Q
 */
void cvxop_moments_reweight(cvxop_moments *opQ, double weight_mod)
{
    opQ->weight *= weight_mod;

    // Scale the Q0 array and copy into Q
    for (int j = 0; j < opQ->Nm; j++) {
        for (int i = 0; i < opQ->N; i++) {
            double temp = cvxmat_get(&(opQ->Q0), j, i);
            cvxmat_set(&(opQ->Q), j, i, opQ->weight * temp / opQ->norms.vals[j]);
        }
        opQ->moment_tol.vals[j] = opQ->weight *  opQ->moment_tol0.vals[j] / opQ->norms.vals[j];
    }

    // Calculate sigQ as inverse of sum(abs(row of Q))
    for (int j = 0; j < opQ->Nm; j++) {
        opQ->sigQ.vals[j] = 0.0;
        for (int i = 0; i < opQ->N; i++) {
            double temp = cvxmat_get(&(opQ->Q), j, i);
            opQ->sigQ.vals[j] += fabs(temp);
        }
        opQ->sigQ.vals[j] = 1.0 / opQ->sigQ.vals[j];
    }

    // zQ is usually reset anyways, but this is needed to maintain the existing relaxation
    for (int i = 0; i < opQ->zQ.N; i++) {
        opQ->zQ.vals[i] *= weight_mod;
    }
}


/**
 * Add absolute value of columns to the tau matrix 
 */
void cvxop_moments_add2tau(cvxop_moments *opQ, cvx_mat *tau_mat)
{
    for (int j = 0; j < opQ->Nm; j++) {
        if (opQ->moment_tol.vals[j] >= 0) {
            for (int i = 0; i < opQ->N; i++) {
                double temp = cvxmat_get(&(opQ->Q), j, i);
                tau_mat->vals[i] += fabs(temp);
            }
        }
    }
}


/**
 * Step the gradient waveform (taumx)
 */
void cvxop_moments_add2taumx(cvxop_moments *opQ, cvx_mat *taumx)
{   
    // MATH: taumx += Q*zQ
    for (int j = 0; j < opQ->Nm; j++) {
        if (opQ->moment_tol.vals[j] >= 0) {
            for (int i = 0; i < opQ->N; i++) {
                double temp = cvxmat_get(&(opQ->Q), j, i);
                taumx->vals[i] += (temp * opQ->zQ.vals[j]);
            }
        }
    }

}

/**
 * Primal dual update
 */
void cvxop_moments_update(cvxop_moments *opQ, cvx_mat *txmx, double rr)
{

    cvxmat_setvals(&(opQ->Qx), 0.0);

    // MATH: Qx = Q * txmx
    for (int j = 0; j < opQ->Nm; j++) {
        for (int i = 0; i < opQ->N; i++) {
            double temp = cvxmat_get(&(opQ->Q), j, i) * txmx->vals[i];
            opQ->Qx.vals[j] += temp;
        }
    }

    // MATH: Qx = Qx * sigQ
    for (int j = 0; j < opQ->Nm; j++) {
        opQ->Qx.vals[j] *= opQ->sigQ.vals[j];
    }

    // MATH: zQbuff  = zQ + Qx = zQ + sigQ.*(Q*txmx) 
    for (int j = 0; j < opQ->Nm; j++) {
        opQ->zQbuff.vals[j] = opQ->zQ.vals[j] + opQ->Qx.vals[j];
    }

    // MATH: zQbar = clip( zQbuff/sigQ , [-moment_tol, moment_tol])
    for (int j = 0; j < opQ->Nm; j++) {
        double temp = opQ->zQbuff.vals[j] /  opQ->sigQ.vals[j];
        if (temp> 0.99*opQ->moment_tol.vals[j]) {
            opQ->zQbar.vals[j] = 0.99*opQ->moment_tol.vals[j];
        } else if (temp < -0.99*opQ->moment_tol.vals[j]) {
            opQ->zQbar.vals[j] = -0.99*opQ->moment_tol.vals[j];
        } else {
            opQ->zQbar.vals[j] = temp;
        }
    }

    // MATH: zQbar = zQbuff - sigQ*zQbar
    for (int j = 0; j < opQ->Nm; j++) {
        opQ->zQbar.vals[j] = opQ->zQbuff.vals[j] - opQ->sigQ.vals[j] * opQ->zQbar.vals[j];
    }

    // MATH: zQ = rr*zQbar + (1-rr)*zQ
    for (int i = 0; i < opQ->Nm; i++) {
        opQ->zQ.vals[i] = rr * opQ->zQbar.vals[i] + (1 - rr) * opQ->zQ.vals[i];
    }
}

/*
 * Check if moments are larger than a fixed tolerance
 */
int cvxop_moments_check(cvxop_moments *opQ, cvx_mat *G)
{
    cvxmat_setvals(&(opQ->Qx), 0.0);
    for (int j = 0; j < opQ->Nm; j++) {
        for (int i = 0; i < opQ->N; i++) {
            double temp = cvxmat_get(&(opQ->Q0), j, i) * G->vals[i];
            opQ->Qx.vals[j] += temp;
        }
    }

    int moments_bad = 0;
    double moment_tol = 1.0e-7;

    for (int j = 0; j < opQ->Nm; j++) {
        if (opQ->moment_tol0.vals[j] >= 0) {
            if (fabs(opQ->Qx.vals[j]) > moment_tol) {
                moments_bad = 1;
            }
        }
    }

    if (opQ->verbose>0) {   
        printf("    Moments check:  (%d)  %.2e  %.2e  %.2e \n", moments_bad, opQ->Qx.vals[0], opQ->Qx.vals[1], opQ->Qx.vals[2]);
    }

    return moments_bad;
}

/*
 * Free memory
 */
void cvxop_moments_destroy(cvxop_moments *opQ)
{
    free(opQ->moment_tol.vals);
    free(opQ->moment_tol0.vals);
    free(opQ->norms.vals);
    free(opQ->Q.vals);
    free(opQ->Q0.vals);
    free(opQ->sigQ.vals);
    free(opQ->zQ.vals);
    free(opQ->zQbuff.vals);
    free(opQ->zQbar.vals);
    free(opQ->Qx.vals);
}



