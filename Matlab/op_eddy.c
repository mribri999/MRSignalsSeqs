#include "op_eddy.h"

/**
 * Initialize the opE struct
 * This is the operator that nulls eddy currents
 */
void cvxop_eddy_init(cvxop_eddy *opE, int N, int ind_inv, double dt,
                     double init_weight, int verbose)
{
    opE->active = 1;
    opE->N = N;
    opE->ind_inv = ind_inv;
    opE->dt = dt;
    opE->verbose = verbose;

    opE->Nrows = 0; // # of eddy current rows

    cvxmat_alloc(&opE->E0, MAXROWS, N);
    cvxmat_alloc(&opE->E, MAXROWS, N);

    cvxmat_alloc(&opE->norms, MAXROWS, 1);
    cvxmat_alloc(&opE->weights, MAXROWS, 1);
    cvxmat_alloc(&opE->checks, MAXROWS, 1);
    cvxmat_alloc(&opE->tolerances, MAXROWS, 1);
    cvxmat_alloc(&opE->goals, MAXROWS, 1);
    cvxmat_alloc(&opE->sigE, MAXROWS, 1);

    cvxmat_alloc(&opE->zE, MAXROWS, 1);
    cvxmat_alloc(&opE->zEbuff, MAXROWS, 1);
    cvxmat_alloc(&opE->zEbar, MAXROWS, 1);
    cvxmat_alloc(&opE->Ex, MAXROWS, 1);

    cvxmat_setvals(&opE->weights, init_weight);
}


/**
 * Add a lambda constraint
 * lambda units are seconds (like dt)
 */
void cvxop_eddy_addrow(cvxop_eddy *opE, double lambda, double goal, double tol, double mode)
{
    
    if (mode < 0.5) {
        for (int i = 0; i < opE->N; i++) {
            double ii = opE->N - 1 - i; 
            double val;
            if (i == 0) {
                val = -exp(-ii*opE->dt/lambda);
            } else {
                val = exp(-(ii+1.0)*opE->dt/lambda) - exp(-ii*opE->dt/lambda);
            }
            cvxmat_set(&(opE->E0), opE->Nrows, i, -val);
        }
    } else {
        double val = 0.0;
        for (int i = 0; i < opE->N; i++) {
            double ii = opE->N - 1 - i;
            val += -exp(-ii*opE->dt/lambda);
        }
        cvxmat_set(&(opE->E0), opE->Nrows, 0, val*1.0e3*opE->dt);

        for (int i = 1; i < opE->N; i++) {
            double ii = opE->N - i;
            val = -exp(-ii*opE->dt/lambda);
            cvxmat_set(&(opE->E0), opE->Nrows, i, val*1.0e3*opE->dt);
        }
    }

    opE->tolerances.vals[opE->Nrows] = tol;
    opE->goals.vals[opE->Nrows] = goal;
    opE->Nrows += 1;
}


/**
 * Scale E to have unit norm rows, and calculate sigE
 */
void cvxop_eddy_finishinit(cvxop_eddy *opE)
{
    // Calculate the row norms of the eddy current array and store
    for (int j = 0; j < opE->Nrows; j++) {
        for (int i = 0; i < opE->N; i++) {
            double temp = cvxmat_get(&(opE->E0), j, i);
            opE->norms.vals[j] += temp*temp;
        }
        opE->norms.vals[j] = sqrt(opE->norms.vals[j]);
    }

    // Scale row norms to 1.0
    for (int j = 0; j < opE->Nrows; j++) {
        opE->weights.vals[j] /= opE->norms.vals[j];
    }

    // Make E as weighted E0
    for (int j = 0; j < opE->Nrows; j++) {
        for (int i = 0; i < opE->N; i++) {
            double temp = cvxmat_get(&(opE->E0), j, i);
            cvxmat_set(&(opE->E), j, i, opE->weights.vals[j] * temp);
        }
    }


    // Calculate sigQ as inverse of sum(abs(row of Q))
    for (int j = 0; j < opE->Nrows; j++) {
        for (int i = 0; i < opE->N; i++) {
            double temp = cvxmat_get(&(opE->E), j, i);
            opE->sigE.vals[j] += fabs(temp);
        }
        opE->sigE.vals[j] = 1.0 / opE->sigE.vals[j];
    }
}



/*
 * Reweight the constraint and update all the subsequent weightings, and also the current descent direction zE
 * basically weight_mod * E
 */
void cvxop_eddy_reweight(cvxop_eddy *opE, double weight_mod)
{
    double ww;
    for (int j = 0; j < opE->Nrows; j++) {
        ww = 1.0;
        if (opE->checks.vals[j] > 0) {
            ww = weight_mod;
            // ww = 1.0 * opE->checks.vals[j];
            // if (ww > weight_mod) {
            //     ww = weight_mod;
            // }
        }
        opE->weights.vals[j] *= ww;
        opE->zE.vals[j] *= ww;
    }

    // Make E as weighted E0
    for (int j = 0; j < opE->Nrows; j++) {
        for (int i = 0; i < opE->N; i++) {
            double temp = cvxmat_get(&(opE->E0), j, i);
            cvxmat_set(&(opE->E), j, i, opE->weights.vals[j] * temp);
        }
    }

    // Calculate sigE as inverse of sum(abs(row of Q))
    for (int j = 0; j < opE->Nrows; j++) {
        for (int i = 0; i < opE->N; i++) {
            double temp = cvxmat_get(&(opE->E), j, i);
            opE->sigE.vals[j] += fabs(temp);
        }
        opE->sigE.vals[j] = 1.0 / opE->sigE.vals[j];
    }
}


/**
 * Add absolute value of columns to the tau matrix 
 */
void cvxop_eddy_add2tau(cvxop_eddy *opE, cvx_mat *tau_mat)
{
    for (int j = 0; j < opE->Nrows; j++) {
        for (int i = 0; i < opE->N; i++) {
            double temp = cvxmat_get(&(opE->E), j, i);
            tau_mat->vals[i] += fabs(temp);
        }
    }    
}



/**
 * Step the gradient waveform (taumx)
 */
void cvxop_eddy_add2taumx(cvxop_eddy *opE, cvx_mat *taumx)
{   
    // MATH: taumx += E*zE
    for (int j = 0; j < opE->Nrows; j++) {
        for (int i = 0; i < opE->N; i++) {
            double temp = cvxmat_get(&(opE->E), j, i);
            taumx->vals[i] += (temp * opE->zE.vals[j]);
        }
    }

}


/**
 * Primal dual update
 */
void cvxop_eddy_update(cvxop_eddy *opE, cvx_mat *txmx, double rr)
{
    if (opE->Nrows > 0) {
        cvxmat_setvals(&(opE->Ex), 0.0);

        // MATH: Ex = E * txmx
        for (int j = 0; j < opE->Nrows; j++) {
            for (int i = 0; i < opE->N; i++) {
                double temp = cvxmat_get(&(opE->E), j, i) * txmx->vals[i];
                opE->Ex.vals[j] += temp;
            }
        }

        // MATH: Ex = Ex * sigE
        for (int j = 0; j < opE->Nrows; j++) {
            opE->Ex.vals[j] *= opE->sigE.vals[j];
        }

        // MATH: zEbuff  = zE + Ex = zE + sigE.*(E*txmx) 
        for (int j = 0; j < opE->Nrows; j++) {
            opE->zEbuff.vals[j] = opE->zE.vals[j] + opE->Ex.vals[j];
        }


        // MATH: zEbar = clip( zEbuff/sigE , [-upper_tol, lower_tol])
        double cushion = 0.99;
        for (int j = 0; j < opE->Nrows; j++) {
            double low =  (opE->goals.vals[j] - cushion*opE->tolerances.vals[j]) * opE->weights.vals[j];
            double high = (opE->goals.vals[j] + cushion*opE->tolerances.vals[j]) * opE->weights.vals[j];
            double val = opE->zEbuff.vals[j] / opE->sigE.vals[j];
            if (val < low) {
                opE->zEbar.vals[j] = low;
            } else if (val > high) {
                opE->zEbar.vals[j] = high;
            } else {
                opE->zEbar.vals[j] = val;
            }
            
        }


        // MATH: zEbar = zEbuff - sigE*zEbar
        for (int j = 0; j < opE->Nrows; j++) {
            opE->zEbar.vals[j] = opE->zEbuff.vals[j] - opE->sigE.vals[j] * opE->zEbar.vals[j];
        }

        // MATH: zE = rr*zEbar + (1-rr)*zE
        for (int j = 0; j < opE->Nrows; j++) {
            opE->zE.vals[j] = rr * opE->zEbar.vals[j] + (1 - rr) * opE->zE.vals[j];
        }
    }
}


/*
 * Check if moments are larger than a fixed tolerance
 */
int cvxop_eddy_check(cvxop_eddy *opE, cvx_mat *G)
{
    cvxmat_setvals(&(opE->Ex), 0.0);

    // MATH: Ex = E * txmx
    for (int j = 0; j < opE->Nrows; j++) {
        for (int i = 0; i < opE->N; i++) {
            double temp = cvxmat_get(&(opE->E0), j, i) * G->vals[i];
            opE->Ex.vals[j] += temp;
        }
    }

    // Set checks to be 0 if within tolerance, otherwise set to the ratio of eddy current to tolerance
    cvxmat_setvals(&(opE->checks), 0.0);
    for (int j = 0; j < opE->Nrows; j++) {
        if (fabs(opE->Ex.vals[j]) > opE->tolerances.vals[j]) {
            opE->checks.vals[j] = fabs(opE->Ex.vals[j]) / opE->tolerances.vals[j];
        }
    }

    // Set checks to be 0 if within tolerance, otherwise set to the ratio of eddy current to tolerance
    cvxmat_setvals(&(opE->checks), 0.0);
    for (int j = 0; j < opE->Nrows; j++) {
        double tol = opE->tolerances.vals[j];
        double low =  opE->goals.vals[j] - tol;
        double high = opE->goals.vals[j] + tol;
        double diff;
        if (opE->Ex.vals[j] < low) {
            diff = opE->goals.vals[j] - opE->Ex.vals[j];
            opE->checks.vals[j] = diff / tol;
        } else if (opE->Ex.vals[j] > high) {
            diff = opE->Ex.vals[j] - opE->goals.vals[j];
            opE->checks.vals[j] = diff / tol;
        }
    }

    int eddy_bad = 0;

    for (int j = 0; j < opE->Nrows; j++) {
        if (opE->checks.vals[j] > 0) {
             eddy_bad = 1;
        }
    }

    if (opE->verbose>0) {   
        printf("    Eddy check:     (%d)  [%.2e %.2e %.2e]  %.2e  %.2e  %.2e \n", eddy_bad, 
        opE->weights.vals[0], opE->weights.vals[1], opE->weights.vals[2], opE->Ex.vals[0], opE->Ex.vals[1], opE->Ex.vals[2]);

    }

    return eddy_bad;
}


/*
 * Free memory
 */
void cvxop_eddy_destroy(cvxop_eddy *opE)
{

    free(opE->norms.vals);
    free(opE->weights.vals);
    free(opE->checks.vals);
    free(opE->tolerances.vals);
    free(opE->goals.vals);


    free(opE->E0.vals);
    free(opE->E.vals);

    free(opE->sigE.vals);

    free(opE->zE.vals);
    free(opE->zEbuff.vals);
    free(opE->zEbar.vals);
    free(opE->Ex.vals);
}

