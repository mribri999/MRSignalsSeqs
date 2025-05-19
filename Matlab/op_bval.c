#include "op_bval.h"

/**
 * Initialize the opB struct
 * This is the operator that maximizes b-value for diffusion purposes
 */
void cvxop_bval_init(cvxop_bval *opB, int N, int ind_inv, double dt, double init_weight, int verbose) {
    
    opB->N = N;
    opB->dt = dt;
    opB->ind_inv = ind_inv;
    opB->verbose = verbose;
    opB->weight = init_weight;
    
    cvxmat_alloc(&opB->sigBdenom, N, 1);
    cvxmat_alloc(&opB->sigB, N, 1);

    cvxmat_alloc(&opB->csx, N, 1);
    cvxmat_alloc(&opB->Btau, N, 1);

    cvxmat_alloc(&opB->zB, N, 1);
    cvxmat_alloc(&opB->zBbuff, N, 1);
    cvxmat_alloc(&opB->Bx, N, 1);

    if (opB->active > 0) {

            // These are all derived from the matrix B, which isn't calculated anymore
            // TODO:  Add some matlab/python code to show this matrix

            double mat_norm = 0.0;
            for (int i = 0; i < N; i++) {
                mat_norm += (2.0*i + 1.0) * pow((double)(N-i), 2.0);
            }
            mat_norm = sqrt(mat_norm);

            opB->mat_norm = mat_norm;
            
            opB->weight /= mat_norm;

            for (int i = 0; i < opB->N; i++) {
                opB->sigBdenom.vals[opB->N-1-i] = opB->weight * ((double)i + 1.0) * ((double)opB->N - ((double)i/2.0));
            } 

            for (int i = 0; i < opB->N; i++) {
                opB->sigB.vals[i] = 1.0/opB->sigBdenom.vals[i];
            }

    }
}

/**
 * Reweight the constraint and update all the subsequent weightings, and also the current descent direction zB
 * basically weight_mod * B
 */
void cvxop_bval_reweight(cvxop_bval *opB, double weight_mod)
{
    opB->weight *= weight_mod;

    for (int i = 0; i < opB->N; i++) {
        opB->sigBdenom.vals[opB->N-1-i] = opB->weight * ((double)i + 1.0) * ((double)opB->N - ((double)i/2.0));
    } 

    for (int i = 0; i < opB->N; i++) {
        opB->sigB.vals[i] = 1.0/opB->sigBdenom.vals[i];
    }

    for (int i = 0; i < opB->zB.N; i++) {
        opB->zB.vals[i] *= weight_mod;
    }
    
}

/**
 * Add sigBdenom to the tau matrix 
 */
void cvxop_bval_add2tau(cvxop_bval *opB, cvx_mat *tau_mat)
{
    if (opB->active > 0) {
        for (int i = 0; i < opB->N; i++) {
            tau_mat->vals[i] += fabs(opB->sigBdenom.vals[i]);
        }
    }
}


/**
 * This replaces the old matrix multiplication, O(2N) instead of O(N^2)
 * out = B*in
 */
void cvxmat_bval_multB(cvx_mat *out, cvxop_bval *opB, cvx_mat *in) {
    
    // Integration
    double gt = 0.0;
    for (int i = 0; i < opB->N; i++) {
        if (i < opB->ind_inv) {
            gt += in->vals[i];
        } else {
            gt -= in->vals[i];
        }
        opB->csx.vals[i] = gt;
    }  

    // Integration transpose
    gt = 0.0;
    for (int i = (opB->N-1); i >= 0; i--) {
        gt += opB->csx.vals[i];
        if (i < opB->ind_inv) {
            out->vals[i] = gt;
        } else {
            out->vals[i] = -gt;
        }
    }

    // Apply weight
    for (int i = 0; i < out->N; i++) {
        out->vals[i] *= opB->weight;
    }
}

/**
 * Step the gradient waveform (taumx)
 */
void cvxop_bval_add2taumx(cvxop_bval *opB, cvx_mat *taumx)
{
    if (opB->active > 0) {
        
        // MATH: Btau = B*zB
        cvxmat_setvals(&(opB->Btau), 0.0);
        cvxmat_bval_multB(&opB->Btau, opB, &opB->zB);

        // MATH: taumx -= Btau
        for (int i = 1; i < taumx->N; i++) {
            taumx->vals[i] -= opB->Btau.vals[i];
        }

    }
}

/**
 * Primal dual update
 */
void cvxop_bval_update(cvxop_bval *opB, cvx_mat *txmx, double rr)
{
    if (opB->active > 0) {

        // MATH: Bx = (B*txmx)
        cvxmat_setvals(&(opB->Bx), 0.0);
        cvxmat_bval_multB(&opB->Bx, opB, txmx);

        // MATH: Bx = sigB*Bx
        for (int i = 0; i < opB->Bx.N; i++) {
            opB->Bx.vals[i] *= opB->sigB.vals[i];
        }

        // MATH: zBbuff  = = zB + Bx = zB + sigB.*(B*txmx) 
        for (int i = 0; i < opB->zBbuff.N; i++) {
            opB->zBbuff.vals[i] = opB->zB.vals[i] + opB->Bx.vals[i];
        }

        // rr is the overrelaxation parameter
        // MATH: zB = rr*zBbuff + (1-rr)*zB
        for (int i = 0; i < opB->zB.N; i++) {
            opB->zB.vals[i] = rr * opB->zBbuff.vals[i] + (1 - rr) * opB->zB.vals[i];
        }
    }
}


void cvxop_bval_destroy(cvxop_bval *opB)
{
    free(opB->sigBdenom.vals);
    free(opB->sigB.vals);

    free(opB->csx.vals);
    free(opB->Btau.vals);
    
    free(opB->zB.vals);
    free(opB->zBbuff.vals);
    free(opB->Bx.vals);
}
