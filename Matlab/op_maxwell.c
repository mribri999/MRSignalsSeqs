#include "op_maxwell.h"

/**
 * Initialize the opX struct
 */
void cvxop_maxwell_init(cvxop_maxwell *opX, int N, double dt, int ind_inv, double tol, int verbose) 
{
    opX->active = 1;
    opX->N = N;
    opX->dt = dt;
    opX->ind_inv = ind_inv; 
    opX->verbose = verbose;
    opX->tol = tol;
   
    opX->sign = 1.0;
    opX->sigX = 1.0;

    cvxmat_alloc(&opX->x_store, N, 1);
}

/**
 * Add absolute volume of columns to tau
 */
void cvxop_maxwell_add2tau(cvxop_maxwell *opX, cvx_mat *tau_mat)
{
    if (opX->active > 0) {

        int N = tau_mat->rows;
        for (int i = 0; i < N; i++) {
            tau_mat->vals[i] += 1.0;
        }

    }
}


void cvxop_maxwell_add2taumx(cvxop_maxwell *opX, cvx_mat *taumx)
{   
    if (opX->active > 0) {
        // MATH: taumx += E*zE
        for (int i = 0; i < opX->N; i++) {
            double temp = opX->x_store.vals[i];
            if (i > opX->ind_inv) {temp *= -1.0;}
            if (opX->sign < 0) {temp *= -1.0;}
            taumx->vals[i] -= 1.0*(temp * opX->zX);
        }
    }
}


/**
 * Primal dual update
 */
void cvxop_maxwell_update(cvxop_maxwell *opX, cvx_mat *txmx, double rr)
{
    if (opX->active > 0) {
        
        double x0 = 0.0;
        double x1 = 0.0;
        for (int i = 0; i < (opX->N); i++) {
            if (i < opX->ind_inv) {
                x0 += txmx->vals[i] * txmx->vals[i];
            } else {
                x1 += txmx->vals[i] * txmx->vals[i];
            }
        }
    
        opX->Xx = sqrt(x0) - sqrt(x1);
        if (opX->Xx < 0) {
            opX->sign = -1.0;
        } else {
            opX->sign = 1.0;
        }
        opX->Xx = fabs(opX->Xx);
        opX->zXbuff = opX->zX + opX->Xx;


        // MATH: zEbar = clip( zEbuff/sigE , [-upper_tol, lower_tol])
        opX->zXbar = (opX->zXbuff - opX->tol) / opX->zXbuff;
        if (opX->zXbar < 0.0) {
            opX->zXbar = 0.0;
        }

        // MATH: zEbar = zEbuff - sigE*zEbar
        opX->zXbar = opX->zXbuff - opX->sigX * opX->zXbar;

        
        // MATH: zE = rr*zEbar + (1-rr)*zE
        rr = 1.0;
        opX->zX = rr * opX->zXbar + (1 - rr) * opX->zX;

        // for (int i = 0; i < (opX->N); i++) {
        //     opX->x_store.vals[i] = txmx->vals[i];
        // }

    }
}



int cvxop_maxwell_check(cvxop_maxwell *opX, cvx_mat *G)
{

    double x0 = 0.0;
    double x1 = 0.0;
    for (int i = 0; i < (opX->N); i++) {
        if (i < opX->ind_inv) {
            x0 += G->vals[i] * G->vals[i];
        } else {
            x1 += G->vals[i] * G->vals[i];
        }
    }
    
    double x_diff = sqrt(x0) - sqrt(x1);

    if (opX->verbose>0) {  
        printf("    Maxwell:        (%d)  [%.2e]  %.2f \n", 0, 1.0, x_diff);
    }

    return 0;
}