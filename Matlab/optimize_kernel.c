#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

//*******************
// Everything in the next section is used to time functions in Windows
// Annoyingly, in pure C this is difficult to do in a cross-platform manner
//*******************

// #include <time.h>
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
    #include <Windows.h>

    typedef struct {
        long tv_sec;
        long tv_usec;
    } timeval;

    typedef struct {
        int tz_minuteswest; 
        int tz_dsttime; 
    } timezone;

    int gettimeofday(struct timeval * tp, int *tzp)
    {
        // FILETIME Jan 1 1970 00:00:00
        // Note: some broken versions only have 8 trailing zero's, the correct epoch has 9 trailing zero's
        static const uint64_t EPOCH = ((uint64_t) 116444736000000000ULL); 

        SYSTEMTIME  nSystemTime;
        FILETIME    nFileTime;
        uint64_t    nTime;

        GetSystemTime( &nSystemTime );
        SystemTimeToFileTime( &nSystemTime, &nFileTime );
        nTime =  ((uint64_t)nFileTime.dwLowDateTime )      ;
        nTime += ((uint64_t)nFileTime.dwHighDateTime) << 32;

        tp->tv_sec  = (long) ((nTime - EPOCH) / 10000000L);
        tp->tv_usec = (long) (nSystemTime.wMilliseconds * 1000);
        return 0;
    }
#else
    #include <sys/time.h>
#endif

#include "cvx_matrix.h"
#include "op_slewrate.h"
#include "op_moments.h"
#include "op_eddy.h"
#include "op_beta.h"
#include "op_bval.h"
#include "op_gradient.h"
#include "op_pns.h"
#include "op_maxwell.h"
#include "te_finder.h"


#define EDDY_PARAMS_LEN 4
#define MOMENTS_PARAMS_LEN 7

void test_TE_finders();
void test_timer();

void cvx_optimize_kernel(cvx_mat *G, cvxop_gradient *opG, cvxop_slewrate *opD, 
                        cvxop_moments *opQ, cvxop_eddy *opE, cvxop_beta *opC, 
                        cvxop_bval *opB, cvxop_pns *opP, cvxop_maxwell *opX,
                         int N, double relax, int verbose, double bval_reduction, 
                         double *ddebug, int N_converge, double stop_increase,
                         int diffmode, double search_bval)
{    
    int max_iter = 30000;
    int check_amount = 100;
    int i_check = 0;
    int N_backlog = max_iter / check_amount;
    double *bval_backlog = (double *)malloc(N_backlog*sizeof(double));
    for (int i = 0; i < N_backlog; i++) {
        bval_backlog[i] = 0.0;
    } 
    
    int NDcheck = 1;
    
    int NDsum = 100;

    int NNi = 0;
    int NDsum_i = 0;
    int NDD = 4;

    int converge_count = 0;
    int limit_count = 0;

    int rebalance_count = 0;

    for (int i = 0; i < opB->zB.N; i++) {
        opB->zB.vals[i] = 0.0;
    } 
    for (int i = 0; i < opD->zD.N; i++) {
        opD->zD.vals[i] = 0.0;
    }
    for (int i = 0; i < opQ->zQ.N; i++) {
        opQ->zQ.vals[i] = 0.0;
    }

    cvx_mat xbar;
    copyNewMatrix(G, &xbar);
    cvxmat_setvals(&xbar, 0.0);

    cvx_mat taumx;
    copyNewMatrix(G, &taumx);
    cvxmat_setvals(&taumx, 0.0);

    cvx_mat txmx;
    copyNewMatrix(G, &txmx);
    cvxmat_setvals(&txmx, 0.0);

    // tau scaling vector
    cvx_mat tau;
    copyNewMatrix(G, &tau);
    cvxmat_setvals(&tau, 0.0);

    cvx_mat Ghist0;
    copyNewMatrix(G, &Ghist0);
    cvxmat_setvals(&Ghist0, 0.0);

    cvx_mat Ghist1;
    copyNewMatrix(G, &Ghist1);
    cvxmat_setvals(&Ghist1, 0.0);

    cvx_mat Gfullhist;
    cvxmat_alloc(&Gfullhist, NDsum*2*N, 1);
    cvxmat_setvals(&Gfullhist, 0.0);

    cvxop_slewrate_add2tau(opD, &tau);
    cvxop_moments_add2tau(opQ, &tau);
    cvxop_eddy_add2tau(opE, &tau);
    cvxop_bval_add2tau(opB, &tau);
    cvxop_beta_add2tau(opC, &tau);
    cvxop_pns_add2tau(opP, &tau);
    cvxop_maxwell_add2tau(opX, &tau);
    cvxmat_EWinvert(&tau);

    double obj0 = 1.0;
    double obj1 = 1.0;    

    int count = 0;
    while (count < max_iter) {
    
        // xbar = G-tau.*((D'*zD)+(Q'*zQ)+C'+(B'*zB))
        cvxmat_setvals(&taumx, 0.0);

        cvxop_slewrate_add2taumx(opD, &taumx);
        cvxop_moments_add2taumx(opQ, &taumx);
        cvxop_eddy_add2taumx(opE, &taumx);
        cvxop_beta_add2taumx(opC, &taumx);
        cvxop_bval_add2taumx(opB, &taumx);
        cvxop_pns_add2taumx(opP, &taumx);
        cvxop_maxwell_add2taumx(opX, &taumx);

        cvxmat_EWmultIP(&taumx, &tau);

        cvxmat_subractMat(&xbar, G, &taumx);

        // xbar = gradient_limits(xbar)
        cvxop_gradient_limiter(opG, &xbar);


        // txmx = 2*xbar-G;
        cvxmat_subractMatMult1(&txmx, 2.0, &xbar, G);
        for (int i = 0; i < N; i++) {
            opX->x_store.vals[i] = xbar.vals[i];
        }

  
        // zDbuff  = zD + sigD.*(D*txmx);
        // zQbuff  = zQ + sigQ.*(Q*txmx);

        // zDbar = zDbuff - sigD.*min(SRMAX,max(-SRMAX,zDbuff./sigD));
        // zQbar = zQbuff - sigQ.*min(mvec,max(-mvec,zQbuff./sigQ));

        // zD=p*zDbar+(1-p)*zD;
        // zQ=p*zQbar+(1-p)*zQ;
        cvxop_slewrate_update(opD, &txmx, relax);
        cvxop_moments_update(opQ, &txmx, relax);
        cvxop_eddy_update(opE, &txmx, relax);
        cvxop_bval_update(opB, &txmx, relax, ddebug, count);
        cvxop_pns_update(opP, &txmx, relax);
        cvxop_maxwell_update(opX, &txmx, relax);

        // G=p*xbar+(1-p)*G;
        cvxmat_updateG(G, relax, &xbar);

        // int ind_start = 100 + N * count;
        // if (ind_start < 990000) {
        //     for (int i = 0; i < N; i++) {
        //         ddebug[ind_start + i] = G->vals[i];
        //     }
        // }

        /*
        for (int i = (Gfullhist.N - G->N - 1); i >= 0; i--) {
            Gfullhist.vals[i + G->N] = Gfullhist.vals[i];
        }
        for (int i = 0; i < N; i++) {
            Gfullhist.vals[i] = G->vals[i];
            Ghist0.vals[i] = 0.0;
            Ghist1.vals[i] = 0.0;
        }
        for (int j = 0; j < NDsum; j++) {
            for (int i = 0; i < N; i++) {
                Ghist0.vals[i] += Gfullhist.vals[i + j * G->N];
                Ghist1.vals[i] += Gfullhist.vals[i + j * G->N + NDsum * G->N];
            }
        }
        for (int i = 0; i < N; i++) {
            Ghist0.vals[i] /= (double)NDsum;
            Ghist1.vals[i] /= (double)NDsum;
        }

        if ( count % NDcheck == 0 ) {
            double diff_norm = 0.0;
            for (int i = 0; i < N; i++) {
                diff_norm += pow((Ghist0.vals[i] - Ghist1.vals[i]), 2.0);
            }
            ddebug[100 + NDD*NNi] = diff_norm;
            NNi += 1;
        }
        */

        // Need checks here
        if ( count % check_amount == 0 ) {


            obj1 = cvxop_gradient_getbval(opG, G);
            
            bval_backlog[i_check] = obj1;
            
            double ii_backlog = 0;
            double bval_backlog0 = 0.0;
            for (int i = 0; i < 5; i++) {
                if ((i_check - i) >= 0) {
                    bval_backlog0 += bval_backlog[ (i_check - i) ];
                } else {
                    bval_backlog0 += 999999999999.0;
                }
                ii_backlog += 1.0;
            }
            bval_backlog0 /= ii_backlog;
            

            ii_backlog = 0;
            double bval_backlog1 = 0.0;
            for (int i = 5; i < 10; i++) {
                if ((i_check - i) >= 0) {
                    bval_backlog1 += bval_backlog[ (i_check - i) ];
                } else {
                    bval_backlog1 += 9999999999.0;
                }
                ii_backlog += 1.0;
            }
            bval_backlog1 /= ii_backlog;

            i_check += 1;

            double backlog_diff = (sqrt(bval_backlog0) - sqrt(bval_backlog1)) / sqrt(bval_backlog1);

            bval_backlog0 = (sqrt(obj1) - sqrt(bval_backlog1)) / sqrt(obj1);
            bval_backlog1 = (sqrt(obj1) - sqrt(bval_backlog1)) / sqrt(obj1);

            int is_converged = 0;
            if (fabs(backlog_diff) < stop_increase) {
                is_converged = 1;
            }

            rebalance_count += 1;


            if (verbose>0) {printf("\ncount = %d   rc = %d   obj = %.1f  backlog_diff = %.2e\n", count, rebalance_count, obj1, backlog_diff);}
            int bad_slew = cvxop_slewrate_check(opD, G);
            int bad_moments = cvxop_moments_check(opQ, G);
            int bad_gradient = cvxop_gradient_check(opG, G);
            int bad_eddy = cvxop_eddy_check(opE, G);
            int bad_pns = cvxop_pns_check(opP, G);
            int bad_maxwell = cvxop_maxwell_check(opX, G);

            int limit_break = 0;
            limit_break += bad_slew;
            limit_break += bad_moments;
            limit_break += bad_gradient;
            limit_break += bad_eddy;
            limit_break += bad_pns;

            if (diffmode == 0) {
                is_converged = 1;
                if ( (count > 0) && (limit_break == 0)) {
                    if (verbose > 0) {
                        printf("** Early termination at count = %d for free mode all constraints met\n", count);
                    }
                    break;
                }
            }

            if ( (search_bval > 0) && (obj1 > search_bval) && (limit_break == 0)) {
                if (verbose > 0) {
                    printf("** Early termination for search bval at count = %d   bval = %.1f > %.1f\n", count, obj1, search_bval);
                }
                break;
            }
            
            int finish_mod = 1;
            if (diffmode == 0) {
                finish_mod = 2;
            }

            if ( (count > 0) && (rebalance_count > finish_mod*N_converge)  && (is_converged > 0) && (limit_break == 0)) {
                if (verbose > 0) {
                    printf("** Early termination at count = %d   bval = %.1f\n", count, obj1);
                }
                break;
            }


           int needs_rebalancing = 0;
           if ( (diffmode > 0) && (is_converged > 0) && (rebalance_count > N_converge) ) {
               needs_rebalancing = 1;
           }

           if ( (diffmode == 0) && (rebalance_count > N_converge) && (limit_break > 0)) {
               needs_rebalancing = 1;
           }

        //    int gmax_overflow = cvxop_gradient_overflowcheck(opG, G);
        //    if (gmax_overflow > 0) {
        //        needs_rebalancing = 1;
        //    }

            if ( (bval_reduction > 0.0) && (count > 0) && (needs_rebalancing > 0) ) {

            
                rebalance_count = 0;

                if (verbose > 0) {
                    printf("\n\n !-!-!-!-!-!-! Converged to an inadequate waveform, reweighting !-!-!-!-!-!-! \n");
                } 

                // cvxop_slewrate_reweight(opD, 2.0);

                if (bad_moments > 0) {
                    cvxop_moments_reweight(opQ, 1.0 + bval_reduction);
                    // cvxop_slewrate_reweight(opD, 1.0 + bval_reduction/4.0);
                    if (verbose > 0) {printf("  ^^ moments ^^  ");}
                }
                if (bad_slew > 0) {
                    cvxop_slewrate_reweight(opD, 1.0 + bval_reduction);
                    if (verbose > 0) {printf("  ^^ slew ^^  ");}
                }
                if (bad_eddy > 0) {
                    cvxop_eddy_reweight(opE, 1.0 + 10.0*bval_reduction);
                    if (verbose > 0) {printf("  ^^ eddy ^^  ");}
                }
                if (bad_pns > 0) {
                    cvxop_pns_reweight(opP, 1.0 + 1.0*bval_reduction);
                    if (verbose > 0) {printf("  ^^ PNS ^^  ");}
                }     

                /*
                // if ((bad_slew < 1) && (bad_moments < 1) && (bad_eddy < 1) && (bad_pns < 1)) {
                    cvxop_bval_reweight(opB, 1.0 + bval_reduction/10.0);
                    // cvxop_beta_reweight(opC, 1.0 + bval_reduction/10.0);
                    if (verbose > 0) {printf("  ^^ bval ^^  ");}
                    for (int i = 0; i < opB->zB.N; i++) {
                        opB->zB.vals[i] = 0.0; 
                    }
                // }
                */
                           
                                   for (int i = 0; i < opB->zB.N; i++) {
                        opB->zB.vals[i] = 0.0; 
                    }
                if (verbose > 0) { printf("\n");}

                cvxmat_setvals(&tau, 0.0);
                cvxop_slewrate_add2tau(opD, &tau);
                cvxop_moments_add2tau(opQ, &tau);
                cvxop_eddy_add2tau(opE, &tau);
                cvxop_bval_add2tau(opB, &tau);
                cvxop_beta_add2tau(opC, &tau);
                cvxop_pns_add2tau(opP, &tau);
                cvxmat_EWinvert(&tau);
                
                
                for (int i = 0; i < opD->zD.N; i++) {
                    opD->zD.vals[i] = 0.0;
                }
                for (int i = 0; i < opQ->zQ.N; i++) {
                    opQ->zQ.vals[i] = 0.0;
                }
                for (int i = 0; i < opE->zE.N; i++) {
                    opE->zE.vals[i] = 0.0;
                }
                for (int i = 0; i < opP->zP.N; i++) {
                    opP->zP.vals[i] = 0.0;
                }
            }

            obj0 = obj1;
        }
    
        count++;
        ddebug[0] = count;
        fflush(stdout);
    }
    
    ddebug[13] = cvxop_gradient_getbval(opG, G);

    free(xbar.vals);
    free(taumx.vals);
    free(txmx.vals);
    free(tau.vals);
    free(bval_backlog);

    if (verbose > 0) {printf("  \n --- Final Checks --- \n");}
    int bad_slew = cvxop_slewrate_check(opD, G);
    int bad_moments = cvxop_moments_check(opQ, G);
    int bad_gradient = cvxop_gradient_check(opG, G);
    int bad_eddy = cvxop_eddy_check(opE, G);
    int bad_pns = cvxop_pns_check(opP, G);
    int bad_maxwell = cvxop_maxwell_check(opX, G);

    int limit_break = 0;
    limit_break += bad_slew;
    limit_break += bad_moments;
    limit_break += bad_gradient;
    limit_break += bad_eddy;
    limit_break += bad_pns;
    ddebug[14] = limit_break;

    ddebug[7] = bad_slew;
    ddebug[8] = bad_moments;
    ddebug[9] = bad_gradient;
}

void interp(cvx_mat *G, double dt_in, double dt_out, double TE, double T_readout) {
    int N1 = round((TE-T_readout) * 1.0e-3/dt_out);

    double *new_vals;
    double *temp_free = G->vals;
    new_vals = malloc(N1 * sizeof(double));


    double tt;
    double ti;
    int i0, i1;
    double d0, d1;
    double v0, v1;
    for (int i = 0; i < N1; i++) {
        ti = (dt_out * i) / dt_in;
        
        i0 = floor(ti);
        if (i0 < 0) {i0 = 0;} // This shouldn't happen unless some weird rounding and floor?
        
        i1 = i0+1;

        if (i1 < G->N) {
            d0 = fabs(ti-i1);
            d1 = 1.0 - d0;

            v0 = d0 * temp_free[i0];
            v1 = d1 * temp_free[i1];

            new_vals[i] = v0 + v1;
        } else {
            d0 = fabs(ti-i1);
            v0 = d0 * temp_free[i0]; 
            new_vals[i] = v0;
        }
    }

    G->vals = new_vals;
    G->N = N1;
    G->rows = N1;
    free(temp_free);
}



void run_kernel_diff(double **G_out, int *N_out, double **ddebug, int verbose,
                            int N, double dt, double gmax, double smax, double TE, 
                            int N_moments, double *moments_params, double PNS_thresh,  
                            double T_readout, double T_90, double T_180, int diffmode,
                            double bval_weight, double slew_weight, double moments_weight, 
                            double bval_reduce,  double dt_out,
                            int N_eddy, double *eddy_params,
                            int is_Gin, double *G_in, 
                            double search_bval,
                            int N_gfix, double *gfix,
                            double slew_reg,
                            int Naxis)
{
    double relax = 1.8;

    if (verbose > 0) {
        printf ("\nFirst pass, N = %d    dt = %.2e\n\n", N, dt);
        fflush(stdout);
    }

    /*
    // This is the old style, but I don't think its right for when inversion time doesn't line up with dt
    int ind_inv, ind_end90, ind_start180, ind_end180;
    if (diffmode > 0) {
        ind_inv = round((N + T_readout/(dt*1.0e3))/2.0);
        ind_end90 = floor(T_90*(1e-3/dt));
        ind_start180 = ind_inv - floor(T_180*(1e-3/dt/2));
        ind_end180 = ind_inv + floor(T_180*(1e-3/dt/2));
    } else {
        ind_inv = 9999999;
        ind_end90 = 0;
        ind_start180 = ind_inv;
        ind_end180 = ind_inv;
    }
    */
    
    
    // Calculate times of inversion and rf dead times
    double t_inv = (N*dt + 1e-3 * T_readout) / 2.0;
    double t_end90 = 1e-3 * T_90;
    double t_start180 = t_inv - T_180*1e-3/2.0;
    double t_stop180 = t_inv + T_180*1e-3/2.0;
    
    int ind_inv;
    int ind_end90;
    int ind_start180;
    int ind_end180;

    // Get indices from times, always rounding in the most conservative direction
    if (diffmode > 0) {
        ind_inv = round(t_inv/dt);
        ind_end90 = ceil(t_end90/dt);
        ind_start180 = floor(t_start180/dt);
        ind_end180 = ceil(t_stop180/dt);
    } else {
        ind_inv = 9999999;
        ind_end90 = 0;
        ind_start180 = ind_inv;
        ind_end180 = ind_inv;
    }
    

    if (verbose > 0) {
        printf ("\nN = %d  ind_inv = %d\n90_zeros = %d:%d    180_zeros = %d:%d\n\n", N, ind_inv, 0, ind_end90, ind_start180, ind_end180);
    }

    cvxop_beta opC;
    cvxop_bval opB;
    
    int N_converge = 1;
    double stop_increase = 1;

    if (diffmode == 1) {
        opB.active = 0; 
        opC.active = 1;
        N_converge = 24; 
        stop_increase = 1.0e-4;
    } else if (diffmode == 2) {
        opC.active = 0; 
        opB.active = 1;
        N_converge = 8;
        stop_increase = 1.0e-3; 
    } else {
        opC.active = 0; 
        opB.active = 0;
        N_converge = 8;
        stop_increase = 1.0e-3; 
    }

    cvxop_beta_init(&opC, N, dt, bval_weight, verbose);
    cvxop_bval_init(&opB, N, ind_inv, dt, bval_weight, verbose);

    
    cvxop_gradient opG;
    cvxop_gradient_init(&opG, N, Naxis, dt, gmax, ind_inv, verbose);
    cvxop_gradient_setFixRange(&opG, 0, ind_end90, 0.0);
    cvxop_gradient_setFixRange(&opG, ind_start180, ind_end180, 0.0);
    if (N_gfix > 0) {
        cvxop_gradient_setGfix(&opG, N_gfix, gfix);
    } 

    cvxop_slewrate opD;
    cvxop_slewrate_init(&opD, N, Naxis, dt, smax, slew_weight, slew_reg, verbose);

    cvxop_pns opP;
    cvxop_pns_init(&opP, N, Naxis, dt, ind_inv, PNS_thresh, 1.0, verbose);

    cvxop_maxwell opX;
    cvxop_maxwell_init(&opX, N, dt, ind_inv, .01, verbose);
    opX.active = 0;

    cvxop_moments opQ;
    cvxop_moments_init(&opQ, N, Naxis, ind_inv, dt, moments_weight, verbose);
    for (int i = 0; i < N_moments; i++) {
        cvxop_moments_addrow(&opQ, moments_params + (MOMENTS_PARAMS_LEN*i));
    }
    cvxop_moments_finishinit(&opQ);
    
    
    cvxop_eddy opE;
    cvxop_eddy_init(&opE, N, ind_inv, dt, .01, verbose);
    for (int i = 0; i < N_eddy; i++) {
        cvxop_eddy_addrow(&opE, (eddy_params[EDDY_PARAMS_LEN*i] * 1.0e-3), 
                                eddy_params[EDDY_PARAMS_LEN*i+1], 
                                eddy_params[EDDY_PARAMS_LEN*i+2], 
                                eddy_params[EDDY_PARAMS_LEN*i+3]);
    }
    cvxop_eddy_finishinit(&opE);
    
    cvx_mat G;
    
    cvxmat_alloc(&G, N*Naxis, 1);
    cvxmat_setvals(&G, 0.0);
    if (is_Gin == 0) {   
        cvxop_init_G(&opG, &G);
    } else {
        for (int i = 0; i < N*Naxis; i++) {
            G.vals[i] = G_in[i];
        }
    }

	*ddebug = (double *)malloc(100*sizeof(double));
    for (int i = 0; i < 100; i++) {
        (*ddebug)[i] = 0.0;
    }

    cvx_optimize_kernel(&G, &opG, &opD, &opQ, &opE, &opC, &opB, &opP, &opX, N, relax, verbose, bval_reduce, 
                        *ddebug, N_converge, stop_increase, diffmode, search_bval);

    cvxop_gradient_limiter(&opG, &G);
    
    if (verbose > 0) {
        printf ("\n****************************************\n");
        printf ("--- Finished diff kernel4 #1 in %d iterations  bvalue = %.2f", (int)(*ddebug)[0], (*ddebug)[13]);
        printf ("\n****************************************\n");
    }

    if (dt_out > 0) {
        interp(&G, dt, dt_out, TE, T_readout);
    }

    *N_out = G.rows;
    *G_out = G.vals;

    (*ddebug)[2] = opB.weight;
    (*ddebug)[3]  = opQ.weight;
    (*ddebug)[4]  = opD.weight;

    (*ddebug)[10] = opQ.norms.vals[0];
    (*ddebug)[11] = opQ.norms.vals[1];
    (*ddebug)[12] = opQ.norms.vals[2];

    cvxop_gradient_destroy(&opG);
    cvxop_slewrate_destroy(&opD);
    cvxop_moments_destroy(&opQ);
    cvxop_beta_destroy(&opC);
    cvxop_bval_destroy(&opB);

    fflush(stdout);

}


void run_kernel_diff_fixedN(double **G_out, int *N_out, double **ddebug, int verbose,
                            int N0, double gmax, double smax, double TE, 
                            int N_moments, double *moments_params, double PNS_thresh,  
                            double T_readout, double T_90, double T_180, int diffmode, double dt_out,
                            int N_eddy, double *eddy_params, double search_bval, double slew_reg)
{
    int N = N0;
    double dt = (TE-T_readout) * 1.0e-3 / (double) (N-1);

    // struct timeval start, end;
    // double diff;
    // gettimeofday(&start, NULL);

    run_kernel_diff(G_out, N_out, ddebug, verbose, 
                        N, dt, gmax, smax, TE, 
                        N_moments, moments_params, PNS_thresh,
                        T_readout, T_90, T_180, diffmode,
                        1.0, 1.0, 100.0, 
                        10.0,  dt_out,
                        N_eddy, eddy_params,
                        0, NULL, search_bval,
                        0, NULL, slew_reg,
                        1);
}


void run_kernel_diff_fixedN_Gin(double **G_out, int *N_out, double **ddebug, int verbose,
                                double *G_in, int N0, double gmax, double smax, double TE, 
                                int N_moments, double *moments_params, double PNS_thresh,  
                                double T_readout, double T_90, double T_180, int diffmode, double dt_out,
                                int N_eddy, double *eddy_params, double search_bval, double slew_reg)
{
    int N = N0;
    double dt = (TE-T_readout) * 1.0e-3 / (double) N;

    run_kernel_diff(G_out, N_out, ddebug, verbose,
                        N, dt, gmax, smax, TE, 
                        N_moments, moments_params, PNS_thresh,
                        T_readout, T_90, T_180, diffmode,
                        1.0, 1.0, 100.0, 
                        10.0,  dt_out,
                        N_eddy, eddy_params,
                        1, G_in, search_bval,
                        0, NULL, slew_reg, 1);
}




void run_kernel_diff_fixeddt(double **G_out, int *N_out, double **ddebug, int verbose,
                            double dt0, double gmax, double smax, double TE,
                            int N_moments, double *moments_params, double PNS_thresh, 
                            double T_readout, double T_90, double T_180, int diffmode, double dt_out,
                            int N_eddy, double *eddy_params, double search_bval, double slew_reg, int Naxis)
{
    int N = round((TE-T_readout) * 1.0e-3/dt0);
    if (N < 5) {
        printf ("\nWARNING: N = %d looks too small, setting to 5\n\n", N);
        N = 5;
    }

    // double dt = (TE-T_readout) * 1.0e-3 / (double) N;
    double dt = dt0;

    // struct timeval start, end;
    // double diff;
    // gettimeofday(&start, NULL);

    run_kernel_diff(G_out, N_out, ddebug, verbose,
                            N, dt, gmax, smax, TE, 
                            N_moments, moments_params, PNS_thresh,
                            T_readout, T_90, T_180, diffmode,
                            10.0, 1.0, 10.0, 
                            10.0,  dt_out,
                            N_eddy, eddy_params,
                            0, NULL, search_bval,
                            0, NULL, slew_reg, Naxis);

    // gettimeofday(&end, NULL);
    // diff = (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
    // if (verbose > 0) {
    //     printf("\nOperation took %f ms (%f ms total)\n", (1.0e3*diff), 1.0e3*diff);
    // }
    // (*ddebug)[15] = diff;
}


void run_kernel_diff_fixeddt_fixG(double **G_out, int *N_out, double **ddebug, int verbose,
                            double dt0, double gmax, double smax, double TE,
                            int N_moments, double *moments_params, double PNS_thresh, 
                            double T_readout, double T_90, double T_180, int diffmode, double dt_out,
                            int N_eddy, double *eddy_params, double search_bval,
                            int N_gfix, double *gfix, double slew_reg)
{
    int N = round((TE-T_readout) * 1.0e-3/dt0);
    if (N < 5) {
        printf ("\nWARNING: N = %d looks too small, setting to 5\n\n", N);
        N = 5;
    }

    // double dt = (TE-T_readout) * 1.0e-3 / (double) N;
    double dt = dt0;

    run_kernel_diff(G_out, N_out, ddebug, verbose,
                            N, dt, gmax, smax, TE, 
                            N_moments, moments_params, PNS_thresh,
                            T_readout, T_90, T_180, diffmode,
                            1.0, 1.0, 10.0, 
                            10.0,  dt_out,
                            N_eddy, eddy_params,
                            0, NULL, search_bval,
                            N_gfix, gfix, slew_reg, 1);

}


int main (void)
{
    printf ("In optimize_kernel.c main function\n");
    
    // test_timer();
    test_TE_finders();
    return 0;
}


// This example shows the fast TE finders, which can scan different TE values quickly by
// performing quick initial searches and using multiple threads to seach multiple TE at once
void test_TE_finders()
{
    // 1 = betamax
    // 2 = bval max
    int diffmode = 2;

    double *G;
    int N;
    double *debug;

    int N_eddy = 0; 
    double *eddy_params;

    int N_moments = 3; 
    double *moment_params = (double *)malloc(N_moments*7*sizeof(double));
    for (int i = 0; i < N_moments; i++) {
        int ii = i*7;
        moment_params[ii+0] = 0; moment_params[ii+1] = i; moment_params[ii+2] = 0; 
        moment_params[ii+3] = -1; moment_params[ii+4] = -1; 
        moment_params[ii+5] = 0; moment_params[ii+6] = 1.0e-3;
    }

    double PNS_thresh = -1.0;

    double dt = 400e-6;
    double gmax = .05;
    double smax = 50.0;

    double T_readout = 12.0;
    double T_90 = 4.0;
    double T_180 = 8.0;

    // double dt_out = 10.0e-6;
    double dt_out = -1.0;

    double bval = 300;

    struct timeval start, end;
    double diff;


    // openMP TE finder example
    gettimeofday(&start, NULL);
    minTE_diff_par(&G, &N, &debug, 0, dt, gmax, smax, bval, 
                        N_moments, moment_params, PNS_thresh, 
                            T_readout, T_90, T_180, diffmode, dt_out, N_eddy, eddy_params, 1.0);

    gettimeofday(&end, NULL);
    diff = (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
    printf("OpenMP operation took %.2f ms (TE = %.2f ms)\n", (1.0e3*diff), 1.0e3*N*dt);


    // Single thread TE finder example
    gettimeofday(&start, NULL);
    minTE_diff(&G, &N, &debug, 0, dt, gmax, smax, bval, 
                        N_moments, moment_params, PNS_thresh, 
                            T_readout, T_90, T_180, diffmode, dt_out, N_eddy, eddy_params, 1.0);

    gettimeofday(&end, NULL);
    diff = (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
    printf("Non-threaded operation took %.2f ms (TE = %.2f ms)\n", (1.0e3*diff), 1.0e3*N*dt);

    return;
}

void test_timer()
{
    // 1 = betamax
    // 2 = bval max
    int diffmode = 2;

    double *G;
    int N;
    double *debug;

    // int N_eddy = 1; 
    // double *eddy_params;
    
    int N_eddy = 1; 
    double *eddy_params = (double *)malloc(4*sizeof(double));
    eddy_params[0] = 50.0;
    eddy_params[1] = 0.0;
    eddy_params[2] = 1.0e-4;
    eddy_params[3] = 0.0;


    int N_moments = 3; 
    double *moment_params = (double *)malloc(21*sizeof(double));
    int ii = 0;
    moment_params[ii+0] = 0; moment_params[ii+1] = 0; moment_params[ii+2] = 0; 
    moment_params[ii+3] = -1; moment_params[ii+4] = -1; 
    moment_params[ii+5] = 0; moment_params[ii+6] = 1.0e-3;
    
    ii = 7;
    moment_params[ii+0] = 0; moment_params[ii+1] = 1; moment_params[ii+2] = 0; 
    moment_params[ii+3] = -1; moment_params[ii+4] = -1; 
    moment_params[ii+5] = 0; moment_params[ii+6] = 1.0e-3;

    ii = 14;
    moment_params[ii+0] = 0; moment_params[ii+1] = 2; moment_params[ii+2] = 0; 
    moment_params[ii+3] = -1; moment_params[ii+4] = -1; 
    moment_params[ii+5] = 0; moment_params[ii+6] = 1.0e-3;

    double PNS_thresh = 1.0;

    double dt = 100e-6;
    double gmax = .05;
    double smax = 50.0;
    double TE = 80.0;

    double T_readout = 12.0;
    double T_90 = 4.0;
    double T_180 = 8.0;

    // double dt_out = 10.0e-6;
    double dt_out = -1.0;

    int N_time = 10;
    // struct timeval start, end;
    // double diff;
    // gettimeofday(&start, NULL);

    for (int i = 0; i < N_time; i++) {
        run_kernel_diff_fixeddt(&G, &N, &debug, 1, dt, gmax, smax, TE, 
                            N_moments, moment_params, PNS_thresh, 
                            T_readout, T_90, T_180, diffmode, dt_out, N_eddy, eddy_params, 100.0, 1.0, 1);
    }

    // gettimeofday(&end, NULL);
    // diff = (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
    // printf("\nOperation took %f ms (%f ms total)\n", (1.0e3*diff/N_time), 1.0e3*diff);

    return;
}