/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
 */


// standard
#include <stdio.h>
#include <stdlib.h>
// acados
#include "acados/utils/print.h"
#include "acados/utils/math.h"
#include "acados_c/ocp_nlp_interface.h"
#include "acados_c/external_function_interface.h"
#include "acados_solver_tricycle_frenet.h"

// blasfeo
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"

#define NX     TRICYCLE_FRENET_NX
#define NZ     TRICYCLE_FRENET_NZ
#define NU     TRICYCLE_FRENET_NU
#define NP     TRICYCLE_FRENET_NP
#define NBX    TRICYCLE_FRENET_NBX
#define NBX0   TRICYCLE_FRENET_NBX0
#define NBU    TRICYCLE_FRENET_NBU
#define NSBX   TRICYCLE_FRENET_NSBX
#define NSBU   TRICYCLE_FRENET_NSBU
#define NSH    TRICYCLE_FRENET_NSH
#define NSG    TRICYCLE_FRENET_NSG
#define NSPHI  TRICYCLE_FRENET_NSPHI
#define NSHN   TRICYCLE_FRENET_NSHN
#define NSGN   TRICYCLE_FRENET_NSGN
#define NSPHIN TRICYCLE_FRENET_NSPHIN
#define NSBXN  TRICYCLE_FRENET_NSBXN
#define NS     TRICYCLE_FRENET_NS
#define NSN    TRICYCLE_FRENET_NSN
#define NG     TRICYCLE_FRENET_NG
#define NBXN   TRICYCLE_FRENET_NBXN
#define NGN    TRICYCLE_FRENET_NGN
#define NY0    TRICYCLE_FRENET_NY0
#define NY     TRICYCLE_FRENET_NY
#define NYN    TRICYCLE_FRENET_NYN
#define NH     TRICYCLE_FRENET_NH
#define NPHI   TRICYCLE_FRENET_NPHI
#define NHN    TRICYCLE_FRENET_NHN
#define NH0    TRICYCLE_FRENET_NH0
#define NPHIN  TRICYCLE_FRENET_NPHIN
#define NR     TRICYCLE_FRENET_NR


int main()
{

    tricycle_frenet_solver_capsule *acados_ocp_capsule = tricycle_frenet_acados_create_capsule();
    // there is an opportunity to change the number of shooting intervals in C without new code generation
    int N = TRICYCLE_FRENET_N;
    // allocate the array and fill it accordingly
    double* new_time_steps = NULL;
    int status = tricycle_frenet_acados_create_with_discretization(acados_ocp_capsule, N, new_time_steps);

    if (status)
    {
        printf("tricycle_frenet_acados_create() returned status %d. Exiting.\n", status);
        exit(1);
    }

    ocp_nlp_config *nlp_config = tricycle_frenet_acados_get_nlp_config(acados_ocp_capsule);
    ocp_nlp_dims *nlp_dims = tricycle_frenet_acados_get_nlp_dims(acados_ocp_capsule);
    ocp_nlp_in *nlp_in = tricycle_frenet_acados_get_nlp_in(acados_ocp_capsule);
    ocp_nlp_out *nlp_out = tricycle_frenet_acados_get_nlp_out(acados_ocp_capsule);
    ocp_nlp_solver *nlp_solver = tricycle_frenet_acados_get_nlp_solver(acados_ocp_capsule);
    void *nlp_opts = tricycle_frenet_acados_get_nlp_opts(acados_ocp_capsule);

    // initial condition
    int idxbx0[NBX0];
    idxbx0[0] = 0;
    idxbx0[1] = 1;
    idxbx0[2] = 2;
    idxbx0[3] = 3;
    idxbx0[4] = 4;

    double lbx0[NBX0];
    double ubx0[NBX0];
    lbx0[0] = 0;
    ubx0[0] = 0;
    lbx0[1] = 0;
    ubx0[1] = 0;
    lbx0[2] = 0;
    ubx0[2] = 0;
    lbx0[3] = 0;
    ubx0[3] = 0;
    lbx0[4] = 0;
    ubx0[4] = 0;

    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "idxbx", idxbx0);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "lbx", lbx0);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "ubx", ubx0);

    // initialization for state values
    double x_init[NX];
    x_init[0] = 0.0;
    x_init[1] = 0.0;
    x_init[2] = 0.0;
    x_init[3] = 0.0;
    x_init[4] = 0.0;

    // initial value for control input
    double u0[NU];
    u0[0] = 0.0;
    u0[1] = 0.0;
    // set parameters
    double p[NP];
    p[0] = -12.5;
    p[1] = -0.75;
    p[2] = 1;
    p[3] = 1.570796326794895;
    p[4] = -8;
    p[5] = 10;
    p[6] = 0.7853981633974475;
    p[7] = 0;
    p[8] = 20;
    p[9] = 20;
    p[10] = 1;
    p[11] = 0;
    p[12] = 20;
    p[13] = 20;
    p[14] = 0.5;
    p[15] = 0;

    for (int ii = 0; ii <= N; ii++)
    {
        tricycle_frenet_acados_update_params(acados_ocp_capsule, ii, p, NP);
    }
  

    // prepare evaluation
    int NTIMINGS = 1;
    double min_time = 1e12;
    double kkt_norm_inf;
    double elapsed_time;
    int sqp_iter;

    double xtraj[NX * (N+1)];
    double utraj[NU * N];


    // solve ocp in loop
    int rti_phase = 0;

    for (int ii = 0; ii < NTIMINGS; ii++)
    {
        // initialize solution
        for (int i = 0; i < N; i++)
        {
            ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "x", x_init);
            ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "u", u0);
        }
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, N, "x", x_init);
        ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "rti_phase", &rti_phase);
        status = tricycle_frenet_acados_solve(acados_ocp_capsule);
        ocp_nlp_get(nlp_config, nlp_solver, "time_tot", &elapsed_time);
        min_time = MIN(elapsed_time, min_time);
    }

    /* print solution and statistics */
    for (int ii = 0; ii <= nlp_dims->N; ii++)
        ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, ii, "x", &xtraj[ii*NX]);
    for (int ii = 0; ii < nlp_dims->N; ii++)
        ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, ii, "u", &utraj[ii*NU]);

    printf("\n--- xtraj ---\n");
    d_print_exp_tran_mat( NX, N+1, xtraj, NX);
    printf("\n--- utraj ---\n");
    d_print_exp_tran_mat( NU, N, utraj, NU );
    // ocp_nlp_out_print(nlp_solver->dims, nlp_out);

    printf("\nsolved ocp %d times, solution printed above\n\n", NTIMINGS);

    if (status == ACADOS_SUCCESS)
    {
        printf("tricycle_frenet_acados_solve(): SUCCESS!\n");
    }
    else
    {
        printf("tricycle_frenet_acados_solve() failed with status %d.\n", status);
    }

    // get solution
    ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, 0, "kkt_norm_inf", &kkt_norm_inf);
    ocp_nlp_get(nlp_config, nlp_solver, "sqp_iter", &sqp_iter);

    tricycle_frenet_acados_print_stats(acados_ocp_capsule);

    printf("\nSolver info:\n");
    printf(" SQP iterations %2d\n minimum time for %d solve %f [ms]\n KKT %e\n",
           sqp_iter, NTIMINGS, min_time*1000, kkt_norm_inf);

    // free solver
    status = tricycle_frenet_acados_free(acados_ocp_capsule);
    if (status) {
        printf("tricycle_frenet_acados_free() returned status %d. \n", status);
    }
    // free solver capsule
    status = tricycle_frenet_acados_free_capsule(acados_ocp_capsule);
    if (status) {
        printf("tricycle_frenet_acados_free_capsule() returned status %d. \n", status);
    }

    return status;
}
