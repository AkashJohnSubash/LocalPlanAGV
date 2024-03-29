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

#ifndef ACADOS_SOLVER_tricycle_frenet_H_
#define ACADOS_SOLVER_tricycle_frenet_H_

#include "acados/utils/types.h"

#include "acados_c/ocp_nlp_interface.h"
#include "acados_c/external_function_interface.h"

#define TRICYCLE_FRENET_NX     5
#define TRICYCLE_FRENET_NZ     0
#define TRICYCLE_FRENET_NU     2
#define TRICYCLE_FRENET_NP     16
#define TRICYCLE_FRENET_NBX    0
#define TRICYCLE_FRENET_NBX0   5
#define TRICYCLE_FRENET_NBU    2
#define TRICYCLE_FRENET_NSBX   0
#define TRICYCLE_FRENET_NSBU   0
#define TRICYCLE_FRENET_NSH    7
#define TRICYCLE_FRENET_NSH0   0
#define TRICYCLE_FRENET_NSG    0
#define TRICYCLE_FRENET_NSPHI  0
#define TRICYCLE_FRENET_NSHN   7
#define TRICYCLE_FRENET_NSGN   0
#define TRICYCLE_FRENET_NSPHIN 0
#define TRICYCLE_FRENET_NSPHI0 0
#define TRICYCLE_FRENET_NSBXN  0
#define TRICYCLE_FRENET_NS     7
#define TRICYCLE_FRENET_NS0    0
#define TRICYCLE_FRENET_NSN    7
#define TRICYCLE_FRENET_NG     0
#define TRICYCLE_FRENET_NBXN   0
#define TRICYCLE_FRENET_NGN    0
#define TRICYCLE_FRENET_NY0    7
#define TRICYCLE_FRENET_NY     7
#define TRICYCLE_FRENET_NYN    5
#define TRICYCLE_FRENET_N      100
#define TRICYCLE_FRENET_NH     7
#define TRICYCLE_FRENET_NHN    7
#define TRICYCLE_FRENET_NH0    0
#define TRICYCLE_FRENET_NPHI0  0
#define TRICYCLE_FRENET_NPHI   0
#define TRICYCLE_FRENET_NPHIN  0
#define TRICYCLE_FRENET_NR     0

#ifdef __cplusplus
extern "C" {
#endif


// ** capsule for solver data **
typedef struct tricycle_frenet_solver_capsule
{
    // acados objects
    ocp_nlp_in *nlp_in;
    ocp_nlp_out *nlp_out;
    ocp_nlp_out *sens_out;
    ocp_nlp_solver *nlp_solver;
    void *nlp_opts;
    ocp_nlp_plan_t *nlp_solver_plan;
    ocp_nlp_config *nlp_config;
    ocp_nlp_dims *nlp_dims;

    // number of expected runtime parameters
    unsigned int nlp_np;

    /* external functions */
    // dynamics

    external_function_param_casadi *forw_vde_casadi;
    external_function_param_casadi *expl_ode_fun;




    // cost

    external_function_param_casadi *cost_y_fun;
    external_function_param_casadi *cost_y_fun_jac_ut_xt;
    external_function_param_casadi *cost_y_hess;



    external_function_param_casadi cost_y_0_fun;
    external_function_param_casadi cost_y_0_fun_jac_ut_xt;
    external_function_param_casadi cost_y_0_hess;



    external_function_param_casadi cost_y_e_fun;
    external_function_param_casadi cost_y_e_fun_jac_ut_xt;
    external_function_param_casadi cost_y_e_hess;


    // constraints
    external_function_param_casadi *nl_constr_h_fun_jac;
    external_function_param_casadi *nl_constr_h_fun;






    external_function_param_casadi nl_constr_h_e_fun_jac;
    external_function_param_casadi nl_constr_h_e_fun;

} tricycle_frenet_solver_capsule;

ACADOS_SYMBOL_EXPORT tricycle_frenet_solver_capsule * tricycle_frenet_acados_create_capsule(void);
ACADOS_SYMBOL_EXPORT int tricycle_frenet_acados_free_capsule(tricycle_frenet_solver_capsule *capsule);

ACADOS_SYMBOL_EXPORT int tricycle_frenet_acados_create(tricycle_frenet_solver_capsule * capsule);

ACADOS_SYMBOL_EXPORT int tricycle_frenet_acados_reset(tricycle_frenet_solver_capsule* capsule, int reset_qp_solver_mem);

/**
 * Generic version of tricycle_frenet_acados_create which allows to use a different number of shooting intervals than
 * the number used for code generation. If new_time_steps=NULL and n_time_steps matches the number used for code
 * generation, the time-steps from code generation is used.
 */
ACADOS_SYMBOL_EXPORT int tricycle_frenet_acados_create_with_discretization(tricycle_frenet_solver_capsule * capsule, int n_time_steps, double* new_time_steps);
/**
 * Update the time step vector. Number N must be identical to the currently set number of shooting nodes in the
 * nlp_solver_plan. Returns 0 if no error occurred and a otherwise a value other than 0.
 */
ACADOS_SYMBOL_EXPORT int tricycle_frenet_acados_update_time_steps(tricycle_frenet_solver_capsule * capsule, int N, double* new_time_steps);
/**
 * This function is used for updating an already initialized solver with a different number of qp_cond_N.
 */
ACADOS_SYMBOL_EXPORT int tricycle_frenet_acados_update_qp_solver_cond_N(tricycle_frenet_solver_capsule * capsule, int qp_solver_cond_N);
ACADOS_SYMBOL_EXPORT int tricycle_frenet_acados_update_params(tricycle_frenet_solver_capsule * capsule, int stage, double *value, int np);
ACADOS_SYMBOL_EXPORT int tricycle_frenet_acados_update_params_sparse(tricycle_frenet_solver_capsule * capsule, int stage, int *idx, double *p, int n_update);

ACADOS_SYMBOL_EXPORT int tricycle_frenet_acados_solve(tricycle_frenet_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT int tricycle_frenet_acados_free(tricycle_frenet_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT void tricycle_frenet_acados_print_stats(tricycle_frenet_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT int tricycle_frenet_acados_custom_update(tricycle_frenet_solver_capsule* capsule, double* data, int data_len);


ACADOS_SYMBOL_EXPORT ocp_nlp_in *tricycle_frenet_acados_get_nlp_in(tricycle_frenet_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT ocp_nlp_out *tricycle_frenet_acados_get_nlp_out(tricycle_frenet_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT ocp_nlp_out *tricycle_frenet_acados_get_sens_out(tricycle_frenet_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT ocp_nlp_solver *tricycle_frenet_acados_get_nlp_solver(tricycle_frenet_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT ocp_nlp_config *tricycle_frenet_acados_get_nlp_config(tricycle_frenet_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT void *tricycle_frenet_acados_get_nlp_opts(tricycle_frenet_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT ocp_nlp_dims *tricycle_frenet_acados_get_nlp_dims(tricycle_frenet_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT ocp_nlp_plan_t *tricycle_frenet_acados_get_nlp_plan(tricycle_frenet_solver_capsule * capsule);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_SOLVER_tricycle_frenet_H_
