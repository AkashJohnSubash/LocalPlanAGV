
from bagpy import bagreader
import numpy as np
import pandas as pd
import matplotlib as mpl
from common import *
# from matplotlib import pyplot as plt
from sys import argv
import shutil

from plot_mpl import plotOptVars, plotCosts, plotResiduals

text_usetex = True if shutil.which('latex') else False
params = {
        'text.latex.preamble': r"\usepackage{gensymb} \usepackage{amsmath}",
        'axes.labelsize': 12,
        'axes.titlesize': 12,
        'legend.fontsize': 12,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'text.usetex': text_usetex,
        'font.family': 'serif'
}

SCALE_BOUND = 1.05
mpl.rcParams.update(params)

file1 = argv[1]
file2 = argv[2]

b1 = bagreader(file1)
b2 = bagreader(file2)

csvfiles = []

# get all from MPC node
vars1_csv   = b1.message_by_topic('/mpc_variables')
metric1_csv   = b1.message_by_topic('/mpc_metrics')
cmd1_csv   = b1.message_by_topic('/cmd_vel_auto')
param1_csv   = b1.message_by_topic('/mpc_parameters')

vars2_csv   = b2.message_by_topic('/mpc_variables')
metric2_csv   = b2.message_by_topic('/mpc_metrics')
cmd2_csv   = b2.message_by_topic('/cmd_vel_auto')
param2_csv   = b2.message_by_topic('/mpc_parameters')

offset_k = 30
mpc1_vars = pd.read_csv(vars1_csv)
mpc1_metric = pd.read_csv(metric1_csv)
mpc1_param = pd.read_csv(param1_csv)
cmd1_vel   = pd.read_csv(cmd1_csv)

mpc2_vars = pd.read_csv(vars2_csv)
mpc2_metric = pd.read_csv(metric2_csv)
mpc2_param = pd.read_csv(param2_csv)
cmd2_vel   = pd.read_csv(cmd2_csv)

iter1 = np.shape(mpc1_metric['cost'])[0] - offset_k
t01 = np.ravel(mpc1_vars['sim_time'])[:]#, [iter, 1])

iter2 = np.shape(mpc2_metric['cost'])[0] -offset_k
t02 = np.ravel(mpc2_vars['sim_time'])[:]#, [iter, 1])

'''Plot state and control'''

lift1 = mpc1_param['lifted'][0]
lift2 = mpc2_param['lifted'][0]

# s01 = np.reshape(mpc1_vars['zeta0_0'][offset_k:], (1, iter1))
# s02 = np.reshape(mpc2_vars['zeta0_0'][offset_k:], (1, iter2))

cost1 = np.reshape(mpc1_metric['cost'][offset_k:], (1, iter1))
stat1 = np.reshape(mpc1_metric['residuals_0'][offset_k+1:], (1, iter1-1))
ineq1 = np.reshape(mpc1_metric['residuals_2'][offset_k+1:], (1, iter1-1))

cost2 = np.reshape(mpc2_metric['cost'][offset_k:], (1, iter2))
stat2 = np.reshape(mpc2_metric['residuals_0'][offset_k+1:], (1, iter2-1))
ineq2 = np.reshape(mpc2_metric['residuals_2'][offset_k+1:], (1, iter2-1))

ns1 = 17#int(mpc1_param['ns'][0])
slack1_lower = np.reshape(mpc1_metric[f'slacks_{ns1}'][offset_k:], (1, iter1))
for i in range(ns1 +1 , 2 * ns1):
    sl1_0 = np.reshape(mpc1_metric[f'slacks_{i}'][offset_k:], (1, iter1))
    slack1_lower = np.vstack((slack1_lower, sl1_0))

slack2_lower = np.reshape(mpc2_metric[f'slacks_{ns1}'][offset_k:], (1, iter2))
for i in range(ns1 +1 , 2 * ns1):
    sl2_0 = np.reshape(mpc2_metric[f'slacks_{i}'][offset_k:], (1, iter2))
    slack2_lower = np.vstack((slack2_lower, sl2_0))


if(lift1 == False):
    s01 = np.reshape(mpc1_vars['zeta0_0'][offset_k:], (1, iter1))

else: 
    s01 = np.reshape(mpc1_vars['zeta0_3'][offset_k:], (1, iter1))
    
if(lift2 == False):
    s02 = np.reshape(mpc2_vars['zeta0_0'][offset_k:], (1, iter2))
    
else: 
    s02 = np.reshape(mpc2_vars['zeta0_3'][offset_k:], (1, iter2))

print("DEBUG", np.shape(slack1_lower))
s1 = np.ravel(s01)
s2 = np.ravel(s02)

cost1 = np.ravel(cost1)
cost2 = np.ravel(cost2)

stat1 = np.ravel(stat1)
stat2 = np.ravel(stat2)

ineq1 = np.ravel(ineq1)
ineq2 = np.ravel(ineq2)

# print("\n oscillation 2nd derv:", v1_osc, v2_osc)
fig1 = mpl.pyplot.figure(num='Costs')
costAx = fig1.add_subplot(2, 1, 1)
costAx.stairs(cost1[:-1], s1[:], baseline=None,label="DEF", color="lightcoral" )
costAx.stairs(cost2[:-1], s2[:], baseline=None,label="LF", color="teal")
costAx.set_ylim( np.hstack([cost1,cost2]).max(axis=0) *SCALE_BOUND, 
                 np.hstack([cost1,cost2]).min(axis=0) - 2.5)
costAx.set_xlim(0, np.amax(np.hstack([s1,s2]).max(axis=0)))

# v_agv.set_xlabel('$track\,progress\,s$ (m)')
costAx.set_ylabel("Costs")
costAx.invert_yaxis() 
costAx.legend(loc='upper right')
costAx.grid()
# costAx.set_yscale('symlog')

statAx = fig1.add_subplot(2, 1, 2)
statAx.stairs(stat1[:], s1[:], baseline=None,label="DEF", color="lightcoral" )
statAx.stairs(stat2[:], s2[:], baseline=None,label="LF", color="teal")
statAx.set_ylim( np.hstack([stat1, stat2]).max(axis=0) *SCALE_BOUND , 
                 np.hstack([stat1, stat2]).min(axis=0) *SCALE_BOUND - 0.1)
statAx.set_xlim(0, np.amax(np.hstack([s1, s2]).max(axis=0)))

statAx.set_xlabel('track progress $s\,(\\mathrm{m})$')
statAx.set_ylabel("Stationarity residual")
statAx.legend(loc='upper right')
statAx.grid()
statAx.invert_yaxis()
statAx.set_yscale('symlog')



fig2 = mpl.pyplot.figure(num='constraints')
obstBreachAx = fig2.add_subplot(2, 1, 1)
obstBreachAx.stairs(slack1_lower[0, :-1], s1[:], baseline=None,label="DEF", color="lightcoral")
for i in range(1, obst_constr_len):
    obstBreachAx.stairs(slack1_lower[i, :-1], s1[:], baseline=None, color="lightcoral" )

obstBreachAx.stairs(slack2_lower[0, :-1], s2[:], baseline=None,label="LF", color="teal")
for i in range(1, obst_constr_len):
    obstBreachAx.stairs(slack2_lower[i, :-1], s2[:], baseline=None, color="teal" )

obstBreachAx.set_ylim( np.hstack([np.ravel(slack1_lower[:, :]),
                                  np.ravel(slack2_lower[:, :])]).max(axis=0) * SCALE_BOUND, 
                       np.hstack([np.ravel(slack1_lower[0, :]),
                                  np.ravel(slack2_lower[0, :])]).min(axis=0) * SCALE_BOUND -0.0025 )

obstBreachAx.set_xlim(0, np.amax(np.hstack([s1, s2]).max(axis=0)))
obstBreachAx.set_ylabel("Obstacle breach")
obstBreachAx.set_xlabel('track progress $s\,(\\mathrm{m})$')
obstBreachAx.legend(loc='upper right')
obstBreachAx.grid()
obstBreachAx.invert_yaxis()

ineqAx = fig2.add_subplot(2, 1, 2)
ineqAx.stairs(ineq1[:], s1[:], baseline=None,label="DEF", color="lightcoral" )
ineqAx.stairs(ineq2[:], s2[:], baseline=None,label="LF", color="teal")
ineqAx.set_ylim( np.hstack([ineq1, ineq2]).max(axis=0) * SCALE_BOUND, 
                 np.hstack([ineq1, ineq2]).min(axis=0) * SCALE_BOUND -0.025)
ineqAx.set_xlim(0, np.amax(np.hstack([s1, s2]).max(axis=0)))

ineqAx.set_xlabel('track progress $s\,(\\mathrm{m})$')
ineqAx.set_ylabel("Inequality residual")
ineqAx.legend(loc='upper right')
ineqAx.grid()
ineqAx.invert_yaxis()

fig2.tight_layout()
fig2.savefig('comp_cnstr.pdf', format='pdf', bbox_inches='tight')
fig1.tight_layout()
fig1.savefig('comp_costs.pdf', format='pdf', bbox_inches='tight')
mpl.pyplot.show()