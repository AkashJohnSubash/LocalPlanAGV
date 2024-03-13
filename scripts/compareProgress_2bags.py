
from bagpy import bagreader
import numpy as np
import pandas as pd
import matplotlib as mpl
from local_planner_ocp.common import *
# from matplotlib import pyplot as plt
from sys import argv
import shutil

from local_planner_ocp.plot_mpl import plotOptVars, plotCosts, plotResiduals
from local_planner_ocp.sys_dynamics import SysDyn

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

sysModel = SysDyn()
print("DEBUG", text_usetex)
mpl.rcParams.update(params)

file1 = argv[1]
file2 = argv[2]

b1 = bagreader(file1)
b2 = bagreader(file2)

csvfiles = []
SCALE_BOUND = 1.005
# get all from MPC node
vars1_csv   = b1.message_by_topic('/mpc_variables')
metric1_csv   = b1.message_by_topic('/mpc_metrics')
cmd1_csv   = b1.message_by_topic('/cmd_vel_auto')
param1_csv   = b1.message_by_topic('/mpc_parameters')

vars2_csv   = b2.message_by_topic('/mpc_variables')
metric2_csv   = b2.message_by_topic('/mpc_metrics')
cmd2_csv   = b2.message_by_topic('/cmd_vel_auto')
param2_csv   = b2.message_by_topic('/mpc_parameters')

offset_k = 0
mpc1_vars = pd.read_csv(vars1_csv)
mpc1_metric = pd.read_csv(metric1_csv)
mpc1_param = pd.read_csv(param1_csv)
cmd1_vel   = pd.read_csv(cmd1_csv)

mpc2_vars = pd.read_csv(vars2_csv)[:-1]
mpc2_metric = pd.read_csv(metric2_csv)
mpc2_param = pd.read_csv(param2_csv)
cmd2_vel   = pd.read_csv(cmd2_csv)

iter1 = np.shape(mpc1_metric['soln_time'])[0] - offset_k
iter2 = np.shape(mpc2_metric['soln_time'])[0] - offset_k

print(iter1)
lift1 = mpc1_param['lifted'][0]
lift2 = mpc2_param['lifted'][0]
'''Plot state and control'''

if(lift1 == False):
    s01 = np.reshape(mpc1_vars['zeta0_0'][offset_k:], (1, iter1))

else:
    s01 = np.reshape(mpc1_vars['zeta0_3'][offset_k:], (1, iter1-1))


if(lift2 == False):
    s02 = np.reshape(mpc2_vars['zeta0_0'][offset_k:], (1, iter2))

else:
    s02 = np.reshape(mpc2_vars['zeta0_3'][offset_k:], (1, iter2-1))

t01 = np.reshape(mpc1_vars['sim_time'][offset_k:], (1, iter1))
t02 = np.reshape(mpc2_vars['sim_time'][offset_k:], (1, iter2-1))

s1 = np.ravel(s01)
s2 = np.ravel(s02)

t1 = np.ravel(t01)
t2 = np.ravel(t02)

cost1 = np.reshape(mpc1_metric['cost'][offset_k:], (1, iter1))
cost2 = np.reshape(mpc2_metric['cost'][offset_k:], (1, iter2))
cost1 = np.ravel(cost1)
cost2 = np.ravel(cost2)



fig1 = mpl.pyplot.figure(num='progress')
zetaFAx = fig1.add_subplot(1, 1, 1)
zetaFAx.stairs(s1[:-1], t1[:],baseline=None,label="Direct", color="lightcoral" )
zetaFAx.stairs(s2[:-1], t2[:], baseline=None,label="Lifted", color="teal")
zetaFAx.set_xlim( np.hstack([t1, t2]).max(axis=0) * SCALE_BOUND, 
                 np.hstack([t1, t2]).min(axis=0) * SCALE_BOUND)
zetaFAx.set_ylim(0, np.hstack([s1,s2]).max(axis=0))

zetaFAx.set_ylabel("Track progress $\mathrm{s}$ (m)")
zetaFAx.invert_xaxis() 
zetaFAx.legend(loc='upper right')
zetaFAx.grid()
zetaFAx.set_xlabel('Simulation time (s)')

# # print("\n oscillation 2nd derv:", v1_osc, v2_osc)
# # fig1 = mpl.pyplot.figure(num='Costs')
# costAx = fig1.add_subplot(2, 1, 2)
# costAx.stairs(cost1[:-1], s1[:], baseline=None,label="Direct", color="lightcoral" )
# costAx.stairs(cost2[:-2], s2[:], baseline=None,label="Lifted", color="teal")
# costAx.set_ylim( np.hstack([cost1, cost2]).max(axis=0) * 1.05, 
#                  np.hstack([cost1, cost2]).min(axis=0) )
# costAx.set_xlim(0, np.amax(np.hstack([s1,s2]).max(axis=0)))

# # v_agv.set_xlabel('$track\,progress\,s$ (m)')
# costAx.set_ylabel("Costs")
# costAx.invert_yaxis() 
# costAx.legend(loc='upper right')
# costAx.grid()

fig1.tight_layout()
fig1.savefig('soln_times.pdf', format='pdf', bbox_inches='tight')
mpl.pyplot.show()