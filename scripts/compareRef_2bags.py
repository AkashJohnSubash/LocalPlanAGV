
from bagpy import bagreader
import numpy as np
import pandas as pd
import matplotlib as mpl
# from matplotlib import pyplot as plt
from sys import argv
import shutil

from local_planner_ocp.plot_mpl import plotOptVars, plotCosts, plotResiduals

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

print("DEBUG", text_usetex)
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


mpc1_vars = pd.read_csv(vars1_csv)
mpc1_metric = pd.read_csv(metric1_csv)
cmd1_vel   = pd.read_csv(cmd1_csv)
mpc1_param = pd.read_csv(param1_csv)

mpc2_vars = pd.read_csv(vars2_csv)
mpc2_metric = pd.read_csv(metric2_csv)
cmd2_vel   = pd.read_csv(cmd2_csv)
mpc2_param = pd.read_csv(param2_csv)

iter1 = np.shape(mpc1_vars['zeta0_0'])[0]
t01 = np.ravel(mpc1_vars['sim_time'])[:]#, [iter, 1])

iter2 = np.shape(mpc2_vars['zeta0_0'])[0]
t02 = np.ravel(mpc2_vars['sim_time'])[:]#, [iter, 1])

'''Plot state and control'''


lift1 = mpc1_param['lifted'][0]
lift2 = mpc2_param['lifted'][0]
'''Plot state and control'''

if(lift1 == False):
    s01 = np.reshape(mpc1_vars['zeta0_0'][:], (1, iter1))

else:
    s01 = np.reshape(mpc1_vars['zeta0_3'][:], (1, iter1))


if(lift2 == False):
    s02 = np.reshape(mpc2_vars['zeta0_0'][:], (1, iter2))


else:
    s02 = np.reshape(mpc2_vars['zeta0_3'][:], (1, iter2))

# v1_agv = np.reshape(cmd1_vel['twist.linear.x'], (1, iter1))
# omg1_agv = np.reshape(cmd1_vel['twist.angular.z'], (1, iter1))

# v2_agv = np.reshape(cmd2_vel['twist.linear.x'], (1, iter2))
# omg2_agv = np.reshape(cmd2_vel['twist.angular.z'], (1, iter2))

v1 = np.reshape(mpc1_vars['zeta0_3'][:], (1, iter1))
alpha1 = np.reshape(mpc1_vars['zeta0_4'][:], (1, iter1))
a1 = np.reshape(mpc1_vars['u0_0'][:], (1, iter1))
omg1 = np.reshape(mpc1_vars['u0_1'][:], (1, iter1))

v2 = np.reshape(mpc2_vars['zeta0_3'][:], (1, iter2))
alpha2 = np.reshape(mpc2_vars['zeta0_4'][:], (1, iter2))
a2 = np.reshape(mpc2_vars['u0_0'][:], (1, iter2))
omg2 = np.reshape(mpc2_vars['u0_1'][:], (1, iter2))


n1_sample = np.shape(v1)[1]
sum = 0
for i in range(1, n1_sample-1):
    sum = sum + np.square(v1[0, i+1] - 2*v1[0, i] + v1[0, i-1])
v1_osc = np.sqrt(sum/n1_sample)

sum = 0
for i in range(1, n1_sample-1):
    sum = sum + np.square(a1[0, i+1] - 2*a1[0, i] + a1[0, i-1])
a1_osc = np.sqrt(sum/n1_sample)

sum = 0
for i in range(1, n1_sample-1):
    sum = sum + np.square(a1[0, i+1] - 2*a1[0, i] + a1[0, i-1])
a1_osc = np.sqrt(sum/n1_sample)

sum = 0
for i in range(1, n1_sample-1):
    sum = sum + np.square(omg1[0, i+1] - 2*omg1[0, i] + omg1[0, i-1])
omg1_osc = np.sqrt(sum/n1_sample)

n2_sample = np.shape(v2)[1]
sum = 0
for i in range(1, n2_sample-1):
    sum = sum + np.square(v2[0, i+1] - 2*v2[0, i] + v2[0, i-1])
v2_osc = np.sqrt(sum/n2_sample)

sum = 0
for i in range(1, n2_sample-1):
    sum = sum + np.square(a2[0, i+1] - 2*a2[0, i] + a2[0, i-1])
a2_osc = np.sqrt(sum/n2_sample)

sum = 0
for i in range(1, n2_sample-1):
    sum = sum + np.square(omg2[0, i+1] - 2*omg2[0, i] + omg2[0, i-1])
omg2_osc = np.sqrt(sum/n1_sample)


s1 = np.ravel(s01)
v1 = np.ravel(v1)
alpha1 = np.ravel(alpha1)
s2 = np.ravel(s02)
v2 = np.ravel(v2)
alpha2 = np.ravel(alpha2)


print("\n oscillation v 2nd derv:", v1_osc, v2_osc)
print("\n oscillation a 2nd derv:", a1_osc, a2_osc)
print("\n oscillation omega 2nd derv:", omg1_osc, omg2_osc)
fig1 = mpl.pyplot.figure(num='Trajectory reference')
vAx = fig1.add_subplot(2, 1, 1)
vAx.stairs(v1[:-1], s1[:], baseline=None,label="Ref1", linewidth=1.25, color="teal" )
vAx.stairs(v2[:-1], s2[:], baseline=None,label="Ref2", color="lightcoral")
vAx.set_ylim( 0.82,#np.hstack(0.85[v1,v2]).max(axis=0) + 0.01, 
                0.52)#np.hstack([v1,v2]).min(axis=0) -0.01 )
vAx.set_xlim(0, np.amax(np.hstack([s1[:-65],s2[:-65]]).max(axis=0)))

# v_agv.set_xlabel('$track\,progress\,s$ (m)')
vAx.set_ylabel("$v\,(\\mathrm{m\,s^{-1}})$")
vAx.invert_yaxis() 
vAx.legend(loc='upper right')
vAx.grid()

alphaAx = fig1.add_subplot(2, 1, 2)
alphaAx.stairs(alpha1[:-1], s1[:], baseline=None,label="Ref1", color="teal" )
alphaAx.stairs(alpha2[:-1], s2[:], baseline=None,label="Ref2", color="lightcoral")
alphaAx.set_ylim( np.hstack([alpha1,alpha2]).max(axis=0) *1.05, 
                np.hstack([alpha1,alpha2]).min(axis=0) *1.05 )
alphaAx.set_xlim(0, np.amax(np.hstack([s1[:-65],s2[:-65]]).max(axis=0)))

alphaAx.set_xlabel('track progress $s\,(\\mathrm{m})$')
alphaAx.set_ylabel("$\\alpha\,(\\mathrm{rad})$")
alphaAx.legend(loc='upper right')
alphaAx.grid()
alphaAx.invert_yaxis()
fig1.tight_layout()

fig1.savefig('ucomp_ref.pdf', format='pdf', bbox_inches='tight')
# # plt.plot(y)
# mpl.pyplot.savefig("test1.svg", format="svg")
mpl.pyplot.show()