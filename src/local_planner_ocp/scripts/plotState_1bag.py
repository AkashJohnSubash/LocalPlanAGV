
from bagpy import bagreader
import numpy as np
import pandas as pd
import matplotlib as mpl
from common import *
# from matplotlib import pyplot as plt
from sys import argv
import shutil

from plot_mpl import plotOptVars, plotCosts, plotResiduals
from sys_dynamics import SysDyn

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


b1 = bagreader(file1)

csvfiles = []
SCALE_BOUND = 1.05
# get all from MPC node
vars1_csv   = b1.message_by_topic('/mpc_variables')
metric1_csv   = b1.message_by_topic('/mpc_metrics')
cmd1_csv   = b1.message_by_topic('/cmd_vel_auto')
param1_csv   = b1.message_by_topic('/mpc_parameters')

offset_k = 0
mpc1_vars = pd.read_csv(vars1_csv)
mpc1_metric = pd.read_csv(metric1_csv)
mpc1_param = pd.read_csv(param1_csv)
cmd1_vel   = pd.read_csv(cmd1_csv)

iter1 = np.shape(mpc1_metric['sim_time'])[0] - offset_k
t01 = np.reshape(mpc1_vars['sim_time'][:], [1, iter1])

# iter1 = np.shape(mpc1_metric['sim_time'][:]) - offset_k
# t01 = np.reshape(mpc1_vars['sim_time'][:], [1, iter1])

lift1 = mpc1_param['lifted'][0]

'''Plot state and control'''

if(lift1 == False):
    s01 = np.reshape(mpc1_vars['zeta0_0'][offset_k:], (1, iter1))
    n01 = np.reshape(mpc1_vars['zeta0_1'][offset_k:], (1, iter1))
    beta01 = np.reshape(mpc1_vars['zeta0_2'][offset_k:], (1, iter1))
    
    v01 = np.reshape(mpc1_vars['zeta0_3'][offset_k:], (1, iter1))
    alpha01 = np.reshape(mpc1_vars['zeta0_4'][offset_k:], (1, iter1))
    
    # zeta0 = np.vstack((s0, n0, beta0, v0, alpha0))

    # Project frenet state to cartesian
    zetaC1_hat = np.zeros((3, iter1))
    x_i, y_i, phi_i = sysModel.Fren2CartT(s01, n01, beta01)
    zetaC1_hat[0, : ] = np.ravel(x_i)
    zetaC1_hat[1, : ] = np.ravel(y_i)
    zetaC1_hat[2, : ] = np.ravel(phi_i)
            
    x01 = np.ravel(zetaC1_hat[0, :-1])
    y01 = np.ravel(zetaC1_hat[1 , :-1])
    phi01 = np.ravel(zetaC1_hat[2, :-1])

else:
    x01 = np.reshape(mpc1_vars['zeta0_0'][offset_k:], (1, iter1))
    y01 = np.reshape(mpc1_vars['zeta0_1'][offset_k:], (1, iter1))
    phi01 = np.reshape(mpc1_vars['zeta0_2'][offset_k:], (1, iter1))
    
    s01 = np.reshape(mpc1_vars['zeta0_3'][offset_k:], (1, iter1))
    n01 = np.reshape(mpc1_vars['zeta0_4'][offset_k:], (1, iter1))
    beta01 = np.reshape(mpc1_vars['zeta0_5'][offset_k:], (1, iter1))
    
    v01 = np.reshape(mpc1_vars['zeta0_6'][offset_k:], (1, iter1))
    alpha01 = np.reshape(mpc1_vars['zeta0_7'][offset_k:], (1, iter1))



a01 = np.reshape(mpc1_vars['u0_0'][offset_k:], (1, iter1))
omg01 = np.reshape(mpc1_vars['u0_1'][offset_k:], (1, iter1))

# print("DEBUG", np.shape(slack1_lower))
x1 = np.ravel(x01)
y1 = np.ravel(y01)
phi1 = np.ravel(phi01)


s1 = np.ravel(s01)
n1 = np.ravel(n01)
beta1 = np.ravel(beta01)

v1 = np.ravel(v01)
alpha1 = np.ravel(alpha01)

a1 = np.ravel(a01)
omg1 = np.ravel(omg01)

t1 = np.ravel(t01)

fig1 = mpl.pyplot.figure(figsize=(6.4, 4.8), num='zeta')
zetaCAx = fig1.add_subplot(2, 1, 1)
zetaCAx.stairs(x1[:], t1[:], baseline=None,label="$x\,(\mathrm{m})$", color="lightcoral" )
zetaCAx.stairs(y1[:], t1[:], baseline=None,label="$y\,(\mathrm{m})$", color="teal")
zetaCAx.stairs(phi1[:], t1[:], baseline=None,label="$\\varphi\,(\mathrm{rad})$", color="#220d7eb6", alpha=0.55)
zetaCAx.set_ylim( np.hstack([x1, y1]).max(axis=0) * SCALE_BOUND, 
                  np.hstack([x1, y1]).min(axis=0) * SCALE_BOUND)
zetaCAx.set_xlim(0, t1.max(axis=0))

zetaCAx.set_ylabel("$\\zeta^{c}$")
zetaCAx.invert_yaxis() 
zetaCAx.legend(loc='upper right')
zetaCAx.grid()
# zetafAx.set_yscale('symlog')

zetaFAx = fig1.add_subplot(2, 1, 2)
zetaFAx.stairs(s1[:-1], t1[:], baseline=None,label="$s\,(\mathrm{m})$", color="lightcoral" )
zetaFAx.stairs(n1[:-1], t1[:], baseline=None,label="$n\,(\mathrm{m})$", color="teal" )
zetaFAx.stairs(beta1[:-1], t1[:], baseline=None,label="$\\beta\,(\mathrm{rad})$", color="#220d7eb6", alpha=0.55)
zetaFAx.set_ylim( np.hstack([s1, n1, beta1]).max(axis=0)  + 5, 
                  np.hstack([s1, n1, beta1]).min(axis=0)  * SCALE_BOUND)
zetaFAx.set_xlim(0, t1.max(axis=0))
zetaFAx.set_yscale('symlog')
zetaFAx.set_ylabel("$\\zeta^{f}$")
zetaFAx.set_xlabel('time (s)')
zetaFAx.legend(loc='upper right')
zetaFAx.grid()
zetaFAx.invert_yaxis()

# fig1_1 = mpl.pyplot.figure(figsize=(6.4, 2.65), num='zeta_1')
fig2 = mpl.pyplot.figure(figsize=(6.4, 4.8), num='u')
zetaUAx = fig2.add_subplot(2, 1, 1)
zetaUAx.stairs(v1[:-1], t1[:], baseline=None,label="$v\,(\mathrm{m\,s^{-1}})$", color="lightcoral" )
zetaUAx.stairs(alpha1[:-1], t1[:], baseline=None,label="$\\alpha\,(\mathrm{rad})$", color="teal")
zetaUAx.set_ylim( np.hstack([v1, alpha1]).max(axis=0)  + 0.1, 
                  np.hstack([v1, alpha1]).min(axis=0)  * SCALE_BOUND)
zetaUAx.set_xlim(0, t1.max(axis=0))

# zetaUAx.set_xlabel('time (s)')
zetaUAx.set_ylabel("$\\zeta^{u}$")
zetaUAx.legend(loc='upper right')
zetaUAx.grid()
zetaUAx.invert_yaxis()

fig1.tight_layout()
fig1.savefig('zeta_time.pdf', format='pdf', bbox_inches='tight')

# fig1_1.tight_layout()
# fig1_1.savefig('zeta_time1.pdf', format='pdf', bbox_inches='tight')

zetaUAx = fig2.add_subplot(2, 1, 1)
zetaUAx.stairs(v1[:-1], t1[:], baseline=None,label="$v\,(\mathrm{m\,s^{-1}})$", color="lightcoral" )
zetaUAx.stairs(alpha1[:-1], t1[:], baseline=None,label="$\\alpha\,(\mathrm{rad})$", color="teal")
zetaUAx.set_ylim( np.hstack([v1, alpha1]).max(axis=0)  + 0.1, 
                  np.hstack([v1, alpha1]).min(axis=0)  * SCALE_BOUND)
zetaUAx.set_xlim(0, t1.max(axis=0))

# zetaUAx.set_xlabel('time (s)')
zetaUAx.set_ylabel("$\\zeta^{u}$")
zetaUAx.legend(loc='upper right')
zetaUAx.grid()
zetaUAx.invert_yaxis()


UAx = fig2.add_subplot(2, 1, 2)
UAx.stairs(a1[:-1], t1[:], baseline=None,label="$a\,(\mathrm{m\,s^{-2}})$", color="lightcoral" )
UAx.stairs(omg1[:-1], t1[:], baseline=None,label="$\\omega\,(\mathrm{rad})$", color="teal")
UAx.set_ylim( np.hstack([a1, omg1]).max(axis=0) * SCALE_BOUND, 
             np.hstack([a1, omg1]).min(axis=0) * SCALE_BOUND)
UAx.set_xlim(0, t1.max(axis=0))

UAx.set_ylabel("$u$")
UAx.invert_yaxis() 
UAx.legend(loc='upper right')
UAx.grid()

UAx.set_xlabel('time (s)')
UAx.set_ylabel("$u$")
UAx.legend(loc='upper right')

UAx.invert_yaxis()

fig2.tight_layout()
fig2.savefig('u_time.pdf', format='pdf', bbox_inches='tight')

# fig2.tight_layout()
# fig2.savefig('zeta_u.pdf', format='pdf', bbox_inches='tight')

# fig1.tight_layout()
# fig1.savefig('zeta_f.pdf', format='pdf', bbox_inches='tight')
mpl.pyplot.show()