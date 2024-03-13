
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
file2 = argv[2]

b1 = bagreader(file1)
b2 = bagreader(file2)

csvfiles = []
SCALE_BOUND = 1.05
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

mpc2_vars = pd.read_csv(vars2_csv)
mpc2_metric = pd.read_csv(metric2_csv)
mpc2_param = pd.read_csv(param2_csv)
cmd2_vel   = pd.read_csv(cmd2_csv)

iter1 = np.shape(mpc1_vars['zeta0_0'])[0] - offset_k
t01 = np.ravel(mpc1_vars['sim_time'])[:]#, [iter, 1])

iter2 = np.shape(mpc2_vars['zeta0_0'])[0] -offset_k
t02 = np.ravel(mpc2_vars['sim_time'])[:]#, [iter, 1])

lift1 = mpc1_param['lifted'][0]
lift2 = mpc2_param['lifted'][0]
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

if(lift2 == False):
    s02 = np.reshape(mpc2_vars['zeta0_0'][offset_k:], (1, iter2))
    n02 = np.reshape(mpc2_vars['zeta0_1'][offset_k:], (1, iter2))
    beta02 = np.reshape(mpc2_vars['zeta0_2'][offset_k:], (1, iter2))
    
    v02 = np.reshape(mpc2_vars['zeta0_3'][offset_k:], (1, iter2))
    alpha02 = np.reshape(mpc2_vars['zeta0_4'][offset_k:], (1, iter2))
    
    # Project frenet state to cartesian
    zetaC2_hat = np.zeros((3, iter2))
    x_i, y_i, phi_i = sysModel.Fren2CartT(s01, (n01), beta01)
    zetaC2_hat[0, : ] = np.ravel(x_i)
    zetaC2_hat[1, : ] = np.ravel(y_i)
    zetaC2_hat[2, : ] = np.ravel(phi_i)
            
    x02 = np.ravel(zetaC2_hat[0, :-1])
    y02 = np.ravel(zetaC2_hat[1 , :-1])
    phi02 = np.ravel(zetaC2_hat[2 , :-1])

else:
    x02 = np.reshape(mpc2_vars['zeta0_0'][offset_k:], (1, iter2))
    y02 = np.reshape(mpc2_vars['zeta0_1'][offset_k:], (1, iter2))
    phi02 = np.reshape(mpc2_vars['zeta0_2'][offset_k:], (1, iter2))
    s02 = np.reshape(mpc2_vars['zeta0_3'][offset_k:], (1, iter2))
    n02 = np.reshape(mpc2_vars['zeta0_4'][offset_k:], (1, iter2))
    beta02 = np.reshape(mpc2_vars['zeta0_5'][offset_k:], (1, iter2))
    v02 = np.reshape(mpc2_vars['zeta0_6'][offset_k:], (1, iter2))
    alpha02= np.reshape(mpc2_vars['zeta0_7'][offset_k:], (1, iter2))


#beta_ref = np.reshape(np.zeros(iter1), (1, iter1))
#x_ref, y_ref, phi_ref = sysModel.Fren2CartT(s01, n_ref, beta_ref)
x2_ref, y2_ref, phi2_ref = InterpolLuT( np.ravel(s02))
x2_ref  = np.ravel(x2_ref)
y2_ref  = np.ravel(y2_ref)
phi2_ref  = np.ravel(phi2_ref)
n2_ref = np.ravel(np.zeros(np.shape(x2_ref)))

a01 = np.reshape(mpc1_vars['u0_0'][offset_k:], (1, iter1))
omg01 = np.reshape(mpc1_vars['u0_1'][offset_k:], (1, iter1))

a02 = np.reshape(mpc2_vars['u0_0'][offset_k:], (1, iter2))
omg02 = np.reshape(mpc2_vars['u0_1'][offset_k:], (1, iter2))


x1 = np.ravel(x01)
x2 = np.ravel(x02)

y1 = np.ravel(y01)
y2 = np.ravel(y02)

phi1 = np.ravel(phi01)
phi2 = np.ravel(phi02)

s1 = np.ravel(s01)
s2 = np.ravel(s02)

n1 = np.ravel(n01)
n2 = np.ravel(n02)

beta1 = np.ravel(beta01)
beta2 = np.ravel(beta02)

v1 = np.ravel(v01)
v2 = np.ravel(v02)

alpha1 = np.ravel(alpha01)
alpha2 = np.ravel(alpha02)

a1 = np.ravel(a01)
a2 = np.ravel(a02)

omg1 = np.ravel(omg01)
omg2 = np.ravel(omg02)

print("DEBUG n avg1", round(np.abs(n1).mean(), 4))
print("DEBUG n avg2", round(np.abs(n2).mean(), 4))

fig3 = mpl.pyplot.figure(figsize=(4.8, 3.6), num='zeta_c')
zetaCAx = fig3.add_subplot(2, 1, 1)
zetaCAx.stairs(x1[:], s1[:], baseline=None,label="Direct", color="lightcoral" )
zetaCAx.stairs(x2[:-1], s2[:], baseline=None,label="Lifted", color="cadetblue")
zetaCAx.stairs(x2_ref[:-1], s2[:], baseline=None,label="Track", color="black", linestyle='dashed' )
zetaCAx.set_ylim( np.hstack([x1, x2]).max(axis=0) * SCALE_BOUND, 
                  np.hstack([x1, x2]).min(axis=0) * SCALE_BOUND)
zetaCAx.set_xlim(0, np.amax(np.hstack([s1,s2]).max(axis=0)))

zetaCAx.set_ylabel("$x\,(\\mathrm{m})$")
zetaCAx.invert_yaxis() 
zetaCAx.legend(loc='upper right')
zetaCAx.grid()
# zetafAx.set_yscale('symlog')

zetaCAx = fig3.add_subplot(2, 1, 2)
zetaCAx.stairs(y1[:], s1[:], baseline=None, color="lightcoral")#,label="Direct", color="lightcoral" )
zetaCAx.stairs(y2[:-1], s2[:], baseline=None, color="cadetblue")#,label="Lifted", color="darkseagreen")
zetaCAx.stairs(y2_ref[:-1], s2[:], baseline=None, linestyle='dashed' ,color="black")#label="Track", linestyle='dashed' ,color="black")
zetaCAx.set_ylim( np.hstack([y1, y2]).max(axis=0)  + 1, 
                  np.hstack([y1, y2]).min(axis=0)  * SCALE_BOUND)
zetaCAx.set_xlim(0, np.amax(np.hstack([s1, s2]).max(axis=0)))
zetaCAx.set_xlabel('track progress $s\,(\\mathrm{m})$')
zetaCAx.set_ylabel("$y\,(\\mathrm{m})$")
# zetaCAx.legend(loc='upper right')
# zetaCAx.legend.remove()
zetaCAx.grid()
zetaCAx.invert_yaxis()

# fig3_1 = mpl.pyplot.figure(figsize=(6.4, 2.65), num='zeta_c_1')
# zetaCAx = fig3_1.add_subplot(1, 1, 1)
# zetaCAx.stairs(phi1[:], s1[:], baseline=None,label="Direct", color="lightcoral" )
# zetaCAx.stairs(phi2[:-1], s2[:], baseline=None,label="Lifted", color="cadetblue")
# zetaCAx.stairs(phi2_ref[:-1], s2[:], baseline=None,label="Track", linestyle='dashed' ,color="black")
# zetaCAx.set_ylim( np.hstack([phi1, phi2]).max(axis=0)  + 0.1, 
#                   np.hstack([phi1, phi2]).min(axis=0)  * SCALE_BOUND)
# zetaCAx.set_xlim(0, np.amax(np.hstack([s1, s2]).max(axis=0)))
# zetaCAx.set_xlabel('track progress $s\,(\\mathrm{m})$')
# zetaCAx.set_ylabel("$\\varphi\,(\\mathrm{rad})$")
# zetaCAx.legend(loc='upper right')
# zetaCAx.grid()
# zetaCAx.invert_yaxis()

# # print("\n oscillation 2nd derv:", v1_osc, v2_osc)
# fig1 = mpl.pyplot.figure(num='zeta_phi_n')
# zetaFAx = fig1.add_subplot(2, 1, 1)
# zetaFAx.stairs(phi1[:], s1[:], baseline=None, color="lightcoral" )#,label="Direct"
# zetaFAx.stairs(phi2[:-1], s2[:], baseline=None, color="cadetblue")#,label="Lifted"
# zetaFAx.stairs(phi2_ref[:-1], s2[:], baseline=None,label="Track", linestyle='dashed' ,color="black")
# zetaFAx.set_ylim( np.hstack([phi1, phi2]).max(axis=0)  + 0.1, 
#                   np.hstack([phi1, phi2]).min(axis=0)  * SCALE_BOUND)
# zetaFAx.set_xlim(0, np.amax(np.hstack([s1, s2]).max(axis=0)))
# #zetaFAx.set_xlabel('track progress $s\,(\\mathrm{m})$')
# zetaFAx.set_ylabel("$\\varphi\,(\\mathrm{rad})$")
# zetaFAx.legend(loc='upper right')
# zetaFAx.grid()
# zetaFAx.invert_yaxis()

# zetaFAx = fig1.add_subplot(2, 1, 2)
# zetaFAx.stairs(n1[:-1], s1[:], baseline=None,label="Direct", color="lightcoral" )
# zetaFAx.stairs(n2[:-1], s2[:], baseline=None,label="Lifted", color="cadetblue")
# # zetaFAx.stairs(n1_ref[:-1], s1[:], baseline=None, label="Reference", linestyle='dashed',color="black")
# zetaFAx.set_ylim( np.hstack([n1, n2]).max(axis=0) * SCALE_BOUND, 
#                  np.hstack([n1, n2]).min(axis=0) * SCALE_BOUND)
# zetaFAx.set_xlim(0, np.amax(np.hstack([s1,s2]).max(axis=0)))

# zetaFAx.set_ylabel("$n\,(\\mathrm{m})$")
# zetaFAx.invert_yaxis() 
# zetaFAx.set_xlabel('track progress $s\,(\\mathrm{m})$')
# zetaFAx.legend(loc='upper right')
# zetaFAx.grid()
# # zetafAx.set_yscale('symlog')

# zetaFAx = fig1.add_subplot(2, 1, 2)
# zetaFAx.stairs(beta1[:-1], s1[:], baseline=None,label="Direct", color="lightcoral" )
# zetaFAx.stairs(beta2[:-1], s2[:], baseline=None,label="Lifted", color="teal")
# zetaFAx.set_ylim( np.hstack([beta1, beta2]).max(axis=0) * SCALE_BOUND, 
#                   np.hstack([beta1, beta2]).min(axis=0) * SCALE_BOUND)
# zetaFAx.set_xlim(0, np.amax(np.hstack([s1, s2]).max(axis=0)))

# zetaFAx.set_xlabel('track progress $s\,(\\mathrm{m})$')
# zetaFAx.set_ylabel("$\\beta\,(\\mathrm{rad})$")
# zetaFAx.legend(loc='upper right')
# zetaFAx.grid()
# zetaFAx.invert_yaxis()
#statAx.set_yscale('symlog')

fig2 = mpl.pyplot.figure(figsize=(4.8, 3.6),num='zeta_u')
zetaUAx = fig2.add_subplot(2, 1, 1)
zetaUAx.stairs(v1[:-1], s1[:], baseline=None,label="Direct", color="lightcoral" )
zetaUAx.stairs(v2[:-1], s2[:], baseline=None,label="Lifted", color="cadetblue")
zetaUAx.set_ylim( np.hstack([v1, v2]).max(axis=0) * SCALE_BOUND, 
                  0.5 )
zetaUAx.set_xlim(0, np.amax(np.hstack([s1, s2]).max(axis=0)))

zetaUAx.set_ylabel("$v\,(\\mathrm{m\,s^{-1}})$")
zetaUAx.invert_yaxis() 
zetaUAx.legend(loc='upper right')
zetaUAx.grid()
# zetafAx.set_yscale('symlog')

zetaUAx = fig2.add_subplot(2, 1, 2)
zetaUAx.stairs(alpha1[:-1], s1[:], baseline=None, color="lightcoral" )#label="Direct"
zetaUAx.stairs(alpha2[:-1], s2[:], baseline=None, color="cadetblue") #label="Lifted",
zetaUAx.set_ylim( np.hstack([alpha1, alpha2]).max(axis=0) * SCALE_BOUND, 
                  np.hstack([alpha1, alpha2]).min(axis=0) * SCALE_BOUND)
zetaUAx.set_xlim(0, np.amax(np.hstack([s1, s2]).max(axis=0)))

zetaUAx.set_xlabel('track progress $s\,(\\mathrm{m})$')
zetaUAx.set_ylabel("$\\alpha\,(\\mathrm{rad})$")
# zetaUAx.legend(loc='upper right')
zetaUAx.grid()
zetaUAx.invert_yaxis()

fig4 = mpl.pyplot.figure(figsize=(4.8, 3.6), num='u')
UAx = fig4.add_subplot(2, 1, 1)
UAx.stairs(a1[:-1], s1[:], baseline=None,label="Direct", color="lightcoral" )
UAx.stairs(a2[:-1], s2[:], baseline=None,label="Lifted", color="cadetblue")
UAx.set_ylim( np.hstack([a1, a2]).max(axis=0) * SCALE_BOUND, 
             np.hstack([a1, a2]).min(axis=0) * SCALE_BOUND)
UAx.set_xlim(0, np.amax(np.hstack([s1, s2]).max(axis=0)))

UAx.set_ylabel("$a\,(\\mathrm{m\,s^{-2}})$")
UAx.invert_yaxis() 
UAx.legend(loc='upper right')
UAx.grid()

UAx = fig4.add_subplot(2, 1, 2)
UAx.stairs(omg1[:-1], s1[:], baseline=None, color="lightcoral" )#label="Direct",
UAx.stairs(omg2[:-1], s2[:], baseline=None, color="cadetblue")#label="Lifted",
UAx.set_ylim( np.hstack([omg1, omg2]).max(axis=0) * SCALE_BOUND, 
              np.hstack([omg1, omg2]).min(axis=0) * SCALE_BOUND)
UAx.set_xlim(0, np.amax(np.hstack([s1, s2]).max(axis=0)))

UAx.set_xlabel('track progress $s\,(\\mathrm{m})$')
UAx.set_ylabel("$\\omega\,(\\mathrm{rad\,s^{-1}})$")
# UAx.legend(loc='upper right')
UAx.grid()
UAx.invert_yaxis()

fig4.tight_layout()
fig4.savefig('u.pdf', format='pdf', bbox_inches='tight')

fig3.tight_layout()
fig3.savefig('zeta_c.pdf', format='pdf', bbox_inches='tight')

# fig3_1.tight_layout()
# fig3_1.savefig('zeta_c1.pdf', format='pdf', bbox_inches='tight')

fig2.tight_layout()
fig2.savefig('zeta_u.pdf', format='pdf', bbox_inches='tight')

# fig1.tight_layout()
# fig1.savefig('zeta_phi_n.pdf', format='pdf', bbox_inches='tight')
mpl.pyplot.show()