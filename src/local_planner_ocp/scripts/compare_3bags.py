
from bagpy import bagreader
import numpy as np
import pandas as pd
import matplotlib as mpl
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

print("DEBUG", text_usetex)
mpl.rcParams.update(params)

file1 = argv[1]
file2 = argv[2]
file3 = argv[3]

b1 = bagreader(file1)
b2 = bagreader(file2)
b3 = bagreader(file3)
csvfiles = []
# get the list of topics
#print(b1.topic_table)
SCALE_BOUND = 1.05
offset_k = 20

# get all from MPC node
vars1_csv   = b1.message_by_topic('/mpc_variables')
metric1_csv   = b1.message_by_topic('/mpc_metrics')

vars2_csv   = b2.message_by_topic('/mpc_variables')
metric2_csv   = b2.message_by_topic('/mpc_metrics')

vars3_csv   = b3.message_by_topic('/mpc_variables')
metric3_csv   = b3.message_by_topic('/mpc_metrics')

mpc1_vars = pd.read_csv(vars1_csv)
mpc1_metric = pd.read_csv(metric1_csv)

mpc2_vars = pd.read_csv(vars2_csv)
mpc2_metric = pd.read_csv(metric2_csv)

mpc3_vars = pd.read_csv(vars3_csv)
mpc3_metric = pd.read_csv(metric3_csv)

iter1 = np.shape(mpc1_vars['sim_time'])[0] - offset_k
t01 = np.ravel(mpc1_vars['sim_time'])[:]#, [iter, 1])

iter2 = np.shape(mpc2_vars['sim_time'])[0] - offset_k
t02 = np.ravel(mpc2_vars['sim_time'])[:]#, [iter, 1])

iter3 = np.shape(mpc3_vars['sim_time'])[0] - offset_k
t03 = np.ravel(mpc3_vars['sim_time'])[:]#, [iter, 1])

'''Plot state and control'''
mpc1_iter = np.arange( iter1)
mpc2_iter = np.arange( iter2)
mpc3_iter = np.arange( iter3)

s01 = np.reshape(mpc1_vars['zeta0_0'][offset_k:], (1, iter1))
n01 = np.reshape(mpc1_vars['zeta0_1'][offset_k:], (1, iter1))
alpha01 = np.reshape(mpc1_vars['zeta0_4'][offset_k:], (1, iter1))

s02 = np.reshape(mpc2_vars['zeta0_0'][offset_k:], (1, iter2))
n02 = np.reshape(mpc2_vars['zeta0_1'][offset_k:], (1, iter2))
alpha02 = np.reshape(mpc2_vars['zeta0_4'][offset_k:], (1, iter2))

s03 = np.reshape(mpc3_vars['zeta0_0'][offset_k:], (1, iter3))
n03 = np.reshape(mpc3_vars['zeta0_1'][offset_k:], (1, iter3))
alpha03 = np.reshape(mpc3_vars['zeta0_4'][offset_k:], (1, iter3))
#plotOptVars(t01, zeta0, u0, None, lifted=lift)

#print(f"1: {t01} 2: {t02}, {t03}")
n1 = np.ravel(n01)
n2 = np.ravel(n02)
n3 = np.ravel(n03)
s1 = np.ravel(s01)
s2 = np.ravel(s02)
s3 = np.ravel(s03)
alpha1 = np.ravel(alpha01)
alpha2 = np.ravel(alpha02)
alpha3 = np.ravel(alpha03)
# print(s01, "\n", s02, "\n", s03)figsize=(16,9)figsize=(7.05, 5)#figsize=(6.4, 2.8),
fig1 = mpl.pyplot.figure( num='Delay compensation schemes')
delayComp = fig1.add_subplot(2, 1, 1)
delayComp.stairs(n3[:-1], s3[:], baseline=None,label="Eu 2", color="lightcoral")
delayComp.stairs(n2[:-1], s2[:], baseline=None,label="Eu 1", color="teal")
delayComp.stairs(n1[:-1], s1[:], baseline=None,label="Eu 0",  color="#220d7eb6", alpha=0.55)
delayComp.set_ylim( 0.15,#np.hstack([n1,n2,n3]).max(axis=0), 
                np.hstack([n1,n2,n3]).min(axis=0) * SCALE_BOUND )
delayComp.set_xlim(0, np.amax(np.hstack([s1,s2,s3]).max(axis=0)))
delayComp.invert_yaxis() 
# delayComp.set_xlabel('track progress $s\,(\\mathrm{m})$')
delayComp.set_ylabel("$n\,(\\mathrm{m})$")
delayComp.legend(loc='upper right')

major_ticks = np.arange(-0.2, 0.61, 0.3)

#delayComp.set_yticks(major_ticks)
delayComp.grid()

alphaAx = fig1.add_subplot(2, 1, 2)
alphaAx.stairs(alpha3[:-1], s3[:], baseline=None,label="Eu 2", color="lightcoral" )
alphaAx.stairs(alpha2[:-1], s2[:], baseline=None,label="Eu 1", color="teal")
alphaAx.stairs(alpha1[:-1], s1[:], baseline=None,label="Eu 0", color="#220d7eb6", alpha=0.55)
alphaAx.set_ylim( np.hstack([alpha1, alpha2, alpha3]).max(axis=0) + 0.025, 
                np.hstack([alpha1, alpha2, alpha3]).min(axis=0) * SCALE_BOUND)
alphaAx.set_xlim(0, np.amax(np.hstack([s1, s2, s3]).max(axis=0)))

alphaAx.set_xlabel('track progress $s\,(\\mathrm{m})$')
alphaAx.set_ylabel("$\\alpha\,(\\mathrm{rad})$")
alphaAx.legend(loc='upper right')
alphaAx.grid()
alphaAx.invert_yaxis()


fig1.tight_layout()

fig1.savefig('comp_pred.pdf', format='pdf', bbox_inches='tight')
# # plt.plot(y)
# mpl.pyplot.savefig("test1.svg", format="svg")
fig1.tight_layout()
fig1.savefig('perdictor.pdf', format='pdf', bbox_inches='tight')
mpl.pyplot.show()