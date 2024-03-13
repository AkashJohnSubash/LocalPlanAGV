import bagpy
from bagpy import bagreader
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sys import argv
from local_planner_ocp.plot_mpl import plotOptVars, plotCosts, plotResiduals
from local_planner_ocp.sys_dynamics import SysDyn
sysModel = SysDyn()

file = argv[1]
file2 = argv[2]
#file3 = argv[3]

b = bagreader(file)
csvfiles = []

b2 = bagreader(file2)
csvfiles2 = []

# b3 = bagreader(file3)
# csvfiles3 = [] 

# get all from MPC node
metric_csv   = b.message_by_topic('/mpc_metrics')
metric_csv2   = b2.message_by_topic('/mpc_metrics')
#metric_csv3   = b3.message_by_topic('/mpc_metrics')

mpc_metric = pd.read_csv(metric_csv)
mpc_metric2 = pd.read_csv(metric_csv2)
# mpc_metric3 = pd.read_csv(metric_csv3)

'''Plot state and control'''

mpc_time = np.array(np.ravel(mpc_metric['mpc_time'])).mean() *1000
obsrv_time = np.array(np.ravel(mpc_metric['obsrv_time'])).mean() *1000
tot_acados_time = np.array(np.ravel(mpc_metric['acados_time'])).mean() *1000
qp_time = np.array(np.ravel(mpc_metric['qp_time'])).mean() *1000
sim_time = np.array(np.ravel(mpc_metric['sim_time'])).mean() *1000

acados_time = tot_acados_time - (sim_time + qp_time) 
interface_time = mpc_time - acados_time - obsrv_time


mpc_time2 = np.array(np.ravel(mpc_metric2['mpc_time'])).mean() *1000
obsrv_time2 = np.array(np.ravel(mpc_metric2['obsrv_time'])).mean() *1000
tot_acados_time2 = np.array(np.ravel(mpc_metric2['acados_time'])).mean() *1000
qp_time2 = np.array(np.ravel(mpc_metric2['qp_time'])).mean() *1000
sim_time2 = np.array(np.ravel(mpc_metric2['sim_time'])).mean() *1000

acados_time2 = tot_acados_time2 - sim_time2 - qp_time2 
interface_time2 = mpc_time2 - acados_time2 -obsrv_time2

# mpc_time3 = np.array(np.ravel(mpc_metric3['soln_time'])).mean() *1000
# tot_acados_time3 = np.array(np.ravel(mpc_metric3['tot_time'])).mean() *1000
# qp_time3 = np.array(np.ravel(mpc_metric3['qp_time'])).mean() *1000
# sim_time3 = np.array(np.ravel(mpc_metric3['sim_time'])).mean() *1000

# acados_time3 = tot_acados_time3 - sim_time3 - qp_time3 
# interface_time3 = mpc_time3 - acados_time3

timings = ('direct', 'lifted')
tests = {'observer' : np.array([obsrv_time, obsrv_time2]),
         'predictor' : np.array([sim_time, sim_time2]),
         'MPC' : np.array([qp_time + acados_time, qp_time2 + acados_time2]),
         'cython interface' : np.array([interface_time, interface_time2]),}

colours = ('cadetblue', '#210d7e5f', 'lightcoral', 'lightgray')
# Plot bar graphs
width = 0.8
fig, ax = plt.subplots(figsize=(4.8, 3.6))

bottom = np.zeros(2)
print(timings, tests)
idx = 0
for test, num_get in tests.items():
    print(test, num_get)
    p = ax.bar(timings, num_get, width, color=colours[idx],label=test, bottom=bottom)
    bottom += num_get
    idx = idx + 1
    ax.bar_label(p, label_type='center')
ax.set_ylabel('Timings (ms)')
#ax.set_title('timings for acados python interface')
ax.legend(loc='upper right', ncol=2)

fig.tight_layout()
fig.savefig('timings.pdf', format='pdf', bbox_inches='tight')
plt.show()

