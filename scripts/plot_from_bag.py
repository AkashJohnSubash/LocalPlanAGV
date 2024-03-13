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

b = bagreader(file)
csvfiles = []
# get the list of topics
print(b.topic_table)

# get all from MPC node
param_csv   = b.message_by_topic('/mpc_parameters')
vars_csv   = b.message_by_topic('/mpc_variables')
metric_csv   = b.message_by_topic('/mpc_metrics')
cmd_csv   = b.message_by_topic('/cmd_vel_auto')

param_vars = pd.read_csv(param_csv)
mpc_vars = pd.read_csv(vars_csv)
mpc_metric = pd.read_csv(metric_csv)
cmd_vel   = pd.read_csv(cmd_csv)

iter = np.shape(mpc_vars['sim_time'])[0]
t0 = np.ravel(mpc_vars['sim_time'])[:]#, [iter, 1])
lift = param_vars['lifted'][0]
'''Plot state and control'''

#print('DEBUG', file)
#print('DEBUG', cmd_vel)

v_agv = np.ravel(cmd_vel['twist.linear.x'])
omg_agv = cmd_vel['twist.angular.z']

#print('DEBUG', v_agv, np.shape(v_agv))

n_sample = np.shape(v_agv)[0]
sum = 0
for i in range(1, n_sample-1):
    sum = sum + np.square(v_agv[i+1] - 2*v_agv[i] + v_agv[i-1])
v_osc = np.sqrt(sum/n_sample)

if(lift == False):
    s0 = np.reshape(mpc_vars['zeta0_0'], (1, iter))
    n0 = np.reshape(mpc_vars['zeta0_1'], (1, iter))
    beta0 = np.reshape(mpc_vars['zeta0_2'], (1, iter))
    v0 = np.reshape(mpc_vars['zeta0_3'], (1, iter))
    alpha0 = np.reshape(mpc_vars['zeta0_4'], (1, iter))
    
    zeta0 = np.vstack((s0, n0, beta0, v0, alpha0))

else:
    x0 = np.reshape(mpc_vars['zeta0_0'], (1, iter))
    y0 = np.reshape(mpc_vars['zeta0_1'], (1, iter))
    phi0 = np.reshape(mpc_vars['zeta0_2'], (1, iter))
    s0 = np.reshape(mpc_vars['zeta0_3'], (1, iter))
    n0 = np.reshape(mpc_vars['zeta0_4'], (1, iter))
    beta0 = np.reshape(mpc_vars['zeta0_5'], (1, iter))
    v0 = np.reshape(mpc_vars['zeta0_6'], (1, iter))
    alpha0 = np.reshape(mpc_vars['zeta0_7'], (1, iter))

    zeta0 = np.vstack((x0, y0, phi0, s0, n0, beta0, v0, alpha0))

a0 = np.reshape(mpc_vars['u0_0'], (1, iter))
omega0 = np.reshape(mpc_vars['u0_1'], (1, iter))

u0 = np.vstack((a0, omega0))

ns = int(param_vars['ns'][0])

'''Plot costs and slacks'''
slack_lower = np.reshape(mpc_metric[f'slacks_{0}'], (1, iter))
for i in range(1, ns):
    sl_0 = np.reshape(mpc_metric[f'slacks_{i}'], (1, iter))
    slack_lower = np.vstack((slack_lower, sl_0))

slack_upper = np.reshape(mpc_metric[f'slacks_{ns}'], (1, iter))
for i in range(ns +1 , 2*ns):
    sl_0 = np.reshape(mpc_metric[f'slacks_{i}'], (1, iter))
    slack_upper = np.vstack((slack_upper, sl_0))

slack_traj = np.reshape(np.vstack((slack_lower[:, 0], slack_upper[:, 0])), (17,2))
for j in range(1, np.shape(slack_lower)[1]):
    slack_j = np.reshape(np.vstack((slack_lower[:, j], slack_upper[:, j])), (17,2))
    slack_traj = np.dstack((slack_traj, slack_j))

cost_traj = np.reshape(mpc_metric[f'cost'], (1, iter))
#print("DEBUG", np.shape(slack_traj))

'''Plot residuals'''
stat_res = np.reshape(mpc_metric['residuals_0'], (1, iter))
eq_res = np.reshape(mpc_metric['residuals_1'], (1, iter))
ineq_res = np.reshape(mpc_metric['residuals_2'], (1, iter))
compl_res = np.reshape(mpc_metric['residuals_3'], (1, iter))
res = np.vstack((stat_res, eq_res, ineq_res, compl_res))

#lat_dev = round(np.abs(n0).mean(), 4)
avg_stat = round(np.abs(stat_res).mean(), 4)
avg_eq = round(np.abs(eq_res).mean(), 4)
avg_ineq = round(np.abs(ineq_res).mean(), 4)
avg_comp = round(np.abs(compl_res).mean(), 4)
agv_speed = round((v_agv).mean(), 3)
agv_n = round((n0).mean(), 3)

#print(f'Avg. lateral deviation\t: {lat_dev} m')
print(f'Avg. stationarity residual \t\t: {avg_stat}')
print(f'Avg. equality residual \t\t: {avg_eq}')
print(f'Avg. inequality residual \t\t: {avg_ineq}')
print(f'Avg. complementarity residual \t\t: {avg_comp}')
print(f'2nd. deriv v \t\t: {v_osc}')

print(f'Avg. v_agv \t\t: {agv_speed}')
print(f'Avg. n \t\t: {agv_n}')

plotOptVars(t0, zeta0, u0, None, lifted=lift)
plotCosts(t0, cost_traj, slack_traj, lifted=lift)
plotResiduals(t0, res, lifted=lift)