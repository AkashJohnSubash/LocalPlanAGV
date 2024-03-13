import numpy as np 
import matplotlib as mpl
from matplotlib import pyplot
import shutil

from common import *
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

mpl.rcParams.update(params)


def plotOptVars(time_stamps, traj_ST0, traj_U0, lifted=False ):
    
    '''Plot1 state and control '''
    sysModel = SysDyn()
    if lifted:
        fig1 = pyplot.figure(figsize=(16,9), num='Lifted state')
    else:
        fig1 = pyplot.figure(figsize=(16,9), num='Direct Elimination state')
    
    zetaC = fig1.add_subplot(4, 1, 1)
    zetaF = fig1.add_subplot(4, 1, 2)
    zetaU = fig1.add_subplot(4, 1, 3)  
    u = fig1.add_subplot(4, 1, 4)    
    
    time_stamps = time_stamps[Tstart_offset:]
    traj_ST0 = traj_ST0[:, Tstart_offset:]
    traj_U0 = traj_U0[:, Tstart_offset:]
    #traj_twist0 = traj_twist0[:, Tstart_offset:]

    dim_st = np.shape(traj_ST0)

    if(lifted == True):
        zetaC_idx = 0
        zetaF_idx = 3
        zetaU_idx = 6
    else:
        zetaC_idx = 0
        zetaF_idx = 0
        zetaU_idx = 3

    if(lifted):
        x = np.ravel(traj_ST0[zetaC_idx, :-1])
        y = np.ravel(traj_ST0[zetaC_idx + 1, :-1])
        phi = np.ravel(traj_ST0[zetaC_idx + 2, :-1])
        
        zetaC.set_ylabel("$\\zeta^{c}$")
        zetaC.set_ylim( np.amin(np.ravel(traj_ST0[zetaC_idx: zetaF_idx, :])) - 1, 
                        np.amax(np.ravel(traj_ST0[zetaC_idx: zetaF_idx, :])) + 1)

    else:
        zetaC_hat = np.zeros((3, dim_st[1]))

        # Project frenet state to cartesian
        #for k in range(N+1):   
        s_i = traj_ST0[0, :]
        n_i = traj_ST0[1, :]
        beta_i = traj_ST0[2, :]

        x_i, y_i, phi_i = sysModel.Fren2CartT(s_i, n_i, beta_i)
        zetaC_hat[0, : ] = np.ravel(x_i)
        zetaC_hat[1, : ] = np.ravel(y_i)
        zetaC_hat[2, : ] = np.ravel(phi_i)
                
        x = np.ravel(zetaC_hat[zetaC_idx, :-1])
        y = np.ravel(zetaC_hat[zetaC_idx + 1, :-1])
        phi = np.ravel(zetaC_hat[zetaC_idx + 2, :-1])

        zetaC.set_ylabel("$\\hat{\\zeta}^{c}$")
        zetaC.set_ylim( np.amin(np.ravel(zetaC_hat[:, :])) - 1, 
                        np.amax(np.ravel(zetaC_hat[:, :])) + 1)
        
    s = np.ravel(traj_ST0[zetaF_idx, :-1])
    n = np.ravel(traj_ST0[zetaF_idx + 1, :-1])
    beta = np.ravel(traj_ST0[zetaF_idx + 2, :-1])
    
    v = np.ravel(traj_ST0[zetaU_idx, :-1])
    alpha = np.ravel(traj_ST0[zetaU_idx + 1, :-1])
    
    vdot = np.ravel(traj_U0[0, :-1])
    alphadot = np.ravel(traj_U0[1, :-1])

    # print("DEBUG", type(x[0]), type(time_stamps[0]))
    #print("DEBUG", type(x), type(time_stamps))
    #print("DEBUG", np.shape(x), np.shape(time_stamps))
    zetaC.stairs(x, time_stamps[ :], baseline=None,label="$x$ ($\mathrm{m}$)", color="coral" )
    zetaC.stairs(y, time_stamps[ :], baseline=None,label="$y$ ($\mathrm{m}$)", color="teal")
    zetaC.stairs(phi, time_stamps[ :], baseline=None,label="$\\varphi$ ($\mathrm{rad}$)", color="#6052a2ff")
    zetaC.set_xlim(0, np.amax(time_stamps[ :]) + 0.2)
    # zetaC.set_xlabel('time (s)')
    zetaC.legend(loc='upper right')
    zetaC.grid()

    zetaF.stairs(s, time_stamps[ :], baseline=None,label="$s$ ($m$)", color="coral" )
    zetaF.stairs(n, time_stamps[ :], baseline=None,label="$n$ ($m$)", color="teal")
    zetaF.stairs(beta, time_stamps[ :], baseline=None,label="$\\beta$ ($\mathrm{rad}$)", color="#6052a2ff")
    zetaF.set_ylim( np.amin(np.ravel(traj_ST0[zetaF_idx : zetaU_idx, :-2])) - 0.2, 
                    np.amax(np.ravel(traj_ST0[zetaF_idx : zetaU_idx, :-2])) + 20)
    zetaF.set_xlim(0, np.amax(time_stamps[ :]) + 0.2)
    # zetaF.set_xlabel('time (s)')
    zetaF.set_ylabel("$\\zeta^{f}$")
    zetaF.set_yscale('symlog')
    zetaF.legend(loc='upper right')
    zetaF.grid()

    zetaU.stairs(v, time_stamps[ :], baseline=None,label="$v$ ($m s^{-1}$)", color="coral" )
    zetaU.stairs(alpha, time_stamps[ :], baseline=None,label="$\\alpha$ ($rad$)", color="teal")
    zetaU.set_ylim( np.amin(np.ravel(traj_ST0[zetaU_idx:, :-2])) - 0.2, 
                    np.amax(np.ravel(traj_ST0[zetaU_idx:, :-2])) + 0.2)
    zetaU.set_xlim(0, np.amax(time_stamps[ :]) + 0.2)
    # zetaU.set_xlabel('time (s)')
    zetaU.set_ylabel('$\\zeta^{u}$')
    zetaU.legend(loc='upper right')
    zetaU.grid()
    
    u.stairs(vdot, time_stamps[ :], baseline=None,label="$a$ ($m s^{-2}$)", color="coral" )
    u.stairs(alphadot, time_stamps[ :], baseline=None,label="$\\omega$ ($rad s^{-1}$)", color="teal")
    u.set_ylim( np.amin(traj_U0[:, :]) - 0.2,
                np.amax(traj_U0[:, :]) + 0.2)
    u.set_xlim(0, np.amax(time_stamps[ :]) + 0.2)
    u.set_xlabel('time (s)')
    u.set_ylabel('u')
    u.legend(loc='upper right')
    u.grid()

    fig1.tight_layout()
    pyplot.show()

def plotResiduals(time_stamps, traj_res, lifted= False):
    '''Plot2 residuals '''

    if lifted:
        fig2 = pyplot.figure(figsize=(16,9), num='Lifted Residuals')
    else:
        fig2 = pyplot.figure(figsize=(16,9), num='Direct Elimination Residuals')
    
    statAx = fig2.add_subplot(4, 1, 1)
    eqAx = fig2.add_subplot(4, 1, 2)
    ineqAx = fig2.add_subplot(4, 1, 3)
    compAx = fig2.add_subplot(4, 1, 4)

    time_stamps = time_stamps[Tstart_offset:]
    traj_res = traj_res[:, Tstart_offset:]

    stat_res = np.ravel(traj_res[0, :-1])
    eq_res = np.ravel(traj_res[1, :-1])
    ineq_res = np.ravel(traj_res[2, :-1])
    comp_res = np.ravel(traj_res[3, :-1])
    comp_res = np.ravel(traj_res[3, :-1])

    statAx.stairs(stat_res, time_stamps[ :], baseline=None,label="$stat$", color="coral" )
    statAx.set_ylim( np.amin(traj_res[0, :]) - 0.2, np.amax(traj_res[0, :]) + 0.2)
    statAx.set_xlim(0, np.amax(time_stamps[ :]) + 0.2)
    statAx.legend(loc='upper right')
    statAx.grid()
    
    eqAx.stairs(eq_res, time_stamps[ :], baseline=None,label="$eq$", color="teal")
    eqAx.set_ylim( np.amin(traj_res[1, :]) - 0.2, np.amax(traj_res[1, :]) + 0.2)
    eqAx.set_xlim(0, np.amax(time_stamps[ :]) + 0.2)
    eqAx.legend(loc='upper right')
    eqAx.grid()

    ineqAx.stairs(ineq_res, time_stamps[ :], baseline=None,label="$ineq$", color="plum")
    ineqAx.set_ylim( np.amin(traj_res[2, :]) - 0.2, np.amax(traj_res[2, :]) + 0.2)
    ineqAx.set_xlim(0, np.amax(time_stamps[ :]) + 0.2)
    ineqAx.legend(loc='upper right')
    ineqAx.grid()

    compAx.stairs(comp_res, time_stamps[ :], baseline=None,label="$comp$", color="steelblue")
    compAx.set_ylim( np.amin(traj_res[3, :]) - 0.2, np.amax(traj_res[3, :]) + 0.2)
    compAx.set_xlim(0, np.amax(time_stamps[ :]) + 0.2)
    compAx.set_xlabel('time (s)')
    compAx.legend(loc='upper right')
    compAx.grid()

    fig2.tight_layout()
    pyplot.show()
    
def plotCosts(time_stamps, traj_cost, traj_slack0, lifted=True):#, traj_twist0, 

    ''' Plot3 costs, slacks'''

    if lifted:
        fig3 = pyplot.figure(figsize=(16,9), num='Lifted Costs')
    else:
        fig3 = pyplot.figure(figsize=(16,9), num='Direct Elimination Costs')
    
    costAx = fig3.add_subplot(3, 1, 1) 

    time_stamps = time_stamps[Tstart_offset:]
    traj_cost = traj_cost[:, Tstart_offset:]
    traj_slack0 = traj_slack0[:, :,Tstart_offset:]
    #print("DEBUG slack", np.shape(traj_slack0))
    cost = np.ravel(traj_cost[0, :-1])
    

    costAx.stairs(cost, time_stamps[ :], baseline=None,label="$cost$", color="steelblue")
    costAx.set_ylim( np.amin(traj_cost[ :]) - 0.2, np.amax(traj_cost[ :]) + 0.2)
    costAx.set_xlim(0, np.amax(time_stamps[ :]) + 0.2)
    # costAx.set_xlabel('time (s)')
    costAx.legend(loc='upper right')
    costAx.grid()

    # plot lower slacks at 1,1 (top right)  
    sl = fig3.add_subplot(3, 1, 2)
    sl.set_ylim(np.amin(np.ravel(traj_slack0[:, 0, :])) - 0.1, 
                np.amax(np.ravel(traj_slack0[:, 0, :])) + 0.1)
    sl.set_xlim(0, np.amax(time_stamps[ :]) + 0.2)
    # sl.set_xlabel('time (s)')
    sl.set_ylabel('lower slacks')
    
    # # Stack all obstacle axes, but label only one
    sl.stairs(traj_slack0[0, 0, :-1], time_stamps[ :], baseline=None, label="obstacle breach (m)", color="lightcoral")
    for i in range(1, obst_constr_len):
        sl.stairs(traj_slack0[i, 0, :-1], time_stamps[ :], baseline=None,color="lightcoral")
    sl.stairs(traj_slack0[ obst_constr_len, 0, :-1], time_stamps[ :], baseline=None, label="$v_{gg}$ breach ($m s^{-1}$)" ,color="burlywood" )
    sl.stairs(traj_slack0[ obst_constr_len + 1, 0, :-1], time_stamps[ :], baseline=None, label="$\\omega_{gg}$ breach ($rad s^{-1}$)" ,color="steelblue" )
    sl.grid()
    sl.legend(loc='upper right')

    # Plot upper slacks (dynamics only) at (1,3) mid right  
    su = fig3.add_subplot(3, 1, 3)
    su.set_ylim(np.amin(np.ravel(traj_slack0[obst_constr_len:, 1, :])) - 0.1, 
                np.amax(np.ravel(traj_slack0[obst_constr_len:, 1, :])) + 0.1)
    su.set_xlim(0, np.amax(time_stamps[ :]) + 0.2)
    su.set_xlabel('time (s)')
    su.set_ylabel('upper slacks')

    su.stairs(traj_slack0[ obst_constr_len, 1, : -1], time_stamps[ :], baseline=None, label="$v_{gg}$ breach ($m s^{-1}$)" ,color="burlywood" )
    su.stairs(traj_slack0[ obst_constr_len + 1, 1, : -1], time_stamps[ :], baseline=None, label="$\\omega_{gg}$ breach ($rad s^{-1}$)" ,color="steelblue" )
    #UniSuAx = su.stairs([], [0], baseline=None, label="unique breach ($rad s^{-1}$)" ,color="red" )
    su.grid()
    su.legend(loc='upper right')

    pyplot.show()