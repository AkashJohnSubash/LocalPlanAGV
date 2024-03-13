import numpy as np 
import matplotlib as mpl
from matplotlib import pyplot, animation
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib.patches import Circle, Ellipse

from local_planner_ocp.common import *
from local_planner_ocp.sys_dynamics import SysDyn


def animOptVars(traj_samples, traj_ST, traj_U, traj_slack, traj_twist, lifted=False):
    '''Plot data as animate matplotlib graph'''
    sysModel = SysDyn()
    # (nx , N+1 , mpc_iter)
    dim_st = np.shape(traj_ST)
    # dim_u = np.shape(traj_ST)
    z_st = np.ones((1, dim_st[2])) * z_const
    aW_max = np.ones((1, dim_st[2])) * A_W_MAX / CONSTR_FACT
    aW_min = np.ones((1, dim_st[2])) * A_W_MIN / CONSTR_FACT
    omgW_min = np.ones((1, dim_st[2])) * OMG_W_MIN / CONSTR_FACT
    omgW_max = np.ones((1, dim_st[2])) * OMG_W_MAX / CONSTR_FACT
    vAgv_max = np.ones((1, dim_st[2])) * V_AGV_MAX / CONSTR_FACT
    vAgv_min = np.ones((1, dim_st[2])) * V_AGV_MIN / CONSTR_FACT
    omgAgvAx_max = np.ones((1, dim_st[2])) * OMG_AGV_MAX / CONSTR_FACT
    omgAgvAx_min = np.ones((1, dim_st[2])) * OMG_AGV_MIN / CONSTR_FACT
    anim_running = True

    if(lifted == True):
        zetaC_idx = 0
        zetaF_idx = 3
        zetaU_idx = 6
    else:
        zetaC_idx = 0
        zetaF_idx = 0
        zetaU_idx = 3
    
    # Plot only the original track without repetition
    [_, xref_orig, yref_orig, _, _] = getTrack()
    ref_samples = len(xref_orig)
    z_ref = np.ones(ref_samples) * z_const
    # subsample the trajectries for faster plots
    
    # nx X N X (mpc_iter/freq)
    traj_ST = traj_ST[:, :, ::f_plot]
    # nu X N X (mpc_iter/freq)
    traj_U = traj_U[:, ::f_plot]
    # ns x (mpc_iter/frPeq)
    traj_slack = traj_slack[:, :,::f_plot]
    # Contains discretization (times, mpc_stages) 3x (mpc_iter/freq)
    traj_samples = traj_samples[:, ::f_plot]
    # Generate trajectory of twist
    traj_twist = traj_twist[:, ::f_plot]
    
    if(lifted == False):
        zetaC_hat = np.zeros((3, dim_st[1], dim_st[2]))
        zetaC_hat = zetaC_hat[:, :, ::f_plot]

        # Project frenet state to cartesian (constraint co-ordinate)
        for k in range(N+1):   
            s_i = traj_ST[0, k, :]
            n_i = traj_ST[1, k, :]
            beta_i = traj_ST[2, k, :]

            x_i, y_i, phi_i = sysModel.Fren2CartT(s_i, n_i, beta_i)
            zetaC_hat[0, k, : ] = np.ravel(x_i)
            zetaC_hat[1, k, : ] = np.ravel(y_i)
            zetaC_hat[2, k, : ] = np.ravel(phi_i)

    def init():
        time_text.set_text('')
        # return path, horizon, sphere_i


    def onClick(event):
        nonlocal anim_running
        if anim_running:
            anim.event_source.stop()
            anim_running = False
        else:
            anim.event_source.start()
            anim_running = True

    def animate(iter):
        '''Update animation'''

        # update State plot
        time_text.set_text(f'time = {traj_samples[0, iter]:.2f} s' )

        # update robot covering circles and horizon
        if(lifted):
            for k in range ( -1, 2):
                x_k = traj_ST[0, 0, iter] + k * D_Kc * ca.cos(traj_ST[2, 0, iter])
                y_k = traj_ST[1, 0, iter] + k * D_Kc * ca.sin(traj_ST[2, 0, iter])
                agv[k]._offsets3d = (float(x_k), float(y_k), z_st[0, iter: ])

            horizon.set_data(traj_ST[zetaC_idx: zetaC_idx+2, 1:, iter])
            horizon.set_3d_properties(z_st[0, iter ])

            xAx.set_data(traj_ST[zetaC_idx,  0, : iter], traj_samples[0, :iter +1] )
            yAx.set_data(traj_ST[zetaC_idx + 1,  0, : iter], traj_samples[0, :iter +1] )
            phiAx.set_data(traj_ST[zetaC_idx + 2,  0, : iter], traj_samples[0, :iter +1] )

        else:
            for k in range ( -1, 2):
                x_k = zetaC_hat[0, 0, iter] + k * D_Kc * ca.cos(zetaC_hat[2, 0, iter])
                y_k = zetaC_hat[1, 0, iter] + k * D_Kc * ca.sin(zetaC_hat[2, 0, iter])
                agv[k]._offsets3d = (float(x_k), float(y_k), z_st[0, iter:])

                xAx.set_data(zetaC_hat[zetaC_idx,  0, : iter], traj_samples[0, :iter +1] )
                yAx.set_data(zetaC_hat[zetaC_idx + 1,  0, : iter], traj_samples[0, :iter +1] )
                phiAx.set_data(zetaC_hat[zetaC_idx + 2,  0, : iter], traj_samples[0, :iter +1] )

            
            horizon.set_data(zetaC_hat[zetaC_idx: zetaC_idx+2, 1:, iter])
            horizon.set_3d_properties(z_st[0, iter ])
                
        
        # Update state plots
        
        sAx.set_data(traj_ST[zetaF_idx,  0, : iter], traj_samples[0, :iter +1] )
        nAx.set_data(traj_ST[zetaF_idx + 1,  0, : iter], traj_samples[0, :iter +1] )
        betaAx.set_data(traj_ST[zetaF_idx + 2,  0, : iter], traj_samples[0, :iter +1] )
        velAx.set_data(traj_ST[zetaU_idx,  0, : iter], traj_samples[0, :iter +1] )
        alpAx.set_data(traj_ST[zetaU_idx + 1,  0, : iter], traj_samples[0, :iter +1] )

        # Update control plot (drive command)
        accWAx.set_data(traj_U[0, :iter], traj_samples[0, :iter +1] )
        accWMaxAx.set_data(aW_max[0, :iter], traj_samples[0, :iter +1] )
        accWMinAx.set_data(aW_min[0, :iter], traj_samples[0, :iter +1] )
        
        omgWAx.set_data(traj_U[1, :iter], traj_samples[0, :iter +1] )
        omgWMaxAx.set_data(omgW_max[0, :iter], traj_samples[0, :iter +1] )
        omgWMinAx.set_data(omgW_min[0, :iter], traj_samples[0, :iter +1] )
        
        # Update ROS twist plot and other constraints
        vAgvAx.set_data(traj_twist[0, :iter], traj_samples[0, :iter +1])
        vAgvMaxAx.set_data(vAgv_max[0, :iter], traj_samples[0, :iter +1])
        vAgvMinAx.set_data(vAgv_min[0, :iter], traj_samples[0, :iter +1])
        
        omgAgvAx.set_data(traj_twist[1, :iter], traj_samples[0, :iter +1] )
        omgAgvMaxAx.set_data(omgAgvAx_max[0, :iter], traj_samples[0, :iter +1])
        omgAgvMinAx.set_data(omgAgvAx_min[0, :iter], traj_samples[0, :iter +1])

        # update cost plot
        costAx.set_data(traj_samples[1, :iter], traj_samples[0, :iter +1] )

        # Update lower slack plot (obstacle, dynamics)
        for i in range(obst_constr_len):
            obstSlAx[i].set_data(traj_slack[i, 0, :iter], traj_samples[0, :iter +1] )

        vAgvSlAx.set_data(traj_slack[ obst_constr_len, 0, : iter], traj_samples[0, :iter +1] )
        omgAgvSlAx.set_data(traj_slack[ obst_constr_len + 1, 0, : iter], traj_samples[0, :iter +1] )
            
        # Update upper slack plot (dynamics only)

        vAgvSuAx.set_data(traj_slack[ obst_constr_len, 1, : iter], traj_samples[0, :iter +1] )
        omgAgvSuAx.set_data(traj_slack[ obst_constr_len + 1, 1, : iter], traj_samples[0, :iter +1] )
        #UniSuAx.set_data(traj_slack[ obst_constr_len + 2, 1, : iter], traj_samples[0, :iter +1] )

        return path, horizon, #agv

    # Plotting routine entry point    

    # Create a figure which occupies the full screen
    fig = pyplot.figure(figsize=(15,10))

    # plot state on right (merge top and bottom right. i.e subplots 2, 3, 5, 6, 8, 9)
    ax3d = fig.add_subplot(4, 3, (8, 12), projection='3d')
    ax3d.azim = 120
    ax3d.elev = 87
    fig.add_axes(ax3d)
    
    # time field 
    time_text = ax3d.text2D(0.02, 0.95, '', transform=ax3d.transAxes)

    # reference trajectory
    pyplot.plot(xref_orig, yref_orig, z_ref, linestyle='dashed', marker = 'x', c='gray', dashes=(5, 15))
    # path
    path = ax3d.plot([], [], 'b', alpha=0.5, linewidth=0.5)[0]
    # horizon
    horizon, = ax3d.plot([], [],'x-g', alpha=0.5)

    cage_x = [-20, 25]
    cage_y = [-25, 10]
    cage_z = [0, .3]

    ax3d.set_aspect('equal')
    ax3d.set_xlim3d(left = cage_x[0], right = cage_x[1])
    ax3d.set_ylim3d(bottom = cage_y[0], top = cage_y[1])
    ax3d.set_zlim3d(bottom = cage_z[0], top = cage_z[1])

    # multiple circles around agv
    agv = [ None, None, None ]
    for k in range ( -1, 2):
        x, y, phi = init_zeta[0], init_zeta[1], init_zeta[2]
        x_k = x + k * D_Kc * ca.cos(phi)
        y_k = y + k * D_Kc * ca.sin(phi)
        agv[k] = ax3d.scatter(x_k, y_k, z_const, s = PI * rob_el_a**2 * sphere_scale, c='cornflowerblue', alpha=0.01)

    obstacle_obj = []
    # Ellipse around obstacle position
    for i in range (N_obst_max):
        for j in range(15):
            ellipse = Ellipse(xy = (float(obst_constr[i * obst_dim]), float(obst_constr[i * obst_dim + 1])), 
                              width =  2 * float(obst_constr[i * obst_dim + 2]), 
                              height =  2 * SCALE_LAM* float(obst_constr[i * obst_dim + 2]), 
                              angle = float(obst_constr[i * obst_dim + 3]) * 180 / PI, 
                              alpha = 0.04, color='lightcoral')
            
            ax3d.add_patch(ellipse)
            # stack several ellipses to add depth to obstacle viz.
            ell_o = art3d.pathpatch_2d_to_3d(ellipse, z = z_const + 0.001 * j)
            
        obstacle_obj.append(ell_o)

    ax3d.set_xlabel('X (m)')
    ax3d.set_ylabel('Y (m)')
    ax3d.set_zlabel('Z (m)')

    # State zeta_frenet at (1,1) top left
    zetaF = fig.add_subplot(4, 3, 1)
    zetaF.set_ylim( np.amin(np.ravel(traj_ST[zetaF_idx : zetaU_idx, 0, :-2])) - 0.2, 
                    np.amax(np.ravel(traj_ST[zetaF_idx : zetaU_idx, 0, :-2])) + 0.2)
    zetaF.set_xlim(0, np.amax(traj_samples[0, :]) + 0.2)
    zetaF.set_xlabel('time (s)')
    zetaF.set_ylabel("$\\zeta^{f}$")
    zetaF.set_yscale('symlog')
    sAx = zetaF.stairs([], [0], baseline=None,label="$s$ ($m$)", color="burlywood" )
    nAx = zetaF.stairs([], [0], baseline=None,label="$n$ ($m$)", color="plum")
    betaAx = zetaF.stairs([], [0], baseline=None,label="$\\beta$ ($rad$)", color="steelblue")
    zetaF.legend(loc='upper right')

    # State zeta_cartesian at (1,2) top left
    zetaC = fig.add_subplot(4, 3, 2)
    if(lifted):
        
        zetaC.set_ylim( np.amin(np.ravel(traj_ST[zetaC_idx: zetaF_idx, 0, :])) - 0.2, 
                        np.amax(np.ravel(traj_ST[zetaC_idx: zetaF_idx, 0, :])) + 0.2)
        zetaC.set_ylabel("$\\zeta^{c}$")
    else:
        zetaC.set_ylim(np.amin(np.ravel(zetaC_hat[zetaC_idx: zetaU_idx, 0, :])) - 1, 
                       np.amax(np.ravel(zetaC_hat[zetaC_idx: zetaU_idx, 0, :])) + 1)
        zetaC.set_ylabel("$\\hat{\\zeta}^{c}$")

    zetaC.set_xlim(0, np.amax(traj_samples[0, :]) + 0.2)
    zetaC.set_xlabel('time (s)')
    
    #zetaC.set_yscale('symlog')
    xAx = zetaC.stairs([], [0], baseline=None,label="$x$ ($m$)", color="burlywood" )
    yAx = zetaC.stairs([], [0], baseline=None,label="$y$ ($m$)", color="plum")
    phiAx = zetaC.stairs([], [0], baseline=None,label="$\\varphi$ ($rad$)", color="steelblue")
    zetaC.legend(loc='upper right')

    # State zeta_u at (2, 1) mid left
    zetaU = fig.add_subplot(4, 3, 4)
    zetaU.set_ylim(np.amin(np.ravel(traj_ST[zetaU_idx:, :-2])) - 0.2, 
               np.amax(np.ravel(traj_ST[zetaU_idx:, :-2])) + 0.2)
    zetaU.set_xlim(0, np.amax(traj_samples[0, :]) + 0.2)
    zetaU.set_xlabel('time (s)')
    zetaU.set_ylabel('$\\zeta^{u}$')
    velAx = zetaU.stairs([], [0], baseline=None,label="$v$ ($m s^{-1}$)", color="lightcoral" )
    alpAx = zetaU.stairs([], [0], baseline=None,label="$\\alpha$ (rad)", color="darkturquoise")
    zetaU.legend(loc='upper right')

    # Control, Twist, Slack 2D subplots
    # plot ROS cmd, Controls at (1,1) top mid
    u = fig.add_subplot(4, 3, 7)
    u.set_ylim(np.amin(np.ravel(traj_U[:, :-2])) - 0.2, 
               np.amax(np.ravel(traj_U[:, :-2])) + 0.2)
    u.set_xlim(0, np.amax(traj_samples[0, :]) + 0.2)
    u.set_xlabel('time (s)')
    u.set_ylabel('u')
    
    accWAx = u.stairs([], [0], baseline=None,label="$a$ ($m s^{-2}$)", color="salmon" )
    accWMaxAx = u.stairs([], [0], baseline=None, color="maroon")
    accWMinAx = u.stairs([], [0], baseline=None, color="maroon" )
    
    omgWAx = u.stairs([], [0], baseline=None,label="$\\omega$ ($rad s^{-1}$)", color="teal")
    omgWMaxAx = u.stairs([], [0], baseline=None,  color="darkslategray")
    omgWMinAx = u.stairs([], [0], baseline=None, color="darkslategray" )
    u.legend(loc='upper right')

    # plot ROS Twist at (2,2) mid mid 
    uDer = fig.add_subplot(4, 3, 5)
    uDer.set_ylim(np.amin(np.ravel(traj_twist[:, :-2])) - 0.2, 
                  np.amax(np.ravel(traj_twist[:, :-2])) + 0.2)
    uDer.set_xlim(0, np.amax(traj_samples[0, :]) + 0.2)
    uDer.set_xlabel('time (s)')
    uDer.set_ylabel('twist')

    vAgvAx = uDer.stairs([], [0], baseline=None,label="$v_{gg}$ ($m s^{-1}$)", color="burlywood" )
    vAgvMaxAx = uDer.stairs([], [0], baseline=None,color="darkorange")
    vAgvMinAx = uDer.stairs([], [0], baseline=None,color="darkorange" )
    
    omgAgvAx = uDer.stairs([], [0], baseline=None, label="$\\omega_{gg}$ ($rad s^{-1}$)", color="steelblue" )
    omgAgvMaxAx = uDer.stairs([], [0], baseline=None, color="navy")
    omgAgvMinAx = uDer.stairs([], [0], baseline=None,color="navy" )

    uDer.legend(loc='upper right')

    # plot costs at (4,1) bottom left 
    misc = fig.add_subplot(4, 3, 10)
    misc.set_ylim(np.amin(np.ravel(traj_samples[1, :-2])) - 0.2, 
                  np.amax(np.ravel(traj_samples[1, :-2])) + 0.2)
    misc.set_xlim(0, np.amax(traj_samples[0, :]) + 0.2)
    misc.set_xlabel('time (s)')
    misc.set_ylabel('cost')
    misc.set_yscale('symlog')

    costAx = misc.stairs([], [0], baseline=None, color="burlywood" )    

    # plot lower slacks at 1,1 (top right)  
    sl = fig.add_subplot(4, 3, 3)
    sl.set_ylim(np.amin(np.ravel(traj_slack[:, 0, :])) - 0.1, 
                np.amax(np.ravel(traj_slack[:, 0, :])) + 0.1)
    sl.set_xlim(0, np.amax(traj_samples[0, :]) + 0.2)
    sl.set_xlabel('time (s)')
    sl.set_ylabel('lower slacks')
    
    # # Stack all obstacle axes, but label only one
    obstSlAx = [None] * obst_constr_len
    obstSlAx[0] = sl.stairs([], [0], baseline=None, label="obstacle breach (m)", color="lightcoral")
    for i in range(1, obst_constr_len):
        obstSlAx[i] = sl.stairs([], [0], baseline=None,color="lightcoral")
    vAgvSlAx = sl.stairs([], [0], baseline=None, label="$v_{gg}$ breach ($m s^{-1}$)" ,color="burlywood" )
    omgAgvSlAx = sl.stairs([], [0], baseline=None, label="$\\omega_{gg}$ breach ($rad s^{-1}$)" ,color="steelblue" )
    sl.legend(loc='upper right')

    # Plot upper slacks (dynamics only) at (1,3) mid right  
    su = fig.add_subplot(4, 3, 6)
    su.set_ylim(np.amin(np.ravel(traj_slack[obst_constr_len:, 1, :])) - 0.1, 
                np.amax(np.ravel(traj_slack[obst_constr_len:, 1, :])) + 0.1)
    su.set_xlim(0, np.amax(traj_samples[0, :]) + 0.2)
    su.set_xlabel('time (s)')
    su.set_ylabel('upper slacks')

    vAgvSuAx = su.stairs([], [0], baseline=None, label="$v_{gg}$ breach ($m s^{-1}$)" ,color="burlywood" )
    omgAgvSuAx = su.stairs([], [0], baseline=None, label="$\\omega_{gg}$ breach ($rad s^{-1}$)" ,color="steelblue" )
    #UniSuAx = su.stairs([], [0], baseline=None, label="unique breach ($rad s^{-1}$)" ,color="red" )
    su.legend(loc='upper right')

    fig.canvas.mpl_connect('button_press_event', onClick)
    anim = animation.FuncAnimation(fig=fig, func=animate, 
                                   init_func=init, 
                                   frames=len(traj_samples[0, :]), 
                                   interval=refresh_ms, 
                                   repeat=False,
                                   blit=False)

    fig.tight_layout()
    pyplot.show()