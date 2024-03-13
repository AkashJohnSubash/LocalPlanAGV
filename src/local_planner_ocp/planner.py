import rospy
from nav_msgs.msg import Odometry
from geometry_msgs.msg import TwistStamped

import numpy as np
from time import time, sleep
import casadi as ca

from local_planner_ocp.common import *
from local_planner_ocp.measurement import Measurements
from local_planner_ocp.sys_dynamics import Cart2Fren_num
from local_planner_ocp.visualize_fox import PublishPoseArray, PublishTrackWaypoints
from local_planner_ocp.visualize_fox import PublishMpcParams, PublishMpcStats, PublishMpcVars, PublishMpcMetrics
from obstacle_detector.msg import Obstacles

cmd_vel_auto_msg = TwistStamped()

def plan_ocp( ocp_wrapper):
    
    '''Motion control of AGV using OCPs'''

    rospy.loginfo(" Simulation started.")
    agv = Measurements()
    rospy.Subscriber('odom_pred', Odometry, agv.odom_callback)
    rospy.Subscriber('obstacles', Obstacles, agv.obstacle_callback)
    pub = rospy.Publisher('cmd_vel_auto', TwistStamped, queue_size=1)
    
    # dimensions
    nx = ocp_wrapper.nx
    nu = ocp_wrapper.nu
    ns = ocp_wrapper.ns

    PublishMpcParams(ocp_wrapper.lifted, ocp_wrapper.ns)
    # initialize iteration variables
    t_sim = 0
    times = np.array([[0]])
    mpc_iter = 0
    cost = 0
    fail = 0
    v_agv = 0
    phi_dt = 0

    # Define states, controls, and slacks as coloumn vectors
    zeta_0 = np.copy((ocp_wrapper.zeta_0))
    u_0 = np.copy(ref_u)

    # # Convert 1D array (nx,) to 3D (nx,1, N+1), and (nu,) to 3D (nu,1, N)
    zeta_N = ca.repmat(np.reshape(zeta_0, (nx,1)), 1, N+1)
    # nx x N+1 x mpc_iter (to plot state during entire Horizon)
    state_traj = ca.repmat(zeta_0, 1, N+1)
    # nu x mpc_iter
    control_traj = ca.repmat(u_0, 1, 1)
    # 3 x mpc_iter (iteration time, target index, cost)
    sample_traj = ca.repmat(np.array([0, 0]).T, 1)
    # 1 x mpc_iter
    cost_traj  = ca.repmat(cost, 1, 1)
    
    PublishTrackWaypoints()
    
    # Control loop entry point
    for i in range(Nsim):

        # Store previous iterate data for plots
        if LOGGING:
            print(f'\nSoln.{mpc_iter}: \t \u03b6F1: {np.round(ocp_wrapper.zeta_N[:, 1], 2).T} at iter {i}')
            cost_traj = np.concatenate( (cost_traj, np.reshape((cost), (1,1))), axis = 1) 
            state_traj = np.dstack((state_traj, ocp_wrapper.zeta_N))
            control_traj = np.concatenate((control_traj, 
                                            np.reshape(ocp_wrapper.u_mem[:, 0], (nu, 1))), 
                                            axis = 1)
            PublishPoseArray(ocp_wrapper.zeta_N, ocp_wrapper.sysModel, lifted = ocp_wrapper.lifted)
            PublishMpcVars(t_sim, np.reshape(ocp_wrapper.zeta_N[:, 1], (nx, 1)), np.reshape(ocp_wrapper.u_mem[:, 0], (nu, 1) ))

        t0 = time()

        # update track reference
        end = ocp_wrapper.cost_update_ref(zeta_N, ref_u)
        
        # update obstacle position, radius
        ocp_wrapper.update_parameters()
        t1 = time()

        if (end or fail > 10):
            print(f"\nMPC Exit {fail}!")
            break
       
        # Observe Frenet state vector from Cartesian measurement vector
        zetaC_meas = np.copy(agv.st_meas[0 : 3])
        zetaU_meas = np.copy(agv.st_meas[3 : 5])
        zetaF_est = Cart2Fren_num(np.copy(zetaC_meas), zeta_N, lifted = ocp_wrapper.lifted )
        t2 = time()

        # Forward simulate kinematics
        ocp_wrapper.delay_compensate(zetaC = zetaC_meas, zetaF = zetaF_est, zetaU = zetaU_meas)

        # Solve the OCP with updated state and controls
        status =  ocp_wrapper.solve_ocp()
        if(status !=0):
            fail = fail + 1
            print("Optima not found")
        
        acados_t, qp_t, sim_t  = ocp_wrapper.get_timings()
        t_sim = round(t_sim + T_del, 3)
        t3 = time()

        mpc_time = t3 - t0
        obsrv_time = t2 - t1
        ref_time = t1 - t0

        times = np.vstack(( times, mpc_time))
        mpc_iter = mpc_iter + 1

        zeta_N = ocp_wrapper.zeta_N

        v_agv, phi_dt = ocp_wrapper.compute_twist(zeta_N[:, 1])
                   
        # activate For Sim-i-L and H-i-L
        publish_cmd_auto(v_agv, phi_dt, pub)
        PublishMpcMetrics(mpc_time, obsrv_time, ref_time, acados_t, qp_t, sim_t, cost)
        delay_ms = SLEEP_SEC - (mpc_time)
        if(delay_ms > 0):
            sleep(delay_ms)
    
    sqp_max_sec = round(np.array(times).max(), 3)
    sqp_avg_sec = round(np.array(times).mean(), 3)

    if(ocp_wrapper.lifted):
        lat_dev = round(np.abs(state_traj[4, 0, :]).mean(), 4)
        w_speed = round(np.abs(state_traj[6, 0, :]).mean(), 4)
    else:
        lat_dev = round(np.abs(state_traj[1, 0, :]).mean(), 4)
        w_speed = round(np.abs(state_traj[3, 0, :]).mean(), 4)

    print(f'\nNumber of NLP failures {fail}')
    print(f'Max. SQP time\t\t: {sqp_max_sec * 1000} ms')
    print(f'Avg. SQP time\t\t: {sqp_avg_sec * 1000} ms')
    print(f'Avg. speed \t\t: {w_speed} m/s')
    lat_dev =0 
    agv_speed =0
    PublishMpcStats(fail, sqp_avg_sec, sqp_max_sec, lat_dev, agv_speed)

    return sample_traj, state_traj, control_traj, cost_traj



def publish_cmd_auto(v_agv, phi_dt, publisher):

    # compute and publish the twist from wheel velocity and orientation
    cmd_vel_auto_msg.twist.linear.x = v_agv 
    cmd_vel_auto_msg.twist.linear.y = 0
    cmd_vel_auto_msg.twist.angular.z = phi_dt
    publisher.publish(cmd_vel_auto_msg) 