import rospy
import tf
from tf.transformations import euler_from_quaternion
import math, numpy as np

from local_planner_ocp.common import InterpolLuT
from local_planner_ocp.common import d, obst_constr, obst_imagine, N_obst_max, obst_dim, PI, SCALE_LAM, SCALE_R
from local_planner_ocp.visualize_fox import PublishCircleMarkers
from local_planner_ocp.sys_dynamics import findClosestS_num
class Measurements():

    def __init__(self):
        self.listener = tf.TransformListener()
        self.odom_agv = {'pos_x_map'     :0.0,
                         'pos_y_map'     :0.0,
                         'phi_map'     :0.0,
                         'phi_quat'    : [0, 0, 0, 0],
                         'v_agv'          :0.0,
                         'phi_dot'      :0.0}
        
        self.st_meas = [0.0, 0.0, 0.0, # x, y, phi in map frame
                        0.0, 0.0]      # v, alpha in wheel frame
        self.prev_phi = 0
        
        self.st_quat = [0, 0, 0, 0]
        self.pos = [0.0, 0.0, 0.0]
        self.rot = [0.0, 0.0, 0.0, 0.0]
        self.n_circles = N_obst_max
        self.cbk_count = 1
        
    def odom_callback(self, odom_data):
        '''Pre-process recieved odometry'''
        
            
        (self.pos, self.rot) = self.listener.lookupTransform('/map', '/base_link', rospy.Time(0))
        self.odom_agv['pos_x_map'] =  self.pos[0]
        self.odom_agv['pos_y_map'] =  self.pos[1]
        self.odom_agv['phi_quat'] =  self.rot
        (_, _, yaw) = euler_from_quaternion (self.rot)
        self.odom_agv['phi_map'] = yaw

        self.odom_agv['v_agv'] = odom_data.twist.twist.linear.x
        self.odom_agv['phi_dot'] = odom_data.twist.twist.angular.z
        
        self.get_state_estimate(self.odom_agv)
    
    def get_state_estimate(self, odometry):
        '''Compute state from odometry'''
                
        # x, y, phi remains in map frame
        self.st_meas[0] = odometry['pos_x_map']
        self.st_meas[1] = odometry['pos_y_map']
        phi = odometry['phi_map']
        phi_corr = self.denormalize(np.copy(phi))
        self.st_meas[2] = phi_corr

        # compute alpha, v in wheel frame
        vcos_alp = odometry['v_agv']
        vsin_alp = d * odometry['phi_dot']
        self.st_meas[3] = math.sqrt(vcos_alp **2 + vsin_alp **2)
        self.st_meas[4] = math.atan2(vsin_alp, vcos_alp)

        self.prev_phi = np.copy(phi_corr)

    def denormalize(self, phi):
        '''correct non continoues angles due to EUL -> QUAT -> EUL normalization'''
        phi_c = 0
        if abs(phi - self.prev_phi) > PI:
            # correct -ve normalization from quat to euler #  check sign
            if(phi - self.prev_phi) < -PI:
                
                phi_c = phi + 2 * PI
            # correct +ve normalization from quat to euler
            else:
                 
                phi_c = phi - 2 * PI

        else:
            phi_c = phi   
        
        return phi_c
    
    def obstacle_callback(self, obst_data):
        global obst_constr
        
        circles = obst_data.circles
        # Number of detected obstacle
        self.n_circles = len(circles)

        circ_list = np.zeros((obst_dim, self.n_circles))

        # for i in range(0, self.n_circles ):
        #     if(i< N_obst_max): #Limits max obstacles from real scan 
        #         self.cbk_count = self.cbk_count +1
        #         x_o , y_o , phi_o = circles[i].center.x, circles[i].center.y, 0
        #         # orient obstacles parallel to track
        #         #s0 = findClosestS_num(x_o, y_o, phi_o)
        #         #_, _, gamma_phi =InterpolLuT(s0)
        #         #gamma_phi = 0
        #         # visualize obstacle
        #         radius = max(circles[i].radius * SCALE_R, 0.7)
        #         circ_list[0, i] = circles[i].center.x
        #         circ_list[1, i] = circles[i].center.y
        #         circ_list[2, i] = radius
        #         circ_list[3, i] = gamma_phi
        #         # restrict constraints to N_obst_max
        #         obst_constr[i * obst_dim ]     = circles[i].center.x
        #         obst_constr[i * obst_dim + 1]  = circles[i].center.y
        #         obst_constr[i * obst_dim + 2 ] = radius
        #         obst_constr[i * obst_dim + 3 ] =  gamma_phi
                
        # Inject imaginary obstacles into detector
        circ_list = np.zeros((4, N_obst_max))
        # visualize imaginary obstacle
        for i in range(N_obst_max):
            circ_list[0, i] = np.copy(obst_imagine[i * obst_dim])
            circ_list[1, i] = np.copy(obst_imagine[i * obst_dim + 1])
            circ_list[2, i] = np.copy(obst_imagine[i * obst_dim + 2])
            circ_list[3, i] = np.copy(obst_imagine[i * obst_dim + 3])

            obst_constr[i * obst_dim ]     = np.copy(obst_imagine[i * obst_dim])
            obst_constr[i * obst_dim + 1]  = np.copy(obst_imagine[i * obst_dim + 1])
            obst_constr[i * obst_dim + 2 ] = np.copy(obst_imagine[i * obst_dim + 2])
            obst_constr[i * obst_dim + 3 ] = np.copy(obst_imagine[i * obst_dim + 3])
            

        PublishCircleMarkers(circ_list, self.st_meas)
        return 