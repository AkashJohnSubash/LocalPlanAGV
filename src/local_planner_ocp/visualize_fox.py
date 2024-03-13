import numpy as np 
from numpy import pi
import casadi as ca

from local_planner_ocp.common import *

import rospy
from visualization_msgs.msg import Marker, MarkerArray
from geometry_msgs.msg import Point
from geometry_msgs.msg import PoseArray, PoseStamped
from local_planner_ocp.msg import mpc_parameters, mpc_statistics, mpc_variables, mpc_metrics
from tf.transformations import quaternion_from_euler

# Set pose
mpc_var_pub = rospy.Publisher('mpc_variables', mpc_variables ,queue_size = 1, latch=True)
mpc_param_pub = rospy.Publisher('mpc_parameters', mpc_parameters ,queue_size = 1, latch=True)
mpc_stat_pub = rospy.Publisher('mpc_statistics', mpc_statistics ,queue_size = 1, latch=True)
mpc_metrics_pub = rospy.Publisher('mpc_metrics', mpc_metrics ,queue_size = 1, latch=True)
markerArr_pub = rospy.Publisher('MarkerArray', MarkerArray, queue_size=1)
marker_pub = rospy.Publisher('Markers', Marker, queue_size=1)
poseArr_pub = rospy.Publisher('ShootingNodes', PoseArray, queue_size=1)
trackArr_pub = rospy.Publisher('WaypointTrack', MarkerArray, queue_size=1, latch=True)

def PublishMpcParams(lifted, ns):
    param = mpc_parameters()

    # Flatten matrices to 1D and publish
    param.Q_cost = np.diag(Q, k=0).tolist()
    param.Qn_cost = np.diag(Qn, k=0).tolist()
    param.R_cost = np.diag(R, k=0).tolist()
    param.Ts = T_del
    param.N = N
    param.Tf = Tf
    param.s_ref = S_REF
    param.v_ref = V_AGV_REF
    param.lifted = lifted
    param.ns = ns

    mpc_param_pub.publish(param)

def PublishMpcVars(t0, zeta0, u0):
    var = mpc_variables()

    # Flatten matrices to 1D and publish
    var.zeta0 = np.ravel(zeta0).tolist()
    var.u0 = np.ravel(u0).tolist()
    var.sim_time = t0
    mpc_var_pub.publish(var)


def PublishMpcStats(fails, sqp_time_avg, sqp_time_max, lat_dev, speed):
    stats = mpc_statistics()

    # Flatten matrices to 1D and publish
    stats.sqp_failures = fails
    stats.sqp_avg_time = sqp_time_avg
    stats.sqp_max_time = sqp_time_max
    stats.agv_avg_lat_dev = lat_dev
    stats.agv_avg_speed = speed

    mpc_stat_pub.publish(stats)

def PublishMpcMetrics(mpc_t, obsrv_t, ref_t, acados_t, qp_t, sim_t, cost):
    metrics = mpc_metrics()

    # Flatten matrices to 1D and publish
    #metrics.sqp_iter = sqp_iter
    metrics.acados_time = acados_t
    metrics.qp_time = qp_t
    metrics.sim_time = sim_t
    metrics.cost = cost
    # metrics.slacks =  np.ravel(slacks).tolist()
    # metrics.residuals = np.ravel(residuals).tolist()
    metrics.mpc_time = mpc_t
    metrics.obsrv_time = obsrv_t
    metrics.ref_time = ref_t

    mpc_metrics_pub.publish(metrics)

def PublishPoseArray(state_list, sysMod, lifted = False):
    '''Visualize shooting nodes by projecting to cartesian every MPC iter'''

    horizonArr = PoseArray()
    time_stamp = rospy.Time()
    horizonArr.header.frame_id = "map"
    horizonArr.header.stamp = time_stamp
    for k in range(0, N, 5):
        states = PoseStamped()
        states.header.frame_id = "map"
        states.header.stamp = rospy.Time()
        if(lifted):
            x, y, phi = state_list[0, k], state_list[1, k], state_list[2, k]
        else:
            x, y, phi = sysMod.Fren2CartT(state_list[0, k], state_list[1, k], state_list[2, k])
        states.pose.position.x = x
        states.pose.position.y = y
        states.pose.position.z = 0.2
        states.pose.orientation.w = 1 

        horizonArr.poses.append(states.pose)

    poseArr_pub.publish(horizonArr)

def PublishTrackWaypoints():
    '''Visualize track way points (updated only on startup)'''

    trackCircList = MarkerArray()
    # Publish waypoints as green spheres
    for i in range (np.shape(x_ref)[0]):
        #print(f"DEBUG waypoints {targ_states}")    
        waypoint = Marker()
        waypoint.header.frame_id = "map"
        waypoint.ns = "waypoint_sphere"
        
        waypoint.type = Marker.SPHERE
        waypoint.action = Marker.ADD
        waypoint.header.stamp = rospy.Time()
        waypoint.lifetime = rospy.Duration()

        # Set the scale
        waypoint.scale.x = 0.2
        waypoint.scale.y = 0.2
        waypoint.scale.z = 0.2
        
        waypoint.pose.orientation.w = 1
        waypoint.pose.position.x = x_ref[i]
        waypoint.pose.position.y = y_ref[i]
        waypoint.pose.position.z = 0.2

        # Set the color
        waypoint.color.r = 0.8
        waypoint.color.g = 0.8
        waypoint.color.b = 0.8
        waypoint.color.a = 0.3

        trackCircList.markers.append(waypoint)
    # Assign IDs to spheres in MarkerArray
    
    id = 0
    for m in trackCircList.markers:
        m.id = id    
        id = id+1
    

    trackArr_pub.publish(trackCircList)

def PublishCircleMarkers(obst_list, state_meas):

    # Publish obstacle marker as red spheres
    CircMarkList = MarkerArray()
    for i in range (np.shape(obst_list)[1]):
        sphere = Marker()
        sphere.header.frame_id = "map"
        sphere.ns = "obstacle_sphere"

        sphere.type = Marker.SPHERE
        sphere.action = Marker.ADD
        sphere.header.stamp = rospy.Time()
        # sphere.lifetime = marker

        # Set the scale of obstacles
        sphere.scale.x = 2 * SCALE_LAM * obst_list[2, i]
        sphere.scale.y = 2 * obst_list[2, i]
        # obstacle circle lies on ground (z = r)
        sphere.scale.z = 2 * obst_list[2, i]
        
        yaw = obst_list[3, i]
        #q = Quaternion()
        #Quaternion.
                        #rot about  x, y, z
        q = quaternion_from_euler(0, 0, yaw)
        #Quaternion myQuaternion    myQuaternion.setRPY( 0, 0, 0 );
        sphere.pose.orientation.x = q[0]
        sphere.pose.orientation.y = q[1]
        sphere.pose.orientation.z = q[2]
        sphere.pose.orientation.w= q[3]
        
        sphere.pose.position.x = obst_list[0, i]
        sphere.pose.position.y = obst_list[1, i]
        sphere.pose.position.z = 0

        # Set the color
        sphere.color.r = 1
        sphere.color.g = 0
        sphere.color.b = 0
        sphere.color.a = 0.8

        CircMarkList.markers.append(sphere)
    
    # Publish blue AGV position marker
    # shape_bounding_marker(CircMarkList, agv_odom) 
    shape_covering_markers(CircMarkList, state_meas)
    
    id = 0
    for m in CircMarkList.markers:
        m.id = id    
        id = id+1
    
    markerArr_pub.publish(CircMarkList)

def PublishLineObstacles(obst_list):
    
    # Publish line segments for walls
    line_list = Marker()
    line_list.header.frame_id = "map"
    line_list.ns = "planner_line"

    line_list.type = Marker.LINE_LIST
    line_list.action = Marker.ADD
    line_list.header.stamp = rospy.Time()

    # Set the scale of line
    line_list.scale.x = 0.1

    # Set the color
    line_list.color.r = 1
    line_list.color.g = 0
    line_list.color.b = 0
    line_list.color.a = 0.8
    
    for i in range (np.shape(obst_list)[1]):
        p = Point()
        p.x = obst_list[0, i]
        p.y = obst_list[1, i]
        p.z = 0
        #print(f"\n line point {i}: {p}")
        line_list.points.append(p)

    marker_pub.publish(line_list)

def shape_bounding_marker(markerList, agv_odom):
    agv = Marker()
    agv.header.frame_id = "map"
    agv.ns = "obstacle_sphere"

    agv.type = Marker.SPHERE
    agv.action = Marker.ADD
    agv.header.stamp = rospy.Time()
    agv.lifetime = rospy.Duration(secs = 1)

    # Scale the bounding figure as an ellipse or circle
    agv.scale.x = 2* rob_el_a
    agv.scale.y = 2* rob_el_b
    agv.scale.z = 1
    
    agv.pose.orientation.x = agv_odom['phi_quat'][0]
    agv.pose.orientation.y = agv_odom['phi_quat'][1]
    agv.pose.orientation.z = agv_odom['phi_quat'][2]
    agv.pose.orientation.w = agv_odom['phi_quat'][3]
    agv.pose.position.x = agv_odom['pos_x_map']
    agv.pose.position.y = agv_odom['pos_y_map']
    agv.pose.position.z = 0

    # Set the color
    agv.color.r = 0
    agv.color.g = 0
    agv.color.b = 1
    agv.color.a = 0.35
    markerList.markers.append(agv)

    
def shape_covering_markers(markerList, state_meas):
    for k in range(-1, 2):
        agv = Marker()
        agv.header.frame_id = "map"
        agv.ns = "obstacle_sphere"

        agv.type = Marker.SPHERE
        agv.action = Marker.ADD
        agv.header.stamp = rospy.Time()
        agv.lifetime = rospy.Duration(secs = 1)

        # Scale the covering circles
        agv.scale.x = 2 *  R_Kc
        agv.scale.y = 2 * R_Kc
        agv.scale.z = 1

        x_map = np.copy(state_meas)[0]
        y_map = np.copy(state_meas)[1]
        phi_map = np.copy(state_meas)[2]
        v = np.copy(state_meas)[3]/2

        agv.pose.position.x = x_map + (C_OFFSET + k * D_Kc) * ca.cos(phi_map)
        agv.pose.position.y = y_map + (C_OFFSET + k * D_Kc) * ca.sin(phi_map)
    
        # Set the color
        agv.color.r = 0
        agv.color.g = 0
        agv.color.b = 1
        agv.color.a = 0.35
        
        markerList.markers.append(agv)
