import rospy

import numpy as np
import os
from pathlib import Path
import casadi as ca
from typing import Union

track="LMS_Track3.txt"

def getTrack():
    track_file = os.path.join(str(Path(__file__).parent), "tracks/",track)
    array=np.loadtxt(track_file)
    sref = array[:,0]
    xref = array[:,1]
    yref = array[:,2]
    phiref = array[:,3]
    kapparef = array[:,4]
    return sref, xref, yref, phiref, kapparef

# helper functions
def DM2Arr(dm):
    return np.array(dm.full())

def get_norm_2(diff):
    '''return euclidean norm '''
    
    norm = ca.sqrt(diff.T @ diff)
    return norm

def get_2norm_2(diff):
    '''return sqared euclidean norm '''
    
    norm = (diff.T @ diff)
    return norm

def get_2norm_W(diff, W):
    '''return weighted norm '''
    
    norm = (diff.T @ W @diff)
    return norm


def get_norm_W(diff, W):
    '''return squared weighted norm '''    
    
    norm = ca.sqrt(diff.T @ W @diff)
    return norm

def InterpolLuT( s: Union[ca.MX, float]):
    '''Interpolate curve x, y, phi based on longitudinal progress s
    <-- xref_s : position (x) on reference curve interpol function
    <-- yref_s : position (y) on reference curve interpol function
    <-- phiref_s : heading (phi) on reference curve interpol function '''

    x_ref_curve = ca.interpolant("x_ref", "bspline", [s_ref], x_ref)
    y_ref_curve = ca.interpolant("y_ref", "bspline", [s_ref], y_ref)
    phi_ref_curve = ca.interpolant("phi_ref", "bspline", [s_ref], phi_ref)

    return x_ref_curve(s), y_ref_curve(s), phi_ref_curve(s)

def normalize_angle(angle):
    '''normalize angle between -PI and PI radians '''

     	
    # normalize if lesser than -pi and greater than pi 
    angle = ca.if_else(angle <= -PI, angle + 2 * PI, angle)
    angle = ca.if_else(angle >= PI, angle - 2 * PI, angle)

    return angle

LOGGING = True



# Global variables
SD_UNIT = rospy.get_param('/wheels')['sd_unit1']
AGV_DYN = rospy.get_param('/speed_limits')
SLEEP_SEC = rospy.get_param('/odom')['fixed_frequency_time']

d  = SD_UNIT['x']
INF = 1e5
PI = 3.1415927

# furthest distance robot allowed to  deviate from track
MIN_DIST = 4

# timing parameters
T_del = 0.06                # time between shooting nodes
N = 100                     # number of shooting nodes
Tf = T_del * 100 * 1

Tsim = 200
Nsim = int(Tsim * N / Tf)

S_END = 62.5

# get params. for AGV dynamics. bounds
CONSTR_FACT = 0.90
CONSTR_FACT2 = 0.80

V_AGV_MAX = AGV_DYN['vx_max'] * CONSTR_FACT2
V_AGV_MIN = AGV_DYN['vx_min'] * CONSTR_FACT
OMG_AGV_MAX = AGV_DYN['vyaw_max'] * CONSTR_FACT
OMG_AGV_MIN = AGV_DYN['vyaw_min'] * CONSTR_FACT

# Control
ref_u = np.array([0, 0])
# get params. for control bounds

A_W_MAX = SD_UNIT['max_acc'] * CONSTR_FACT
A_W_MIN = -SD_UNIT['max_acc'] * CONSTR_FACT
OMG_W_MAX = SD_UNIT['max_steering_angle_rate'] * CONSTR_FACT
OMG_W_MIN = -SD_UNIT['max_steering_angle_rate'] * CONSTR_FACT

# Constraint
obst_constr_len = 6
dyn_constr_len = 3
                    #x, y, phi, s,  n, beta, v, alpa

init_zeta = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]) 
           
# # K covering circle radius and seperation
K = 3
R_Kc = 1
D_Kc = ca.sqrt(2) * R_Kc

C_OFFSET = 0.254
#x_o, y_o, rad_o, phi_o
N_obst_max = 2
obst_constr = ca.repmat([20, 20, 1, 0], N_obst_max)
obst_imagine = [14.0, -0.3, 0.5, 0,
                3.5, 5.5, 0.5,  -4.125,
                15, 15, 1, 0,
                15, 15, 1, 0,
                15, 15, 1, 0] 

# # a = b = r for bounding circle
SCALE_LAM = 1.5
SCALE_R = 1
rob_el_a = 3.914/ca.sqrt(2)
rob_el_b = 1.115/ca.sqrt(2)

obst_dim = 4

# matplotlib parameters
Tstart_offset = 5
f_plot = 2
refresh_ms = 1
sphere_scale = 50 
z_const = 0.130

[s_ref, x_ref, y_ref, phi_ref, kappa_ref] = getTrack()
length = len(s_ref)
pathlength = s_ref[-1]

V_AGV_REF = 0.75
S_REF = 8

#v traj
Q = np.diag([ 1e-8, 2.5e1, 1e-8, 2e2, 1e1])
R =  np.diag([5e0, 2.5e1])
Qn = np.diag([ 1e1, 2.5e1, 1e-8, 1e-8, 5e0])

if LOGGING :
    get_N = N

else:
    get_N = 1