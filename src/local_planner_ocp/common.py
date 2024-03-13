# import rospy
# from geometry_msgs.msg import TwistStamped

import numpy as np
import os
from pathlib import Path
import casadi as ca
from typing import Union
# import visualize_mpl as vizMpl

track="LMS_Track3.txt"

def getTrack():
    track_file = os.path.join(str(Path(__file__).parent), "tracks/", track)
    array=np.loadtxt(track_file)
    sref = array[:,0]
    xref = array[:,1]
    yref = array[:,2]
    psiref = array[:,3]
    kapparef = array[:,4]
    return sref, xref, yref, psiref, kapparef

[s_ref, x_ref, y_ref, phi_ref, kappa_ref] = getTrack()
length = len(s_ref)
pathlength = s_ref[-1]

# helper functions
def DM2Arr(dm):
    return np.array(dm.full())

M_SQRT1_2=0.70710678118654752440

def get_norm_2(diff):
    
    norm = ca.sqrt(diff.T @ diff)

    return norm

def get_2norm_2(diff):
    
    norm = (diff.T @ diff)

    return norm

def get_2norm_W(diff, W):
    
    norm = (diff.T @ W @diff)

    return norm


def get_norm_W(diff, W):
    
    norm = ca.sqrt(diff.T @ W @diff)

    return norm

def InterpolLuT(s: Union[ca.MX, float]):
    '''Interpolate curve x, y, phi based on longitudinal progress
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
    
# Global variables
SLEEP_SEC = 0.06

d  = 1.03
INF = 1e5
PI = 3.14159265358979

# timing parameters
T_del = 0.06               # time between steps in seconds
N = 100                     # number of shooting nodes
Tf = N * T_del * 1

Tsim = 50
Nsim = int(Tsim * N / Tf)

# State 
# get params. for AGV dynamics. bounds
CONSTR_FACT = 0.90
S_MAX = 100
V_AGV_MAX = 1 * CONSTR_FACT
V_AGV_MIN = -1 * CONSTR_FACT
OMG_AGV_MAX =  0.5 * CONSTR_FACT
OMG_AGV_MIN = -0.5 * CONSTR_FACT
C_OFFSET = 0.254
# Control
ref_u = np.array([0, 0])
# get params. for control bounds
A_W_MAX = 0.5 * CONSTR_FACT
A_W_MIN = -0.5 * CONSTR_FACT
OMG_W_MAX = 0.8 * CONSTR_FACT
OMG_W_MIN = -0.8 * CONSTR_FACT

# Constraint
obst_constr_len = 5
dyn_constr_len = 3
                    #x, y, phi, s,  n, beta, v, alpa
init_zeta = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ])
#init_zeta = np.array([0.08, 0.0, 0.0, 0.08, 0.0, 0.0, 0.27, 0.01 ])       

# # # K covering circle radius and seperation
K = 3
R_Kc = 1
D_Kc = ca.sqrt(2) * R_Kc

# # a = b = r for bounding circle
SCALE_LAM = 1.5
rob_el_a = 2.914/ca.sqrt(2)#3.914/ca.sqrt(2)#
rob_el_b = 1.115/ca.sqrt(2)#2.115/ca.sqrt(2)#

#x_o, y_o, rad_o, phi_o
obst_constr = ([-12.5, -0.75, 1, PI/2,
                -8, 10, PI/4, 0,
                20, 20, 1, 0,
                20, 20, 0.5, 0])

obst_dim = 4 # TODO remove hardcode
N_obst_max = int(np.shape(obst_constr)[0]/obst_dim)
# print(f"DEBUG {N_obst_max}")

#  Plot variables
#  sampling frequency of data for faster plotting
Tstart_offset = 0
f_plot = 10
refresh_ms = 1
sphere_scale = 10 #TODO make dependant on map size. (10000/ 20 obst)
z_const = 0.1

V_AGV_REF = 0.75
S_REF = 8

#v traj
Q = np.diag([ 1e-8, 2.5e1, 1e-8, 5e1, 1e1])
R =  np.diag([5e0, 2.5e1])
Qn = np.diag([ 1e1, 2.5e1, 1e-8, 1e-8, 5e0])