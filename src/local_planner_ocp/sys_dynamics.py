import casadi as ca

from local_planner_ocp.common import *

'-------------------Global variables (casADi symbolics)---------------------------'

# State
x = ca.MX.sym('x')                     # in meters
y = ca.MX.sym('y')                     # in meters
phi = ca.MX.sym('phi')                 # in radians
s = ca.MX.sym('s')                     # in meters    
n = ca.MX.sym('n')                     # in meters
beta = ca.MX.sym('beta')               # in radians
v = ca.MX.sym('v')                     # wheel velocity in m/s 
alpha = ca.MX.sym('alpha')             # wheel turning angle in radians

# State derivative
xdot = ca.MX.sym("xdot")
ydot = ca.MX.sym("ydot")
phidot = ca.MX.sym("phidot")
sdot = ca.MX.sym("sdot")
ndot = ca.MX.sym("ndot")
betadot = ca.MX.sym("betadot")
vdot = ca.MX.sym("vdot")
alphadot = ca.MX.sym("alphadot")

# Control u 
aW = ca.MX.sym('aW');             # wheel acc. in m/s^2 
omgW = ca.MX.sym('omgW')          # wheel turning rate in radians/s

u = ca.vertcat( aW, omgW)

zeta_DEF = ca.vertcat(s, n, beta, v, alpha)
zeta_LF = ca.vertcat(x, y, phi, s, n, beta, v, alpha)

class SysDyn():
    
    def __init__(self):
        self.n_samples = 0
        self.solver = None

    def SetupOde(self, lifted = False):

        '''ODEs for Lifted-Frenet (LF)  / Direct-Elimination-Frenet (DEF) kinematics
        <-- kin_expl symbolic kinematic ODE (explicit) vector
        <-- kin_impl symbolic kinematic ODE (implicit) vector
        <-- zeta symbolic state vector
        <-- u symbolic control vector
        <-- dyn_expl symbolic dynamics (explicit) vector '''
        
        kapparef_s = ca.interpolant("kapparef_s", "bspline", [s_ref], kappa_ref)
        
        v_agv = v * ca.cos(alpha) 
        omg_agv = v * ca.sin(alpha) / d
        proj = kapparef_s(s) * n

        # kinematics for zeta_c
        xdot = v_agv * ca.cos(phi)
        ydot = v_agv * ca.sin(phi)        
        phidot = omg_agv

        # kinematics for zeta_f
        sdot = (v_agv * ca.cos(beta)) / (1 - kapparef_s(s) * n)
        ndot = v_agv * ca.sin(beta)
        betadot = (omg_agv) - kapparef_s(s) * sdot

        # zeta_u propogation for control rate limitation 
        vdot = aW
        alphadot = omgW
        
        if(lifted):
            kin_expl = ca.vertcat(xdot, ydot, phidot, sdot, ndot, betadot, vdot, alphadot)
            zeta = zeta_LF
            
        else:
            kin_expl = ca.vertcat(sdot, ndot, betadot, vdot, alphadot)
            zeta = zeta_DEF

        dyn_expl = ca.vertcat(v_agv, omg_agv, proj)
                
        return kin_expl, zeta, u, dyn_expl



    def Fren2CartT(self, s, n, beta):
        ''' Frenet to Cartesian transform
        <-- x : position (x) projection on curve w.r.t map
        <-- y : position (y) projection on curve w.r.t map
        <-- phi : AGV heading (phi) in map '''

        gamma_x, gamma_y, gamma_phi  = InterpolLuT(s)
        x = gamma_x - n * ca.sin(gamma_phi)
        y = gamma_y + n * ca.cos(gamma_phi)
        phi = gamma_phi + beta  

        return x, y, phi

def Cart2Fren_num(zetaC, zeta_N, lifted=False):
    ''' Cartesian to Frenet projection 
    (numerical approx of closed form transform)
    <-- s : longitudinal progress along reference curve
    <-- n : lateral deviation from reference curve
    <-- beta : angular deviation from reference curve '''

    x = zetaC[0]
    y = zetaC[1]
    phi = zetaC[2]

    s0 = findClosestS_num(zetaC, zeta_N, lift=lifted)       

    gamma_x, gamma_y, gamma_phi = InterpolLuT(s0)
    gamma_phi = float(gamma_phi)
    s = s0
    n = (y - float(gamma_y)) * ca.cos(gamma_phi) - (x - float(gamma_x)) * ca.sin(gamma_phi)
    beta = phi - gamma_phi
    
    return np.array([float(s), float(n), float(beta)])
    
def findClosestS_num(zetaC, zeta_N, lift=False):
    '''compute closes point by traversing entire track'''
    
    x = zetaC[0]
    y = zetaC[1]

    if(lift):
        s1 = zeta_N[3, 1]
    else:
        s1 = zeta_N[0, 1]

    # find 2 closest points, and interpolate to approximate s
    idxmindist = findClosestUniquePoint(x, y, s1, x_ref, y_ref, s_ref)
    idxmindist2 = findClosestNeighbour(x, y, x_ref, y_ref, idxmindist)
    t = findProjection(x, y, x_ref,y_ref, s_ref, idxmindist, idxmindist2)
    s0 = (1 - t) * s_ref[idxmindist] + t * s_ref[idxmindist2] 
    
    return s0

def findProjection(x, y, xref, yref, sref, idxmindist, idxmindist2):
    '''compute interpolation coefficients'''

    vabs = abs(sref[idxmindist] - sref[idxmindist2])
    vl = np.empty(2)
    u = np.empty(2)
    vl[0] = xref[idxmindist2] - xref[idxmindist]
    vl[1] = yref[idxmindist2] - yref[idxmindist]
    u[0] = x-xref[idxmindist]
    u[1] = y-yref[idxmindist]
    t = (vl[0] * u[0] + vl[1] * u[1])/vabs/vabs
    
    return t

def findClosestUniquePoint(x, y, s1, xref, yref, sref):
    '''approximation of the optimization problem with track crossing''' 
    
    mindist = MIN_DIST 
    idxmindist = 0
    for i in range(xref.size):
        dist = dist2D(x, xref[i], y, yref[i])
        if dist < mindist and s1 - MIN_DIST <= sref[i] and sref[i] <= s1 +MIN_DIST:
            mindist = dist
            idxmindist = i

    return idxmindist

def findClosestNeighbour(x, y, xref, yref, idxmindist):
    '''identify if vehicle is ahead of behind closest point'''
    
    if idxmindist != xref.size-1:    
        distBefore = dist2D(x, xref[idxmindist - 1], y, yref[idxmindist - 1])
        distAfter = dist2D(x, xref[idxmindist + 1], y, yref[idxmindist + 1])
        if(distBefore < distAfter):
            idxmindist2 = idxmindist - 1
        else:
            idxmindist2 = idxmindist + 1
        if(idxmindist2 < 0):
            idxmindist2 = xref.size - 1
    else:
        idxmindist2 = 0
    return idxmindist2


def dist2D(x1, x2, y1, y2):
    '''return euclidian distance'''
    
    return np.sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2))