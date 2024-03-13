import casadi as ca
from common import *

'-------------------Global symbolic variables---------------------------'

# State
x = ca.MX.sym('x')                     # in meters
y = ca.MX.sym('y')                     # in meters
phi = ca.MX.sym('phi')                 # in radians
s = ca.MX.sym('s')                     # in meters    
n = ca.MX.sym('n')                     # in meters
beta = ca.MX.sym('beta')               # in radians
v = ca.MX.sym('v')                     # wheel velocity in m/s 
alpha = ca.MX.sym('alpha')             # wheel angle in radians

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
omgW = ca.MX.sym('omgW')          # wheel ang vel. in radians/s

u = ca.vertcat(aW, omgW)

zeta_DEF = ca.vertcat(s, n, beta, v, alpha)
zeta_dt_DEF = ca.vertcat( sdot, ndot, betadot, vdot, alphadot)

zeta_LF = ca.vertcat(x, y, phi, s, n, beta, v, alpha)
zeta_dt_LF = ca.vertcat(xdot, ydot, phidot, sdot, ndot, betadot, vdot, alphadot)
class SysDyn():

    def __init__(self):
        
        self.n_samples = 0
        self.solver = None
    
    def SetupOde(self, lifted=False):
        '''ODEs for Lifted-Frenet (LF)  / Direct-Elimination-Frenet (DEF) kinematics
        <-- kin_expl symbolic kinematic ODE (explicit) vector
        <-- kin_impl symbolic kinematic vector
        <-- zeta symbolic state vector
        <-- u symbolic control vector
        <-- dyn_expl symbolic dynamics (explicit) vector '''
        
        kapparef_s = ca.interpolant("kapparef_s", "bspline", [s_ref], kappa_ref)
        
        v_agv = v * ca.cos(alpha) 
        omg_agv = v * ca.sin(alpha) / d
        proj = kapparef_s(s) * n
    
        # Rate of change of cartesian state (algebraic)
        x, y, phi = self.Fren2CartT(s, n, beta)
        xdot = v_agv * ca.cos(phi)
        ydot = v_agv * ca.sin(phi)
        phidot = omg_agv

        # Rate of change frenet state (ODE)
        sdot = (v_agv * ca.cos(beta)) / (1 - kapparef_s(s) * n)
        ndot = v_agv * ca.sin(beta)
        betadot = (omg_agv) - kapparef_s(s) * sdot

        # extended state variable
        vdot = aW
        alphadot = omgW
        
        if(lifted):
            kin_expl = ca.vertcat(xdot, ydot, phidot, sdot, ndot, betadot, vdot, alphadot)
            kin_impl = zeta_dt_LF - kin_expl
            zeta = zeta_LF
            zeta_dt = zeta_dt_LF
           

        else:
            kin_expl = ca.vertcat(sdot, ndot, betadot, vdot, alphadot)
            kin_impl = zeta_dt_DEF - kin_expl
            zeta = zeta_DEF
            zeta_dt = zeta_dt_DEF
        
        kin_fp = ca.Function('f', [zeta, u], [kin_expl])
        dyn_expl = ca.vertcat(v_agv, omg_agv, proj)
                
        
        return kin_expl, zeta, u, dyn_expl, kin_impl, zeta_dt, kin_fp


    def Fren2CartT(self, s, n, beta):
        ''' Frenet to Cartesian transform
        <-- x : position (x) projection w.r.t reference curve
        <-- y : position (y) projection w.r.t reference curve
        <-- phi : heading (phi) projection w.r.t reference curve '''

        gamma_x, gamma_y, gamma_phi  = InterpolLuT(s)
        x = gamma_x - n * ca.sin(gamma_phi)
        y = gamma_y + n * ca.cos(gamma_phi)
        phi = gamma_phi + beta 
        # phi = normalize_angle(phi)
        return x, y, phi
    
    def setupNlp(self):
        ''' setup NLP, solver for Cartesian to Frenet transform 
        <-- solver : solver object for transform '''

        self.n_samples = np.shape(x_ref)[0]

        s = ca.MX.sym("s", 1)
        W = ca.MX.sym("W", self.n_samples)
        cost = 0
        Wi_sum = 0
        SiWi_sum = 0
        
        curve = np.array([x_ref.T, y_ref.T])
        for i in range(self.n_samples):
            cost += W[i] * ca.norm_2(ca.vertcat(x,y) - curve[:,i])
            Wi_sum += W[i]
            SiWi_sum += s_ref[i] * W[i]

        self.solver = ca.nlpsol("solve_s", "ipopt", 
                        {"f": cost, "g": ca.vcat([Wi_sum - 1, SiWi_sum - s]), 
                            "x": ca.vertcat(s, W), "p": ca.vertcat(x,y)}, 
                            {"ipopt.print_level": 0})
        

    def Cart2FrenT(self, zetaC):
        ''' Transform cartesian position and heading to Frenet estimate for reference curve
        --> zetaC  : numpy array of Cartesian state (x, y ,phi)
        '''
        
        s_opt = self.solver( lbx = np.zeros(self.n_samples + 1),
                             ubx = s_ref[-1] * np.ones(self.n_samples + 1),
                             lbg = np.zeros(2),
                             ubg = np.zeros(2),
                             p = zetaC[:2])['x'][0]

        gamma_x, gamma_y, gamma_phi  = InterpolLuT(s_opt)

        n =  (zetaC[1] - gamma_y) * ca.cos(gamma_phi) - (zetaC[0] - gamma_x) * ca.sin(gamma_phi) 
        beta = zetaC[2] - gamma_phi
        #beta = normalize_angle(beta)

        return s_opt, n, beta
    
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
        if dist < mindist and s1 - MIN_DIST <= sref[i] and sref[i] <= s1 +MIN_DIST: #and dist3<mindist: and sref[i] <= s1 +MIN_DIST
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