import casadi as ca

from common import *
from sys_dynamics import SysDyn

class Model:
    def __init__(self):
        self.x = None
        self.u = None


class Integrator:

    def rk4_explicit(kin_fp, state, ctrl, Dt):
        ''' Forward simulate kinematics using explicit runge-kutta 4 integrator'''
        k1 = kin_fp(state, ctrl)
        k2 = kin_fp(state + Dt * k1/2, ctrl)
        k3 = kin_fp(state + Dt * k2/2, ctrl)
        k4 = kin_fp(state + Dt * k3, ctrl)
        zeta_next = state + (Dt /6) * (k1 + 2 * k2 + 2 * k3 + k4)
        
        return zeta_next
class Ocp:
    def __init__(self):
        self.model = Model()
        self.cost = None
        self.constraints = {'lbxu' : None,
                            'ubxu' : None,
                            'lbg'  : None,
                            'ubg'  : None }

class IpoptCustomOcp:
    
    def __init__(self):

        self.nx = 0
        self.nu = 0
        self.ny = 0
        self.ns = N_obst_max  + 2
    
        self.ocp = None, 
        self.solver = None, 
        self.integrator = None
        
        self.zeta_0 = None
        self.zeta_N = None

        self.yref = None
        self.yref_N = None
        self.nxi = 0
        self.lifted = False

    def setup_ipopt_ocp(self, lifted = False):
        '''Formulate IPOPT OCP'''
        sysModel = SysDyn()
        # generates optimal trajectory using an NMPC cost function
        if(lifted):
            kin_expl, zeta, u, _, _, _, kin_fp = sysModel.SetupOde(lifted=True)
            self.nxi = 3
            self.lifted = True

            # construct lifted state from initial guess
            x0, y0, phi0 = sysModel.Fren2CartT(init_zeta[3], init_zeta[4], init_zeta[5])
            zetaC_0 = np.array([x0, y0, phi0])
            self.zeta_0 = np.append(zetaC_0, init_zeta[3:])
        else:
            kin_expl, zeta, u, _, _, _, kin_fp = sysModel.SetupOde()
            self.zeta_0 = np.copy(init_zeta[3:])
       
        # set dimensions
        self.nx = nx = zeta.numel()
        self.nu = nu = u.numel()

        self.zeta_N = ca.repmat(np.reshape(self.zeta_0, (self.nx,1)), 1, N+1)
        self.u_N = ca.repmat(ref_u, 1, N)

        ocp = Ocp()
        ocp.model.x = zeta
        ocp.model.u = u

        # state, control, param. vector over horizon 
        zeta_n = ca.MX.sym('zeta_n', nx, N + 1)
        u_n = ca.MX.sym('u_n', nu, N)
        # state refence only in Ipopt
        self.yref = ca.MX.sym('yref', nx + nu, N + 1)
        self.yref_N = ca.MX.sym('yref_N', nx, 1)
        
        # Formulate cost function (LLS)
        cost_fn = 0
        # Stage cost                 
        for k in range(N):
            zeta_k = zeta_n[:, k]
            u_k = u_n[:, k]
            cost_fn = cost_fn + ((zeta_k[self.nxi:] - self.yref[self.nxi:-2, k]).T @ Q @ (zeta_k[self.nxi:] - self.yref[self.nxi:-2,k]) + 
                                 (u_k - self.yref[-2:]).T @ R @ (u_k - self.yref[-2:]))
        # # Terminal cost
        # zeta_k = zeta_n[:, N]
        # cost_fn = cost_fn + ((zeta_k[self.nxi:] - self.yref_N[self.nxi:, 0]).T @ Qn @ (zeta_k[self.nxi:] - self.yref_N[self.nxi:,0]))

        cost_fn = cost_fn * .25
        
        # Formulate continuity constraint
        cont_constr_eqn = []
        for k in range(N):
            zeta_k = zeta_n[:, k]
            u_k = u_n[:, k]
            zeta_nxt = zeta_n[:, k+1]
            zeta_intgr = Integrator.rk4_explicit(kin_fp, zeta_k, u_k, T_del)
            cont_constr_eqn = ca.vertcat(cont_constr_eqn, zeta_nxt - zeta_intgr) 
        cont_constr_len = cont_constr_eqn.shape[0] 

        # Formulate path inequality constraint (convex ?, Non-Linear)
        obst_constr_eqn = []
        
        # Elliptical bot with circular obstacles
        
        for i in range (N_obst_max):
            for k in range(N):
                phi = zeta_n[2, k]
                Rot_2 = ca.vertcat( ca.horzcat(ca.cos(phi), -ca.sin(phi)),
                                    ca.horzcat(ca.sin(phi), ca.cos(phi)))
                
                obst_pos = obst_constr[i*obst_dim : i*obst_dim +2]
                obst_r =  obst_constr[i*obst_dim +2]

                gamma =  obst_r + rob_el_a 
                beta =   obst_r + rob_el_b

                Ax =  ca.vertcat(   ca.horzcat(1/gamma**2, 0),
                                    ca.horzcat(0,          1/beta**2 ))
                
                Sigma = Rot_2.T @ Ax @ Rot_2
                dist_constr = get_2norm_W((zeta_n[0: 2, k] - obst_pos), Sigma) - 1
                
                obst_constr_eqn = ca.vertcat(obst_constr_eqn, dist_constr)
        
        obst_constr_len = obst_constr_eqn.shape[0]

        
        # Fomulate constraint on AGV velocity (convex ?, Non-linear)
        dyn_constr_eqn = []
        for k in range(N):
            v = zeta_n[3, k]
            alpha = zeta_n[4, k]
            v_agv = v * ca.cos(alpha) 
            dyn_constr_eqn = ca.vertcat(dyn_constr_eqn , (v_agv))
        
        # Fomulate constraint on AGV angular velocity (convex ?, Non-linear)
        for k in range(N):
            v = zeta_n[3, k]
            alpha = zeta_n[4, k]
            phi_dt = (v/d) * ca.sin(alpha) 
            dyn_constr_eqn = ca.vertcat(dyn_constr_eqn , (phi_dt))

        dyn_constr_len = dyn_constr_eqn.shape[0]

        ineq_constr_eqn = []
        ineq_constr_eqn = ca.vertcat(cont_constr_eqn, 
                                  obst_constr_eqn, 
                                  dyn_constr_eqn)

        st_size = nx * (N)
        u_size = nu * N
        
        # Define bounds on decision variables
        lbxu = ca.DM.zeros((st_size +nx + u_size, 1)) 
        ubxu = ca.DM.zeros((st_size +nx + u_size, 1))

        # State bounds
        lbxu[0:  st_size: nx] = -ca.inf;             ubxu[0: st_size: nx] = ca.inf     # x lower, upper bounds
        lbxu[1:  st_size: nx] = -ca.inf;             ubxu[1: st_size: nx] = ca.inf     # y bounds
        lbxu[2:  st_size: nx] = -ca.inf;             ubxu[2: st_size: nx] = ca.inf     # z bounds
        lbxu[3:  st_size: nx] = -ca.inf;             ubxu[3: st_size: nx] = ca.inf     # v bounds
        lbxu[4:  st_size: nx] = -ca.inf;             ubxu[4: st_size: nx] = ca.inf     # alpha bounds

        # Control bounds
        lbxu[: : nu] = -A_W_MAX;      ubxu[ : : nu] = A_W_MAX
        lbxu[1 : : nu] = -OMG_W_MAX;    ubxu[1 : : nu] = OMG_W_MAX
        
        # Bounds on equality, inequality constraints
        lbg = ca.DM.zeros((st_size)+ (N*N_obst_max) + (N*2))
        ubg = ca.DM.zeros((st_size)+ (N*N_obst_max) + (N*2))

        # Bounds on continuity constraints (equality)
        lbg[0 : st_size] = 0;                     ubg[0 : st_size] = 0

        # Bounds on path constraints (inequality)
        lbg[st_size : st_size + (N)*N_obst_max] = 0;  ubg[st_size : st_size + (N)*N_obst_max] = ca.inf

        # AGV velocity constraint: V_AGV_MIN < v_agv < V_AGV_MAX
        lbg[st_size + (N*N_obst_max): st_size + (N*N_obst_max) + N]   = -V_AGV_MAX 
        ubg[st_size + (N*N_obst_max): st_size + (N*N_obst_max) + N]   = V_AGV_MAX

        #AGV Angular velocity constraint: phi_DT_MIN < omega_agv < phi_DT_MAX
        lbg[st_size + (N*N_obst_max) + N: st_size + (N*N_obst_max) + (N*2)]   = -OMG_AGV_MAX     
        ubg[st_size + (N*N_obst_max) + N: st_size + (N*N_obst_max) + (N*2)]   = OMG_AGV_MAX

        ocp.constraints['lbxu'] = lbxu
        ocp.constraints['ubxu'] = ubxu
        ocp.constraints['lbg'] = lbg
        ocp.constraints['ubg'] = ubg

        # Configure NLP solver

        OPT_variables = ca.vertcat( zeta_n.reshape((-1, 1)),  u_n.reshape((-1, 1)) )

        nlp_prob = {'f': cost_fn, 
                    'x': OPT_variables, 
                    'g': ineq_constr_eqn, 
                    'p': ca.reshape(self.yref, ((self.nx+self.nu)*(N+1), 1))}

        opts = {'ipopt' : { 'max_iter': 100,
                            'print_level': 0, 
                            'acceptable_tol': 1e-3,
                            'acceptable_obj_change_tol': 1e-3, 
                            'linear_solver' :'mumps', 
                            'nlp_scaling_method' : 'none'},

                'print_time': 0, 
                'jit' : False,
                'compiler' : 'shell',
                'jit_options' : { 'verbose': True, 'flags' : ['-O2']},
                'jit_cleanup' : True,
                }

        self.ocp = ocp
        self.solver = ca.nlpsol('solver', 'ipopt', nlp_prob, opts)

        return True
    
    def solve_and_sim(self):
        '''Solve the OCP with multiple shooting, and forward sim. with RK4'''

        status = 0

        #print(f"DEBUG {s_N[:, 1]} , \n {self.s0}")
        # TODO Fix closed loop state estim feedback
        zeta_N = np.copy(self.zeta_N)
        u_N = np.copy(self.u_N)
        
        # Optimization vars vector for Direct Multiple Shooting
        w0 = ca.vertcat(    (ca.reshape(zeta_N, self.nx* (N+1), 1)),
                            (ca.reshape(u_N, self.nu* N, 1)))
        
        # print("\nDEBUG optim vars shape!!!\t", w0.shape)
        # print("DEBUG reference shape\t", self.yref.shape)
        # print("DEBUG constraint shape\t", self.ocp.constraints['lbxu'].shape)
        sol = self.solver(x0=w0, 
                          p=ca.reshape(self.yref,(self.nx+self.nu)*(N+1), 1),
                                       #ca.reshape(self.yref_N,(self.nx, 1))),
                          lbx=self.ocp.constraints['lbxu'],
                          ubx=self.ocp.constraints['ubxu'], 
                          lbg=self.ocp.constraints['lbg'],  
                          ubg=self.ocp.constraints['ubg'] )

        zeta_N = ca.reshape(sol['x'][ : self.nx * (N+ 1)],  (self.nx, N+1))
        u_N = ca.reshape(sol['x'][self.nx * (N+ 1) : ],  (self.nu, N))
        
        self.zeta_0 = (zeta_N[:, 1])
        self.zeta_N = ca.horzcat( zeta_N[:, 1:], ca.reshape(zeta_N[:, -1], self.nx, 1)) 

        if self.solver.stats()["success"] != True:
            #raise Exception(f'IPOPT failed. Rerun with \'print_level\' > 3 for verbose log')
            status = 1
        return status

    def cost_update_ref(self, zeta_0, u_ref = None):
        '''Update the target/ current state'''
        

        self.yref = []
        
        if(self.lifted):
            s0 = zeta_0[3]
            sref =  s0 + S_REF
            for j in range(N+1):
                sref_j = s0 + (sref - s0) * j /N
                ref = np.array([0, 0, 0, sref_j, 0, 0, V_AGV_REF, 0, 0, 0]).reshape((self.nx + self.nu, 1))
                self.yref = ca.vertcat(self.yref, ca.DM(ref))
            
            # sref_j = sref
            # ref = np.array([0, 0, 0,sref_j, 0, 0, V_AGV_REF, 0]).reshape((self.nx, 1))
            # self.yref_N = ca.DM(ref)
        
        else:
            s0 = zeta_0[0]
            sref =  s0 + S_REF
            for j in range(N+1):
                sref_j = s0 + (sref - s0) * j /N
                ref = np.array([sref_j, 0, 0, V_AGV_REF, 0, 0, 0]).reshape((self.nx -self.nxi + self.nu, 1))
                self.yref = ca.vertcat(self.yref, ca.DM(ref))

            # sref_j = sref
            # ref = np.array([sref_j, 0, 0, V_AGV_REF, 0]).reshape((self.nx -self.nxi, 1))
            # self.yref_N = ca.DM(ref)

        if s0 >= S_MAX:
            return True
    
    def get_ineq_slack(self): 
    # TODO Not supported currently with Ipopt solver wrapper
        # curr_slack_s_L(), curr_slack_s_U() of Ipopt not available in Casadi release versions.
        Sl = np.zeros((self.ns, 2))
        return ca.reshape(Sl, (self.ns, 2))
    
    def update_parameters(self):
        # TODO Obstacle parameter updates not supported with Ipopt solver wrapper
        pass

    def get_cost(self):
        # TODO Not supported currently with Ipopt solver wrapper
        return 0.0
    
    def get_residuals(self):
        # TODO Not supported currently with Ipopt solver wrapper
        pass
