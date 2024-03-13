#!/usr/bin/env python
import argparse
import rospy

from local_planner_ocp.acados_ocp import AcadosCustomOcp
from local_planner_ocp.planner import plan_ocp
from local_planner_ocp.visualize_mpl import animOptVars
from local_planner_ocp.plot_mpl import plotOptVars, plotCosts, plotResiduals

def argparse_init():
   '''
   Initialization for command line switches
   '''
   parser = argparse.ArgumentParser()
   parser.add_argument("-m", "--MatPlotLib", action='store_true', help="Display animations")
   parser.add_argument("-lf", "--lifted", action = 'store_true', help = "Use Lifted Frenet ocp")
   return parser

if __name__ == '__main__':
   rospy.init_node('ocp_main')
   rospy.loginfo("ROS entry node found")
   parser = argparse_init()
   args = parser.parse_args()

   # specify solver parameters in required format
  
   custom_ocp = AcadosCustomOcp()
   custom_ocp.setup_acados_ocp(lifted=args.lifted)
      
   # # Compute and publish OCP controls as ROS Twist
   traj_sample, traj_ST, traj_U, traj_slack, twist_traj, res_traj, cost_traj = plan_ocp(custom_ocp)

   # Plot controls and state over simulation period
   if args.MatPlotLib:
      #animOptVars(traj_sample, traj_ST, traj_U, traj_slack, twist_traj, lifted=args.lifted)
      plotOptVars(traj_sample[0, :], 
                  traj_ST[:, 0, :], 
                  traj_U, 
                  twist_traj,
                  lifted=args.lifted)
      
      plotCosts(traj_sample[0, :], 
               cost_traj,
               traj_slack,
               lifted=args.lifted)
      
      plotResiduals(traj_sample[0, :],
                  res_traj,
                  lifted=args.lifted)