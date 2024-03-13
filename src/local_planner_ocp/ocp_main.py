#!/usr/bin/env python
import argparse

from acados_ocp import AcadosCustomOcp
from ipopt_ocp import IpoptCustomOcp
from planner import plan_ocp
from visualize_mpl import animOptVars
from plot_mpl import plotOptVars

def argparse_init():
   '''
   Initialization for cmd line args
   '''
   parser = argparse.ArgumentParser()
   parser.add_argument("-m", "--MatPlotLib", action = 'store_true', help = "Display animations")
   parser.add_argument("-ip", "--IPOPT", action = 'store_true', help = "Use IPOPT solver")
   parser.add_argument("-lf", "--lifted", action = 'store_true', help = "Use Lifted Frenet ocp")
   return parser

if __name__ == '__main__':
   parser = argparse_init()
   args = parser.parse_args()

   # Compute and publish OCP controls as ROS Twist
   if args.IPOPT:
      custom_ocp = IpoptCustomOcp()
      custom_ocp.setup_ipopt_ocp(lifted=args.lifted)
      traj_sample, traj_ST, traj_U, traj_slack = plan_ocp(custom_ocp)

   else:
      custom_ocp = AcadosCustomOcp()
      custom_ocp.setup_acados_ocp(lifted=args.lifted)
      traj_sample, traj_ST, traj_U, traj_slack = plan_ocp(custom_ocp)

   # Plot controls and state over simulation period
   if args.MatPlotLib:
      animOptVars(traj_sample, traj_ST, traj_U, lifted=args.lifted)
      plotOptVars(traj_sample[0, :], 
                  traj_ST[:, 0, :], 
                  traj_U, 
                  lifted=args.lifted)