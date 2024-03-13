# Introduction 
This repository focuses solving a local motion control problem with Model predictive control (MPC).
It implements a  numerical optimization_based control strategy, achieved with the software Casadi and Acados toolboxes in Python 

# Setup
Setup instructions for POSIX-compliant OS like Linux, macOS

*Step                           : CLI command*
1. python3, pip                 : ```apt install python3, pip3```
2. install virtual environment  : ```pip install virtualenv```
                                  ```python -m venv <env. name>```
   activate virtual environment : ```source <env. name>/bin/activate```
                                         
3. python package dependencies  : ```pip install -r requirements.txt```
   casadi, numpy, matplotlib
   acados **

4. ROS interface  ***           : Clone this repository in the ```cd src/apps/``` folder of your Catkin workspace 
                                   ```catkin_make```
                                   ```source /devel/setup.bash```

# Branches
1. 'main' runs MPC with agv forward simulation from integrator (point mass) in acados (ROS installation not required)
2. 'main_ros' runs MPC  as a ROS node 'local_planner_ocp'.
   It interacts with the application stack launched with ROS simulation.launch (for Gazebo sim)
   or compact_gg.launch (real vehicle)

** Acados needs to be built from source for all operating systems. 
Refer https://docs.acados.org/installation/

*** Use linux for ROS support and executables

