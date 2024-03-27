# Introduction 
This ROS node (for Noetic) focuses on solving a local path planning problem with Model predictive control (MPC).
It implements a  numerical optimization_based control strategy, achieved with the software toolboxes in Python like Casadi and Acados

# Setup
Setup instructions for POSIX-compliant OS like Linux, macOS

*Step                           : CLI command*
1. python3, pip                 : ```apt install python3, pip3```
2. install virtual environment  : ```pip install virtualenv```
                                  ```python -m venv <env. name>```
   activate virtual environment : ```source <env. name>/bin/activate```
                                         
3. python package dependencies  : ```pip install -r requirements.txt```
   casadi*, numpy, matplotlib
   acados **

4. ROS interface                 : Clone this repository in the ```cd src/apps/```                                folder of your Catkin workspace 
                                   ```catkin_make```
                                   ```source /devel/setup.bash```

*Casadi needs to be built from source for macOs. Refer https://web.casadi.org/get/

** Acados needs to be built from source for all operating systems. Refer
https://docs.acados.org/installation/