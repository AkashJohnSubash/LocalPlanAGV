
from bagpy import bagreader
import csv
import re
import numpy as np
import pandas as pd
import matplotlib as mpl
from common import *
# from matplotlib import pyplot as plt
from sys import argv
import shutil

from plot_mpl import plotOptVars, plotCosts, plotResiduals

text_usetex = True if shutil.which('latex') else False
params = {
        'text.latex.preamble': r"\usepackage{gensymb} \usepackage{amsmath}",
        'axes.labelsize': 12,
        'axes.titlesize': 12,
        'legend.fontsize': 12,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'text.usetex': text_usetex,
        'font.family': 'serif'
}

SCALE_BOUND = 1.05
mpl.rcParams.update(params)

file1 = argv[1]
# file2 = argv[2]

# b1 = bagreader(file1)
# b2 = bagreader(file2)

# csvfiles = []

# # get all from MPC node
# vars1_csv   = b1.message_by_topic('/mpc_variables')
# metric1_csv   = b1.message_by_topic('/mpc_metrics')
# cmd1_csv   = b1.message_by_topic('/cmd_vel_auto')
# param1_csv   = b1.message_by_topic('/mpc_parameters')

# vars2_csv   = b2.message_by_topic('/mpc_variables')
# metric2_csv   = b2.message_by_topic('/mpc_metrics')
# cmd2_csv   = b2.message_by_topic('/cmd_vel_auto')
# param2_csv   = b2.message_by_topic('/mpc_parameters')

offset_k = 30
# mpc1_vars = pd.read_csv(vars1_csv)
# mpc1_metric = pd.read_csv(metric1_csv)
# mpc1_param = pd.read_csv(param1_csv)
# cmd1_vel   = pd.read_csv(cmd1_csv)

# mpc2_vars = pd.read_csv(vars2_csv)
# mpc2_metric = pd.read_csv(metric2_csv)
# mpc2_param = pd.read_csv(param2_csv)
# cmd2_vel   = pd.read_csv(cmd2_csv)

with open(file1, 'r') as file:
  csvreader = csv.reader(file)
  data = []
  for row in csvreader:
    data = row[0].split(": ")
    # numList = re.findall(r"[\d]*[.][\d]+", data[0])
    # requiredString = s.substring(s.indexOf("(") + 1, s.indexOf(")"));
    if(len(data)>8):
        print("\n", float(data[9].split("\n")[0]))