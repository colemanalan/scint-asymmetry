import numpy as np
import math
import scipy.stats as stats
import statistics as st
import matplotlib.pyplot as plt
import pandas as pd

import argparse
import os

ABS_PATH_HERE = str(os.path.dirname(os.path.realpath(__file__)))

parser = argparse.ArgumentParser()
parser.add_argument("--output", type=str, default=ABS_PATH_HERE + "/../plots/Parameters_16_0.pdf", help="Name for the pdf with the plotted data")
args = parser.parse_args()

data = pd.read_csv('Parameters_16_0.csv', sep = ";")
print(type(data))

data = np.array(data)

Zen = data[:,0]
A = data[:,1]
B = data[:,2]
C = data[:,3]

print("done")

fig = plt.figure(figsize=(6 , 5 ))

plt.scatter(Zen, B)
plt.scatter(Zen, C)
plt.scatter(Zen, A, color="blue")
fig.set_yscale("log")

print("done")

plt.grid()
plt.xlabel("Depth")
plt.ylabel("Number of particles")



fig.savefig(args.output, bbox_inches="tight")