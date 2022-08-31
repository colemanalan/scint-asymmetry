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
print(Zen, A, B, C)

minimum = np.amin([np.amin(A),np.amin(C),np.amin(C)])

A= A + np.abs(minimum)
B= B + np.abs(minimum)
C= C + np.abs(minimum)

print(minimum)
print(A,B,C)



fig = plt.figure(figsize=(6 , 5 ))

plt.scatter(Zen, B, color="darkgoldenrod", alpha = 0.8)
plt.scatter(Zen, C, color= "darkslategrey", alpha=0.9)
plt.scatter(Zen, A, color="mediumpurple", alpha=0.7)

plt.legend([r"$A$",r"$B$",r"$C$"])
plt.grid()
plt.xlabel(r"Zenith angle $\theta$ ")
plt.ylabel("Parameter")
plt.yscale("log")
plt.title(r"Parameters of the fitting in terms of the zenith angle for $E=10^{16.0}$ eV")

#Can not use log scale since some parameters are negative.



fig.savefig(args.output, bbox_inches="tight")