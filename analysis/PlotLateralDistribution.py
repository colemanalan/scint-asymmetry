#!/usr/bin/env python

## -------------------------------------------------------------------------------------------------
## Reads in one shower and plots the lateral distribution of signal
## -------------------------------------------------------------------------------------------------

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FixedLocator

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"

from icecube import icetray, dataio, dataclasses, radcube
from icecube.dataclasses import I3Constants

import numpy as np
import os

ABS_PATH_HERE = str(os.path.dirname(os.path.realpath(__file__)))

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("input", help="Input I3 file with scintillator response")
parser.add_argument("--output", type=str, default=ABS_PATH_HERE + "/../plots/LateralDistribution.pdf", help="Name for the pdf with the plotted data")
args = parser.parse_args()

from I3Tray import I3Tray


def GetRadius(particle, pos):
    # Particle is the primary particle and pos is the detector position
    x_c = particle.pos.x
    y_c = particle.pos.y
    z_c = particle.pos.z

    nx = particle.dir.x
    ny = particle.dir.y

    abs_x_sq = (pos[0] - x_c) * (pos[0] - x_c) + (pos[1] - y_c) * (pos[1] - y_c) + (pos[2] - z_c) * (pos[2] - z_c)

    n_prod_x = nx * (pos[0] - x_c) + ny * (pos[1] - y_c) - np.sqrt(1.0 - nx * nx - ny * ny) * (pos[2] - z_c)

    return np.sqrt(abs_x_sq - n_prod_x * n_prod_x)


class PlotLDF(icetray.I3Module):
    def __init__(self, ctx):
        icetray.I3Module.__init__(self, ctx)

        self.i3scintgeo = False

    def Geometry(self, frame):
        self.i3scintgeo = frame["I3ScintGeometry"].scintgeo
        assert self.i3scintgeo

        print("Found scint geometry for", len(self.i3scintgeo), "panels")

    def DAQ(self, frame):
        particle = frame["MCPrimary"]

        amps = []
        radii = []
        times = []
        delays = []

        pulseSeriesMap = frame["SiPMRecoPulses"]
        #Iterate over all the panels in the frame
        for scintkey in pulseSeriesMap.keys():
            pulse = pulses[scintkey]

            #Extract the signal properties for this panel
            signalAmplitude = pulse.charge
            signalArrivalTime = pulse.t

            pos = self.i3scintgeo[scintkey]
            radius = GetRadius(particle, pos)

            #Store things for plotting later
            amps.append(signalAmplitude)
            times.append(signalArrivalTime)
            radii.append(radius)

            #Some math to calculate the expected arrival time of particles
            dR = pos - particle.pos
            arrival = particle.time - dR.Dot(dataclasses.I3Position(particle.dir)) / I3Constants.c
            delays.append(arrival - signalArrivalTime)

        #Convert into numpy arrays to make things easier
        amps = np.array(amps)
        times = np.array(times)
        radii = np.array(radii)
        delays = np.array(delays)

        #Start making the plots
        NRows = 1
        NCols = 2
        gs = gridspec.GridSpec(NRows, NCols, wspace=0.1, hspace=0.3)
        fig = plt.figure(figsize=(6 * NCols, 5 * NRows))

        # LDF Calculation
        ax = fig.add_subplot(gs[0])
        ax.scatter(radii / I3Units.m, amps)
        ax.set_yscale("log")
        ax.set_xlabel("Radius [m]")
        ax.set_yscale("Signal [MIP]")

        # Shower front calculation
        ax = fig.add_subplot(gs[1])
        ax.scatter(radii / I3Units.m, delays / I3Units.ns)
        ax.set_yscale("log")
        ax.set_xlabel("Radius [m]")
        ax.set_yscale("Delay w.r.t. Shower Plane [ns]")

        print("Saving plot as", args.output)
        fig.savefig(args.output, bbox_inches="tight")


tray = I3Tray()
tray.AddSegment("I3Reader", "Reader", input=args.input)
tray.AddModule(PlotLDF, "Plotter")
tray.Execute()
tray.Finish()
