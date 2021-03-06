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
from icecube.icetray import I3Units

import numpy as np
import os

ABS_PATH_HERE = str(os.path.dirname(os.path.realpath(__file__)))

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("input", nargs="+", help="Input I3 file with scintillator response")
parser.add_argument("--output", type=str, default=ABS_PATH_HERE + "/../plots/LateralDistribution.pdf", help="Name for the pdf with the plotted data")
args = parser.parse_args()

from I3Tray import I3Tray

def dR(pos,core):
    return pos-core

def Radius(pos, core, dirUnit):
    dr=np.sqrt(dR(pos, core).x**2+dR(pos, core).y**2+dR(pos, core).z**2)
    a=abs(dR(pos, core).x*dirUnit.x+dR(pos, core).y*dirUnit.y+dR(pos, core).z*dirUnit.z)

    return np.sqrt(dr**2-a**2)

def TimeDelay(pos, time_core, core , dirUnit):
    c=I3Constants.c
    t_plane=time_core-(1/c)*(dR(pos, core).x*dirUnit.x+dR(pos, core).y*dirUnit.y+dR(pos, core).z*dirUnit.z)

    return t_plane


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
        lates = []

        xs = []
        ys = []
        xsall = []
        ysall = []
        xsinsc = []
        ysinsc = []
        


        pulseSeriesMap = frame["SiPMRecoPulses"]
        for scintkey in self.i3scintgeo.keys():
            pos = self.i3scintgeo[scintkey].position
            regularx = pos.x
            regulary = pos.y

            xsall.append(regularx)
            ysall.append(regulary)


        # Iterate over all the panels in the frame
        for scintkey in pulseSeriesMap.keys():
            pulse = pulseSeriesMap[scintkey][0]

            # Extract the signal properties for this panel
            signalAmplitude = pulse.charge
            signalArrivalTime = pulse.time

            # 3D vector descibing the location of the scintillator
            pos = self.i3scintgeo[scintkey].position

            # 3D vector describing the location of the core
            core = particle.pos

            # Unit vector describing the direction of shower propagation
            dirUnit = particle.dir

            # Time at which the shower core gets to the ground
            time_core = particle.time    

            #Position in shower coordinates
            posinsc = radcube.GetShowerFromIC(pos - particle.pos, particle.dir)



            ########################################
            # You will have to fill this part out to calculate the distance from the shower axis and the
            ########################################
            radius = Radius(pos, core, dirUnit)
            delay = signalArrivalTime-TimeDelay(pos, time_core, core , -dirUnit)
            xinsc = posinsc.x
            yinsc = posinsc.y
            late = (posinsc/abs(posinsc)).y

            regularx = pos.x
            regulary = pos.y

            # Store things for plotting later
            amps.append(signalAmplitude)
            times.append(signalArrivalTime)
            radii.append(radius)
            delays.append(delay)
            lates.append(late) #These are ys in sc
            xsinsc.append(xinsc)
            ysinsc.append(yinsc)
            xs.append(regularx)
            ys.append(regulary)

        # Convert into numpy arrays to make things easier
        amps = np.array(amps)
        times = np.array(times)
        radii = np.array(radii)
        delays = np.array(delays)
        lates = np.array(lates) #ys in sc
        xs = np.array(xs)
        ys  = np.array(ys)
        xsall = np.array(xsall)
        ysall  = np.array(ysall)
        xsinsc = np.array(xsinsc)
        ysinsc = np.array(ysinsc)
       

        # Start making the plots
        NRows = 2
        NCols = 2
        gs = gridspec.GridSpec(NRows, NCols, wspace=0.3, hspace=0.3)
        fig = plt.figure(figsize=(6 * NCols, 5 * NRows))

        # LDF Calculation
        ax = fig.add_subplot(gs[0])
        ax.scatter(radii / I3Units.m, amps, c=lates, cmap='coolwarm' )
        ax.set_yscale("log")
        ax.set_xlabel("Radius [m]")
        ax.set_ylabel("Signal [MIP]")

        # Shower front calculation
        ax = fig.add_subplot(gs[1])
        ax.scatter(radii / I3Units.m, delays / I3Units.ns, c=lates, cmap='coolwarm')
        ax.set_xlabel("Radius [m]")
        ax.set_ylabel("Delay w.r.t. Shower Plane [ns]")

        # Shower front calculation
        ax = fig.add_subplot(gs[2])
        ax.scatter(xsall / I3Units.m, ysall / I3Units.m, alpha=0.2, color="k")
        ax.scatter(xs / I3Units.m, ys / I3Units.m, c=lates, cmap='coolwarm')
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_aspect("equal")

        # Shower front calculation
        ax = fig.add_subplot(gs[3])
        ax.scatter(xsinsc / I3Units.m, ysinsc / I3Units.m, c=lates, cmap='coolwarm')
        ax.set_xlabel("x in s.c [m]")
        ax.set_ylabel("y in s.c. [m]")
        ax.set_aspect("equal")



        print("Saving plot as", args.output)
        fig.savefig(args.output, bbox_inches="tight")
        exit()


tray = I3Tray()
tray.AddModule("I3Reader", "Reader", FilenameList=args.input)
tray.AddModule(PlotLDF, "Plotter")
tray.Execute()
tray.Finish()
