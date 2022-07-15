#!/usr/bin/env python

## -------------------------------------------------------------------------------------------------
## Plot the rate between the signal and the average signal.
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
parser.add_argument("--output", type=str, default=ABS_PATH_HERE + "/../plots/SignalRate.pdf", help="Name for the pdf with the plotted data")
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
        slopes = []
        
       
        pulseSeriesMap = frame["SiPMRecoPulses"]

        #Number of rings
        nrings = 5
        narms = 12
        

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
            #print(scintkey)
            xinsc = posinsc.x
            yinsc = posinsc.y
            late = (posinsc/abs(posinsc)).y

            regularx = pos.x #In IceCube coordinates
            regulary =pos.y #In IceCube coordinates


            if (scintkey.station <= nrings):                
                

                # Store things for plotting later
                amps.append(signalAmplitude)
                #times.append(signalArrivalTime)
                delays.append(delay)
                lates.append(late) #These are ys in sc
                xsinsc.append(xinsc)
                ysinsc.append(yinsc)
                xs.append(regularx)
                ys.append(regulary)

                if (scintkey.panel ==1):
                    radii.append(radius)

        
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

        xfittings = np.linspace(-1, 1, 200)

    
        # Start making the plots and fill this in the for
        NRows = 2
        NCols = 3
        gs = gridspec.GridSpec(NRows, NCols, wspace=0.3, hspace=0.3)
        fig = plt.figure(figsize=(6 * NCols, 5 * NRows))
 
        #For doing each plot
        for ii in np.arange(nrings):

            lates_per_ring = lates[ii*narms:(ii+1)*narms]

            #Average of the signal per ring.
            avg = np.average(amps[ii*narms:(ii+1)*narms])

            #Signals normalized to the average signal (in each ring).
            rate_per_ring = amps[ii*narms:(ii+1)*narms]/avg

            #Coeficients for the linear fitting
            coef=np.polyfit(lates_per_ring, rate_per_ring, 1)
            slopes.append(coef[0])

            #Polynomials from the fitting
            fit1 = np.poly1d(coef)
            r_plot = 15*ii

            #For plotting the fittings
            ax = fig.add_subplot(gs[ii])
            ax.scatter(lates_per_ring, rate_per_ring, c=lates_per_ring, cmap='coolwarm' )
            ax.plot(xfittings, fit1(xfittings), color="k")
            #ax.set_yscale("log")
            ax.set_xlabel("$\\sin\\psi$")
            ax.set_ylabel("$S/\\bar{S}$")
            ax.set_title("r= "+str(r_plot)+" m, $\\bar{S}= $"+ str(avg)+".")
                
        slopes = np.array(slopes)

        print(radii)
        print(slopes)

        ax = fig.add_subplot(gs[5])
        ax.scatter(radii, slopes, color="steelblue")
        #ax.set_yscale("log")
        ax.set_xlabel("Radius [m]")
        ax.set_ylabel("Slope")
        ax.set_title("Slope in terms of the radius")

        print("Saving plot as", args.output)
        fig.savefig(args.output, bbox_inches="tight")
        exit()


tray = I3Tray()
tray.AddModule("I3Reader", "Reader", FilenameList=args.input)
tray.AddModule(PlotLDF, "Plotter")
tray.Execute()
tray.Finish()    