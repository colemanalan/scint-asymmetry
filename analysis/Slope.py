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

            #Number of rings
            nrings = 5
            narms = 12


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

    

        #Separate the info for the different rings
        lates1 = lates[0:narms]
        lates2 = lates[narms: 2*narms]
        lates3 = lates[2*narms:3*narms]
        lates4 = lates[3*narms:4*narms]
        lates5 = lates[4*narms:5*narms]

        #Signals normalized to the average signals in each ring  

        rate1 = amps[0:narms]/np.average(amps[0:narms])
        rate2 = amps[narms: 2*narms]/np.average(amps[narms: 2*narms])
        rate3 = amps[2*narms:3*narms]/np.average(amps[2*narms:3*narms])
        rate4 = amps[3*narms:4*narms]/np.average(amps[3*narms:4*narms])
        rate5 = amps[4*narms:5*narms]/np.average(amps[4*narms:5*narms])

        #Getting linear fitting coefficients
        coef1 = np.polyfit(lates1, rate1, 1)
        coef2 = np.polyfit(lates2, rate2, 1) 
        coef3 = np.polyfit(lates3, rate3, 1)
        coef4 = np.polyfit(lates4, rate4, 1)
        coef5 = np.polyfit(lates5, rate5, 1) 

        print(coef1, coef2, coef3,coef4, coef5)

        #Save the slopes in an array
        slopes.append(coef1[0])
        slopes.append(coef2[0])
        slopes.append(coef3[0])
        slopes.append(coef4[0])
        slopes.append(coef5[0])

        slopes = np.array(slopes)

        print(radii)
        print(slopes)

        #Ploynomials from the fitting

        fit1 = np.poly1d(coef1)
        fit2 = np.poly1d(coef2)
        fit3 = np.poly1d(coef3)
        fit4 = np.poly1d(coef4)
        fit5 = np.poly1d(coef5)

        print(fit1)

        #For plotting the fittings
        xfittings = np.linspace(-1, 1, 200)



        ##########################################################################################################
        ##########################################################################################################
        ##Computation for the slope vs the radius.
        ##########################################################################################################
        ##########################################################################################################




        # Start making the plots
        NRows = 2
        NCols = 3
        gs = gridspec.GridSpec(NRows, NCols, wspace=0.3, hspace=0.3)
        fig = plt.figure(figsize=(6 * NCols, 5 * NRows))

        # Signal normalized by the average signal
        ax = fig.add_subplot(gs[0])
        ax.scatter(lates1, rate1, c=lates1, cmap='coolwarm' )
        ax.plot(xfittings, fit1(xfittings), color="k")
        #ax.set_yscale("log")
        ax.set_xlabel("$\\sin\\psi$")
        ax.set_ylabel("$S/\\bar{S}$")
        ax.set_title("r = 15 m")

        ax = fig.add_subplot(gs[1])
        ax.scatter(lates2, rate2, c=lates2, cmap='coolwarm' )
        ax.plot(xfittings, fit2(xfittings), color="k")
        #ax.set_yscale("log")
        ax.set_xlabel("$\\sin\\psi$")
        ax.set_ylabel("$S/\\bar{S}$")
        ax.set_title("r = 30 m")


        ax = fig.add_subplot(gs[2])
        ax.scatter(lates3, rate3, c=lates3, cmap='coolwarm' )
        ax.plot(xfittings, fit3(xfittings), color="k")
        #ax.set_yscale("log")
        ax.set_xlabel("$\\sin\\psi$")
        ax.set_ylabel("$S/\\bar{S}$")
        ax.set_title("r = 45 m")


        ax = fig.add_subplot(gs[3])
        ax.scatter(lates4, rate4, c=lates4, cmap='coolwarm' )
        ax.plot(xfittings, fit4(xfittings), color="k")
        #ax.set_yscale("log")
        ax.set_xlabel("$\\sin\\psi$")
        ax.set_ylabel("$S/\\bar{S}$")
        ax.set_title("r = 60 m")


        ax = fig.add_subplot(gs[4])
        ax.scatter(lates5, rate5, c=lates5, cmap='coolwarm' )
        ax.plot(xfittings, fit5(xfittings), color="k")
        #ax.set_yscale("log")
        ax.set_xlabel("$\\sin\\psi$")
        ax.set_ylabel("$S/\\bar{S}$")
        ax.set_title("r = 75 m")

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