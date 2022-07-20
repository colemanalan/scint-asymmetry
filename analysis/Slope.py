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

import random

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

###Linear fit
def Linearfit(X,Y):
    ##Coefficients for the polynomial
    coef=np.polyfit(X,Y, 1)
    #print(coef)
    #Polynomials from the fitting
    #fit = np.poly1d(coef)
    return coef

def Bootstraping(X,Y):
    num_trials = len(X)
    idxfortrials = np.arange(X.size)
    #print(idxfortrials)

    #Definition of the vectors for whom we are going to take the linear fit
    Xnew = np.zeros(len(X))
    Ynew = np.zeros(len(X))
    fittings = []

    for ii in np.arange(num_trials):
        #Seed for random
        random.seed(ii)

        #Take a sample of k elements with repetition from the set of indexes
        R = random.choices(idxfortrials, k = X.size)
        #print("INDEXES FOR THE TRIALS")
        #print(R)

        for jj in idxfortrials:
            jjnew = R[jj]
            #print(jjnew)

            #Fill the new vectors with the corresponding components of the original ones.
            Xnew[jj] = X[jjnew]
            Ynew[jj] = Y[jjnew]
        # print("X new")
        # print(Xnew)
        # print("Y new")
        # print(Ynew)
        # print("Linear fit")
        # print(Linearfit(Xnew,Ynew))
        #Linear fit for the new vectors
        fittings.append(Linearfit(Xnew,Ynew))

    #print("FITTINGS OF THE BOOTSTRAPING:")
    #print(fittings)
    return np.array(fittings)



class PlotLDF(icetray.I3Module):
    def __init__(self, ctx):
        icetray.I3Module.__init__(self, ctx)

        self.i3scintgeo = False

        self.

    def Geometry(self, frame):
        self.i3scintgeo = frame["I3ScintGeometry"].scintgeo
        assert self.i3scintgeo
        self.narms = 0
        self.nrings = 0
        for scintkey in self.i3scintgeo.keys():
            self.narms = max(self.narms, scintkey.panel)
            self.nrings = max(self.nrings, scintkey.station)

        #Arrays with the data collected for all simulations
        self.ALL_radii = []
        self.ALL_amps = []
        self.ALL_xsinsc = []
        self.ALL_ysinsc = []
        self.ALL_slopes = []
        self.ALL_slopeserr = []


        print("Found scint geometry for", len(self.i3scintgeo), "panels")
        print(self.narms)
        print(self.nrings)

    def DAQ(self, frame):
        particle = frame["MCPrimary"]

        amps = np.zeros((self.narms, self.nrings))
        radii = np.zeros(self.nrings)
        times = np.zeros((self.narms, self.nrings))
        delays = np.zeros((self.narms, self.nrings))
        lates = np.zeros((self.narms, self.nrings))

        xs = np.zeros((self.narms, self.nrings))
        ys = np.zeros((self.narms, self.nrings))
        xsall = np.zeros((self.narms, self.nrings))
        ysall = np.zeros((self.narms, self.nrings))
        xsinsc = np.zeros((self.narms, self.nrings))
        ysinsc = np.zeros((self.narms, self.nrings))
        slopes = np.zeros(self.nrings)
        slopeserr = np.zeros(self.nrings)
        
       
        pulseSeriesMap = frame["SiPMRecoPulses"]

        for scintkey in self.i3scintgeo.keys():

            if scintkey.panel == 1:
                # 3D vector descibing the location of the scintillator
                pos = self.i3scintgeo[scintkey].position

                # 3D vector describing the location of the core
                core = particle.pos

                # Unit vector describing the direction of shower propagation
                dirUnit = particle.dir
                
                radius = Radius(pos, core, dirUnit)

                radii[scintkey.station - 1] = radius

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

            # Store things for plotting later
            amps[scintkey.panel-1, scintkey.station -1] = signalAmplitude
            #times[scintkey.panel-1, scintkey.station -1] = signalArrivalTime
            delays[scintkey.panel-1, scintkey.station -1] = delay
            lates[scintkey.panel-1, scintkey.station -1] = late #These are ys in sc
            xsinsc[scintkey.panel-1, scintkey.station -1] = xinsc
            ysinsc[scintkey.panel-1, scintkey.station -1] = yinsc
            xs[scintkey.panel-1, scintkey.station -1] = regularx
            ys[scintkey.panel-1, scintkey.station -1] = regulary



        # Variables for the fittings on the plots
        xfittings = np.linspace(-1, 1, 200)
        rfittings = np. linspace(0, self.nrings*15, 200)

    
        # Start making the plots and fill this in the for
        NRows = 6
        NCols = 6
        gs = gridspec.GridSpec(NRows, NCols, wspace=0.3, hspace=0.3)
        fig = plt.figure(figsize=(6 * NCols, 5 * NRows))
 
        #For doing each plot
        for ii in np.arange(self.nrings):

            lates_per_ring = lates[:,ii]

            #Average of the signal per ring.
            avg = np.average(amps[:,ii]) #THIS

            #Signals normalized to the average signal (in each ring).
            rate_per_ring = amps[:,ii]/avg #THIS

            #Coeficients for the linear fitting
            coef=Linearfit(lates_per_ring, rate_per_ring)
            slopes[ii] = coef[0]

            #Polynomials from the fitting
            fit = np.poly1d(coef)

            r_plot = 15*ii

            bootstraping1 = Bootstraping(lates_per_ring, rate_per_ring)

            #For plotting the fittings
            ax = fig.add_subplot(gs[ii])
            for pram in bootstraping1:
                plotting = np.poly1d(pram)
                ax.plot(xfittings, plotting(xfittings), color="k", alpha=0.1)

            ax.plot(xfittings, fit(xfittings), color='r')
            slopeserr[ii] = np.std(bootstraping1[:,0])/np.sqrt(len(bootstraping1))

            ax.scatter(lates_per_ring, rate_per_ring, c=lates_per_ring, cmap='coolwarm' )
            #ax.set_yscale("log")
            ax.set_xlabel("$\\sin\\psi$")
            ax.set_ylabel("$S/\\bar{S}$")
            ax.set_title(r"r= {0} m, $\bar{{S}} = ${1:0.1f}".format(r_plot, avg))

        slopes = np.array(slopes)

        #print(radii)
        print("SLOPES TO BE PLOTTED:")
        print(slopes)
        print(radii)
        print(slopeserr)

        ax = fig.add_subplot(gs[35])
        ax.errorbar(radii, slopes, slopeserr,  color="k")
        
        #ax.set_yscale("log")
        ax.set_xlabel("Radius [m]")
        ax.set_ylabel("Slope")
        ax.set_title("Slope in terms of the radius")

        print("Saving plot as", args.output)
        fig.savefig(args.output, bbox_inches="tight")

    def Finish(self):
        #This runs after the last shower is read in
        print("Done!")


tray = I3Tray()
tray.AddModule("I3Reader", "Reader", FilenameList=args.input)
tray.AddModule(PlotLDF, "Plotter")
tray.Execute()
tray.Finish()    