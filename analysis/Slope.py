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
        #Arrays with the data collected for all simulations
        self.ALL_radii = []
        self.ALL_ratio = []
        self.ALL_xsinsc = []
        self.ALL_lates = [] # This is the sin psi. 
        self.ALL_slopes = []
        self.ALL_slopeserr = []

    def Geometry(self, frame):
        self.i3scintgeo = frame["I3ScintGeometry"].scintgeo
        assert self.i3scintgeo
        self.narms = 0
        self.nrings = 0
        for scintkey in self.i3scintgeo.keys():
            self.narms = max(self.narms, scintkey.panel)
            self.nrings = max(self.nrings, scintkey.station)


        print("Found scint geometry for", len(self.i3scintgeo), "panels")
        # print("Number of arms ", self.narms)
        # print("Number of rings ", self.nrings)

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
        avgs = np.zeros(self.nrings)

        rates_per_event = np.zeros((self.narms, self.nrings))
        # rate_per_ring = np.zeros(self.narms)
        # rates_per_event = []
        
       
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

                ####################################################
                #Save the information about radii in the ALL vector.
                self.ALL_radii.append(radii)

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


            radius = Radius(pos, core, dirUnit)
            delay = signalArrivalTime-TimeDelay(pos, time_core, core , -dirUnit)
            #print(scintkey)
            xinsc = posinsc.x
            yinsc = posinsc.y
            late = (posinsc/abs(posinsc)).y

            regularx = pos.x #In IceCube coordinates
            regulary =pos.y #In IceCube coordinates

            # Store things for plotting later, panel 
            amps[scintkey.panel-1, scintkey.station -1] = signalAmplitude
            #times[scintkey.panel-1, scintkey.station -1] = signalArrivalTime
            delays[scintkey.panel-1, scintkey.station -1] = delay
            lates[scintkey.panel-1, scintkey.station -1] = late #These are ys in sc
            xsinsc[scintkey.panel-1, scintkey.station -1] = xinsc
            ysinsc[scintkey.panel-1, scintkey.station -1] = yinsc
            xs[scintkey.panel-1, scintkey.station -1] = regularx
            ys[scintkey.panel-1, scintkey.station -1] = regulary



        # print("Amps ", amps)
        # print("Lates", lates)
        # Variables for the fittings on the plots
        
        #For doing each plot
        for ii in np.arange(self.nrings):

            #Average of the signal per ring.
            avgs[ii] = np.average(amps[:,ii]) #THIS
            # print("average per ring", avgs[ii])

            #Signals normalized to the average signal (in each ring).
            if avgs[ii] != 0:
                rates_per_event [:, ii] = amps[:,ii]/avgs[ii]  
           
        
        
        self.ALL_ratio.append(rates_per_event)
        self.ALL_lates.append(lates)


    def Finish(self):
        #This runs after the last shower is read in

        self.ALL_ratio = np.array(self.ALL_ratio)
        self.ALL_lates = np.array(self.ALL_lates)

        xfittings = np.linspace(-1, 1, 200)
        rfittings = np. linspace(0, self.nrings*15, 200)

    
        # Start making the plots and fill this in the for
        NRows = 6
        NCols = 6
        gs = gridspec.GridSpec(NRows, NCols, wspace=0.3, hspace=0.3)
        fig = plt.figure(figsize=(6 * NCols, 5 * NRows))
 
        for ii in np.arange(self.nrings):

            r_plot = 15*ii

            #For plotting the fittings
            ax = fig.add_subplot(gs[ii])
            ax.scatter(self.ALL_lates[:,:,ii], self.ALL_ratio[:,:,ii], c=self.ALL_lates, cmap='coolwarm' )
            #ax.set_yscale("log")
            ax.set_xlabel("$\\sin\\psi$")
            ax.set_ylabel("$S/\\bar{S}$")
            ax.set_title(r"r= {0} m".format(r_plot))
        
    

        print("Saving plot as", args.output)
        fig.savefig(args.output, bbox_inches="tight")

        print("Done!")


tray = I3Tray()
tray.AddModule("I3Reader", "Reader", FilenameList=args.input)
tray.AddModule(PlotLDF, "Plotter")
tray.Execute()
tray.Finish()   