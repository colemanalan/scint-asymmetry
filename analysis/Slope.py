#!/usr/bin/env python

## -------------------------------------------------------------------------------------------------
## Plot the rate between the signal and the average signal.
## -------------------------------------------------------------------------------------------------

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FixedLocator
from scipy.optimize import curve_fit
from scipy.stats import chisquare

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

#################### Function for fitting the Slope vs r data. A and B need to be fitted.
def fitting_slope(x, A, B, C):
    #y = -B*(np.expm1(-A*x))
    y = A*x/(np.sqrt(1+B*x**2))
    #y =  A*x**3+B*x**2+C*x
    return y

def Bootstraping(Xnd,Ynd): #Considering n dimensional arrays
    # num_trials = len(Xnd)
    # num_trials = Xnd.size
    idxfortrials = np.arange(Xnd.size)
    #print(idxfortrials)

    X = Xnd.flatten()
    Y = Ynd.flatten()
    num_trials = len(X)

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
        self.ALL_ratio = []
        self.ALL_xsinsc = []
        self.ALL_lates = [] # This is the sin psi. 
        self.ALL_slopeserr = []
        self.ALL_chisquare = []

    def Geometry(self, frame):
        self.i3scintgeo = frame["I3ScintGeometry"].scintgeo
        assert self.i3scintgeo
        self.narms = 0
        self.nrings = 0
        for scintkey in self.i3scintgeo.keys():
            self.narms = max(self.narms, scintkey.panel)
            self.nrings = max(self.nrings, scintkey.station)
        self.ALL_radii = np.zeros(self.nrings)
        self.ALL_slopes = np.zeros(self.nrings)
        self.ALL_slopeserr = np.zeros(self.nrings)

        #print("Found scint geometry for", len(self.i3scintgeo), "panels")
        #print("Number of arms ", self.narms)
        #print("Number of rings ", self.nrings)

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
        avgs = np.zeros(self.nrings)

        rates_per_event = np.zeros((self.narms, self.nrings))
        # rate_per_ring = np.zeros(self.narms)
        # rates_per_event = []
        
       
        pulseSeriesMap = frame["SiPMRecoPulses"]

        for scintkey in self.i3scintgeo.keys():
            pos = self.i3scintgeo[scintkey].position


            if scintkey.panel == 1:
                # 3D vector descibing the location of the scintillator
                # 3D vector describing the location of the core
                core = particle.pos

                # Unit vector describing the direction of shower propagation
                dirUnit = particle.dir
                
                radius = Radius(pos, core, dirUnit)

                self.ALL_radii[scintkey.station - 1] = radius

            #Position in shower coordinates
            posinsc = radcube.GetShowerFromIC(pos - particle.pos, particle.dir)

            yinsc = posinsc.y
            late = (posinsc/np.sqrt(posinsc.x**2+posinsc.y**2)).y

            # print(scintkey, posinsc, yinsc, late)

            lates[scintkey.panel-1, scintkey.station -1] = late #These are ys in sc

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

            regularx = pos.x #In IceCube coordinates
            regulary =pos.y #In IceCube coordinates

            # Store things for plotting later, panel 
            amps[scintkey.panel-1, scintkey.station -1] = signalAmplitude
            #times[scintkey.panel-1, scintkey.station -1] = signalArrivalTime
            delays[scintkey.panel-1, scintkey.station -1] = delay
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
        #This runs after the last shower is read.
        self.ALL_ratio = np.array(self.ALL_ratio)
        self.ALL_lates = np.array(self.ALL_lates)
        self.ALL_slopes = np.array(self.ALL_slopes)
        self.ALL_slopeserr = np.array(self.ALL_slopeserr)

        #Chi square considering 1. all data 2. the mean for each value of sin(psi)
        chisquares = []
        chisquares_asso = []

        #Mean for each value of sin(psi) and value of sin(psi) considering only one element per value.
        means2 = np.zeros(7)
        lates_means2 = np.zeros(7)

        ############################ PLOT ############################
        #Configuation for kind of printing

        presentation =  0

        Ini = 50 #Just an initialization

        if presentation == 0:
            NRows = 5
            NCols = 6

        elif presentation == 1:
            NRows = 2
            NCols = 3

        else:
            NRows = 1
            NCols = 1


        gs = gridspec.GridSpec(NRows, NCols, wspace=0.3, hspace=0.3)
        fig = plt.figure(figsize=(6 * NCols, 5 * NRows))
        #fig.suptitle(r"Signal rate in terms of $\sin\psi$ for $E=10^{16.5}$ eV, $\theta=40\degree$")        
        plt.rc('axes', titlesize=13)


      
        for ii in np.arange(self.nrings):
            r_plot = 15*(ii+1)


            ######### SELECTION FOR THE RINGS AND EVENTS IN THE PLOTS#########
            # We only consider rings with at most 3 non-heatted scintillators. 
            sel = np.sum(self.ALL_ratio[:,:,ii] != 0, axis = 1) > 9
            # print("selection ", sel)

            ratio_to_use = self.ALL_ratio[:,:,ii][sel]
            lates_to_use = self.ALL_lates[:,:,ii][sel]

            if  np.any(ratio_to_use):

                xfittings = np.linspace(-1, 1, 50)

                # print("Lates to use", lates_to_use.flatten())
                #This is for the nominal slope for the data
                coef=Linearfit(lates_to_use.flatten(), ratio_to_use.flatten())
                self.ALL_slopes[ii] = coef[0]
                fit = np.poly1d(coef)

                bootstraping1 = Bootstraping(lates_to_use , ratio_to_use)


                # for param in bootstraping1:
                BOOT = np.array([np.poly1d(param)(xfittings) for param in bootstraping1])
                self.ALL_slopeserr[ii] = np.std(bootstraping1[:,0])

                #Means and standard deviations from the bootstraping. These are not the means and stds from the actual data obtained from the events.
                means_boot = np.average(BOOT, axis = 0)
                stds_boot = np.std(BOOT, axis = 0)
                # print("means_boot", means_boot, "stds_boot", stds_boot)

                
                ##Compute the chi square value for this plot
                #Expected values, i.e, the value of the fitting at each arm 
                #ratio_expected = fit(lates_to_use)
                #print("Ratio expected", ratio_to_use)

                #chi_sqr = chisquare(ratio_to_use, ratio_expected, axis = None)
                #chisquares.append(chi_sqr[0])


                #Means and stds for the actual data obtained from the events --> without considering associated arms depending on lateness
                means = []
                for jj in np.arange(self.narms):
                    means.append( np.mean(ratio_to_use[:,jj]) )
                means = np.array(means)

                # ----> Association of arms
                
                for jj in np.arange(7):
                    if jj<3:
                        asso = np.concatenate((ratio_to_use[:,jj],ratio_to_use[ :,(6-jj) ]))
                        means2[jj]=np.mean( asso )
                        lates_means2[jj]= lates_to_use[0,jj]

                    elif jj ==3:
                        means2[jj]=np.mean(ratio_to_use[:,jj])
                        lates_means2[jj]= lates_to_use[0,jj]
                    elif jj ==6:
                        means2[jj]=np.mean(ratio_to_use[:,9])
                        lates_means2[jj]= lates_to_use[0,9]
                    else:
                        asso = np.concatenate((ratio_to_use[:, (3+jj)],ratio_to_use[ :,(15-jj) ]))
                        means2[jj]=np.mean( asso )
                        lates_means2[jj]= lates_to_use[0,3+jj]


                ##Chi square value considering mean 2
                chi_sqr_asso = chisquare(means2, fit(lates_means2), axis = None)

                chisquares_asso.append(chi_sqr_asso[0])
                self.ALL_chisquare.append(chi_sqr_asso)

                Ini = 50
                ############################ PLOT ############################

                if presentation == 0:
                    Ini = ii

                elif presentation == 1:

                    if r_plot == 60:
                        Ini = 0
                    elif r_plot == 120:
                        Ini = 1
                    elif r_plot == 180:
                        Ini = 2
                    #elif r_plot == 300:
                     #   Ini = 3
                    #print("Hello", r_plot, Ini, NRows, NCols)



                ####PLOT
                if Ini != 50: #Different that the initialization

                    
                    #Definition of the subplot
                    ax = fig.add_subplot(gs[Ini])

                    #Bootstraping
                    ax.fill_between(xfittings, means_boot-stds_boot, means_boot+stds_boot, color = "k", alpha= 0.4)

                    #Linear fit for the data
                    ax.plot(xfittings, fit(xfittings), color='r')

                    #Data
                    ax.scatter(lates_to_use, ratio_to_use, c=lates_to_use, cmap='coolwarm' )
                    
                    #Mean of data
                    ax.scatter(lates_to_use[0,:], means, marker="x")
                    ax.scatter(lates_means2, means2, marker="d", color="seagreen")

                    #ax.set_yscale("log")
                    plt.rc('axes', labelsize=13)
                    #ax.set_ylim(-0.2,5)

                    ax.set_xlabel("$\\sin\\psi$")
                    ax.set_ylabel("$S/\\bar{S}$")
                    ax.set_title(r"r= {0} m, $\chi^2$={1}, slope={2}.".format( r_plot, np.round(chi_sqr_asso[0],3), np.round(self.ALL_slopes[ii],3) ))      
                
        self.ALL_chisquare = np.array(self.ALL_chisquare)

        ##Selection for the creation of the subarrays
        sel_slopes = self.ALL_slopes != 0
        sel_chi_sqr = self.ALL_chisquare[:,0] < 0.5 
        #print(sel_chi_sqr)

        plot_slopes = self.ALL_slopes[sel_slopes][sel_chi_sqr]
        plot_err = self.ALL_slopeserr[sel_slopes][sel_chi_sqr]
        plot_radii = self.ALL_radii[sel_slopes][sel_chi_sqr]


        # print(plot_slopes)
        # print(plot_err)
        #print(plot_radii)

        plot_rmax = np.amax(plot_radii)
        rfittings = np. linspace(0, plot_rmax, 200)


        ########### Fitting the slope vs radius data 

        #Uncomment print("Making the fitting")

        param_fit_slope, covariance = curve_fit(fitting_slope, plot_radii, plot_slopes, maxfev = 20000)

        ## Fit for B(1-exp(-Ar))
        fit_A = param_fit_slope[0]
        fit_B = param_fit_slope[1]
        fit_C = param_fit_slope[2]

        fit_slope = fitting_slope(rfittings, fit_A, fit_B, fit_C)

        print(fit_A,";", fit_B,";", fit_C)

        #Uncomment print("Done")

        #Plot of amplitudes vs radius

        if presentation == 0 or presentation ==2:
        
            ax = fig.add_subplot(gs[NCols*NRows-1])

            ax.errorbar(plot_radii, plot_slopes, plot_err,  color="k")
            ax.tick_params(axis='both', which='both', direction='in')
            ax.yaxis.set_ticks_position('both')
            ax.xaxis.set_ticks_position('both')

            ax.scatter(plot_radii, plot_slopes, color="steelblue")
            ax.plot(rfittings,fit_slope, label = 'fit', color ='purple' )

            #ax.set_yscale("log")
            ax.set_xlabel("Radius [m]")
            ax.set_ylabel("Slope")
            ax.set_title(r"Slope in terms of the radius for $E=10^{16.0}$ eV, $\theta=30\degree$")

        #Uncomment print("Saving plot as", args.output)
        fig.savefig(args.output, bbox_inches="tight")

        #Uncomment print("Done!")


tray = I3Tray()
tray.AddModule("I3Reader", "Reader", FilenameList=args.input)
tray.AddModule(PlotLDF, "Plotter")
tray.Execute()
tray.Finish()