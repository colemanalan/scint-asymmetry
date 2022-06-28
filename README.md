# scint-asymmetry

Set of code to study the signal asymmetry within showers detected by the scintillator array.

### Shower Production

You must already have a set of CORSIKA showers to be used to study these effects.
Propagate the particles into the detector by runnning the script `production/SimulateShower.py`.
This is more easily done by submitting the jobs to the cluster using the script `production/submission/Simulate_Submission.py`.


### Analysis

#### LDF Plotting

To make a plot of the lateral distribution of signal for an individual air shower, run the script `analysis/PlotLateralDistribution.py`, supplying an I3File from the steps above. It should make one plot with the lateral distribution and the timing profile for that shower.
