#!/usr/bin/env python

## -------------------------------------------------------------------------------------------------
## Script for running Scintillator simulation on a star-shaped geometry to study the second-order
## radial effects. This will only process one corsika file at a time run one shower at a time.
## -------------------------------------------------------------------------------------------------

from icecube.surface_sim_scripts.segments.segmentITExDefault import SimIceTopScint
from icecube.surface_sim_scripts import radiotools
from icecube.surface_sim_scripts.parser import StandardSurfaceArrayParser
from icecube.icetray.i3logging import *
from icecube import icetray, dataio, dataclasses, radcube
from I3Tray import I3Tray
from icecube.icetray import I3Units

import random
import numpy as np
import os

ABS_PATH_HERE = str(os.path.dirname(os.path.realpath(__file__)))
icetray.I3Logger.global_logger.set_level(icetray.I3LogLevel.LOG_INFO)


def GetRadius(particle, pos):
    # Particle is the primary particle and pos is the detector position
    x_c = particle.pos.x
    y_c = particle.pos.y
    z_c = particle.pos.z

    nx = particle.dir.x
    ny = particle.dir.y

    abs_x_sq = (pos[0] - x_c) ** 2 + (pos[1] - y_c) ** 2 + (pos[2] - z_c) ** 2

    n_prod_x = nx * (pos[0] - x_c) + ny * (pos[1] - y_c) - np.sqrt(1.0 - nx ** 2 - ny ** 2) * (pos[2] - z_c)

    return np.sqrt(abs_x_sq - n_prod_x ** 2)

def MakeI3Geometry(simfile, narms, dr, maxR):
    baseDir = os.path.dirname(simfile)
    showerID = os.path.basename(baseDir)

    print("simdir:", simfile)
    print("ID: ", showerID)

    inpFile = baseDir + "/SIM" + showerID + ".inp"

    with open(inpFile) as f:
        for line in f:
            if "THETAP" in line:
                zenith = float(line.split()[1]) * I3Units.degree
            elif "PHIP" in line:
                azimuth = float(line.split()[1]) * I3Units.degree
            elif "OBSLEV" in line:
                obslev = float(line.split()[1]) * I3Units.cm

    # Convert to I3
    azimuth -= np.pi  # Points to where it came from
    obslev -= dataclasses.I3Constants.OriginElev
    direction = dataclasses.I3Direction(zenith, azimuth)
    direction.rotate_z(radcube.GetMagneticRotation())  # Corsika in mag coords

    print("InI3Coords direction: {} obslev: {}".format(direction, obslev / I3Units.m))

    particle = dataclasses.I3Particle()
    particle.dir = dataclasses.I3Direction(zenith, azimuth)
    particle.pos = dataclasses.I3Position(0, 0, 0)

    def MakeScintGeo(position, orientation=0.0, heightAboveSnow=0.5 * I3Units.m, name="scintillator_test"):
        scintGeo = dataclasses.I3ScintGeo()
        scintGeo.position = position
        up = dataclasses.I3Direction(0, 0, 1)
        d_dir = dataclasses.I3Direction(-1, 0, 0)
        d_dir.rotate_z(orientation * I3Units.rad)
        scintGeo.orientation = dataclasses.I3Orientation(d_dir, up)
        scintGeo.heightAboveSnow = heightAboveSnow
        scintGeo.scintName = name
        scintGeo.scintType = dataclasses.I3ScintGeo.ScintType.MADKIT

        return scintGeo

    # Make the panels in the geometry
    radii = np.linspace(dr, maxR, round((maxR - dr) / dr) + 1)
    arms = np.linspace(0, np.pi * 2, narms + 1)[:-1]

    scintGeometry = dataclasses.I3Geometry()

    for ir, r in enumerate(radii):
        for iang, angle in enumerate(arms):
            scintkey = dataclasses.ScintKey(ir + 1, iang + 1)

            posInShCS = dataclasses.I3Position(r * np.cos(angle), r * np.sin(angle), 0)
            posInIC = radcube.GetICFromShower(posInShCS, particle.dir)
            posInIC += dataclasses.I3Direction(zenith, azimuth) * (posInIC.z - 10. * I3Units.cm) / np.cos(particle.dir.zenith)  # Project the panel a bit below obslev
            posInIC.z += obslev

            scintGeo = MakeScintGeo(posInIC, np.random.rand() * 2 * np.pi, 2.0 * I3Units.m, "arm{0}_rad{1:0.1f}".format(iang + 1, r))
            scintGeometry.scintgeo[scintkey] = scintGeo

    return scintGeometry


def SimulateSingleEvent(corsikaFile, output, options):
    tray = I3Tray()

    tray.AddSegment(SimIceTopScint, "SimIceTopScint", input=[corsikaFile], **options)

    tray.AddModule(
        "Delete",
        "Delete",
        Keys=["HitBinWidth", "BeaconLaunches", "MCTimeIncEventID", "RNGState", "IceTopCalibratedWaveformRange", "InIceRawData", "ZeroPaddedMap"],
        KeyStarts=["CalibratedIceTop", "IceTopComponentWaveforms"],
    )

    ##################Writing output of simulations####################
    tray.AddModule("I3Writer", "i3-writer", Filename=output, streams=[icetray.I3Frame.TrayInfo, icetray.I3Frame.Geometry, icetray.I3Frame.DAQ, icetray.I3Frame.Physics])

    tray.AddModule("TrashCan", "trashcan")
    tray.Execute()
    tray.Finish()


infiles, output, options = StandardSurfaceArrayParser()

if len(infiles) != 1:
    print("Error: This script is only designed to process a single CORSIKA file at a time")
    print("You have this script", len(infiles), "files")
    exit()

for corsikaFile in infiles:

    # Make a file with a Geometry Frame
    parentOutputDir = os.path.dirname(output)
    gFileName = parentOutputDir + "/temp_gcd_{}.i3.gz".format(output.split("_")[-1])
    print("Making GCD", gFileName)

    i3File = dataio.I3File(gFileName, "w")
    i3Geometry = MakeI3Geometry(corsikaFile, narms=12, dr=20 * I3Units.m, maxR=500 * I3Units.m)
    gFrame = icetray.I3Frame(icetray.I3Frame.Geometry)
    gFrame["I3ScintGeometry"] = i3Geometry
    gFrame["I3Geometry"] = dataclasses.I3Geometry()
    i3File.push(gFrame)

    cFrame = icetray.I3Frame(icetray.I3Frame.Calibration)
    cFrame["I3Calibration"] = dataclasses.I3Calibration()
    i3File.push(cFrame)

    dFrame = icetray.I3Frame(icetray.I3Frame.DetectorStatus)
    dFrame["I3DetectorStatus"] = dataclasses.I3DetectorStatus()
    i3File.push(dFrame)

    i3File.close()

    options["gcd"] = gFileName
    options["r"] = 0  # For these simulations, always fix the core at the center
    options["x"] = 0
    options["y"] = 0
    options["samples"] = 1
    options["no_ice_top"] = True

    random.seed(corsikaFile)
    options["seed"] = random.randint(1, 10000000)

    SimulateSingleEvent(corsikaFile, output, options)
    log_info(str(corsikaFile) + "-->" + str(output))

    os.remove(gFileName)
