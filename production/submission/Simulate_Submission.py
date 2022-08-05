#!/bin/env python3

import glob
import os
import getpass
import subprocess

ABS_PATH_HERE = str(os.path.dirname(os.path.realpath(__file__)))
whoami = getpass.getuser()

simDir = "/home/ddazamarroquin/work/simulations/data/discrete" #This is the pathfor the showeres "I" submitted
i3Dir = ABS_PATH_HERE + "/../../data/new-i3-files"
logDir = ABS_PATH_HERE + "/logfiles"

prims = ["proton", "iron"]
energies = ["{0:0.1f}".format(15.5 + i / 10) for i in range(1)]
startid = 20

def MakeCondorSubmission(inputDir, outputDir, subfile):

    if not os.path.isdir("/scratch/{}".format(whoami)):
        print("Please ensure you are on the submit node and make the dir:")
        print("/scratch/{}".format(whoami))
        exit()

    file = open(subfile, "w")

    file.write("#!/bin/bash\n")
    file.write("Executable = {}/Simulate_RunControl.sh\n".format(ABS_PATH_HERE))
    file.write("Error = {0}/{1}_{2}_{3}_$(Process).err\n".format(logDir, prim, eng, zen))
    file.write("Output = {0}/{1}_{2}_{3}_$(Process).out\n".format(logDir, prim, eng, zen))
    file.write("Log = /scratch/{}/log.log\n".format(whoami))
    file.write("request_memory = 4GB\n")
    file.write("Arguments = {} {} $(Process) {}\n".format(inputDir, outputDir, startid))
    file.write("Queue 20")  # This is the number of CORSIKA files that will be processed (per eng/zen dir)
    file.close()


for prim in prims:
    directoryPrim = simDir + "/" + prim + "/"

    for eng in energies:
        directoryE = directoryPrim + "lgE_" + eng + "/"

        zens = list(glob.glob(directoryE + "Zen*"))
        zens.sort()
        zens = [zen.split("/")[-1] for zen in zens]
        zens = [zen.split("_")[-1] for zen in zens]

        for zen in zens:
            directoryZ = directoryE + "Zen_" + str(zen) + "/"
            outputDir = i3Dir + "/" + prim + "/" + "lgE_" + eng + "/" + "Zen_" + str(zen) + "/"
            MakeCondorSubmission(directoryZ, outputDir, "TempSub.sub")
            subprocess.call(["condor_submit", "TempSub.sub", "-batch-name", "{0}_{1}_{2}".format(prim, eng, zen)])

subprocess.call(["rm", "TempSub.sub"])
