# README
#
# Dependencies:
#  - ParaView-4.2.0
#  - ffmpeg release/2.8 (git://source.ffmpeg.org/ffmpeg.git)
#  - mesa
#
# Usage:
#  - Use pvpython to run this script
#  - Use ParaView to set up the scene as you like
#  - Make sure to rename the filter containing the data as 'flow' and 'part'
#  - Save the ParaView state file (File->Save State)
#  - Create a directory for storing the image files
#  - Run animate-cgns.py and follow the instructions to generate the animation

import sys
from paraview.simple import *
import glob, os
import subprocess
import math
import re

def sorted_nicely(path):
    flow_files = glob.glob(path + "/output/flow-*.cgns")
    part_files = glob.glob(path + "/output/part-*.cgns")
    flow_files.sort(key=lambda x:float(x[x.index("flow-")+5:x.index(".cgns")]))
    part_files.sort(key=lambda x:float(x[x.index("part-")+5:x.index(".cgns")]))
    times = [float(x[x.index("flow-")+5:x.index(".cgns")]) for x in flow_files]
    return (times, flow_files, part_files)

# outputs folder
root = "/home-4/yzhan175@jhu.edu/scratch/fluidized_bed/fluidization"
# Paraview state file
state = root + "/visualizer/" + "state.pvsm"
# time start and end
ts = 4990.
te = 5000.5
# image folder
img = "/home-4/yzhan175@jhu.edu/scratch/fluidized_bed/fluidization/visualizer/img"
# animation folder
anim = "/home-4/yzhan175@jhu.edu/scratch/fluidized_bed/fluidization/visualizer/anim"

# load state file
# make sure to build state file using the filter names 'flow' and 'part'
LoadState(state)
flow = FindSource("flow")
part = FindSource("part")
flow_filter = FindSource("flow_filter")
part_filter = FindSource("part_filter")
if flow == None or part == None or flow_filter == None or part_filter == None:
    print "Make sure to build the state file using right names"
view = GetRenderView()

#view.WriteImage(img + "/tmp.png", "vtkPNGWriter", 1)
#os.remove(img + "tmp.png")

if not os.path.exists(img):
    os.makedirs(img)
if not os.path.exists(anim):
    os.makedirs(anim)

(times, flow_files, part_files) = sorted_nicely(root)
mag = int(math.floor(math.log10(float(te))))

# go through all files
for i in range(len(times)):
    if times[i] >= ts and times[i] <= te:
        # change to file given by time
        flow.FileName = flow_files[i]
        flow_filter.FileNameChanged()
        part.FileName = part_files[i]
        part_filter.FileNameChanged()
        print("Saving image for t = " + str(times[i]))
        
        # pad image output time stamp for ffmpeg
        ztime = str(int(times[i])).zfill(mag + 5)
        
        # save screen shot
        #SaveScreenshot(img + "/img-" + ztime + ".png")
        view.WriteImage(img + "/img-" + ztime + ".png", "vtkPNGWriter", 1)

