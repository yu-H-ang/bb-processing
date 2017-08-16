import subprocess

img = "/home-4/yzhan175@jhu.edu/scratch/fluidized_bed/fluidization/visualizer/img"
anim = "/home-4/yzhan175@jhu.edu/scratch/fluidized_bed/fluidization/visualizer/anim/output.mp4"
fps = 60

# stitch together using ffmpeg
subprocess.call(["ffmpeg",
    "-framerate", str(fps),
    "-f", "image2",
    "-start_number", "338",
    "-i", img+"/img-%08d.png",
    anim])
