import sys, os, glob
import h5py as h5
import numpy

# Initialize the reader by passing the directory containing the CGNS files. This
# returns a list containing sorted file pathes available for reading.
def init(basedir):
  global base
  base = basedir
  files = glob.glob(base + "/part-*.cgns")
  if(len(files) == 0):
    print("cannot find any part-*.cgns files in", base)
    sys.exit()
  else:
    times = [i[i.find("part-")+5:-5] for i in files]
  return sorted(times, key=float)

# Open a particular CGNS file using a time value in the list returned by init().
def open(time):
  infile = base + "/part-" + time + ".cgns"
  try:
    f = h5.File(infile, 'r')
    return f
  except OSError:
    f = None
    print("file", infile, "does not exist")
    return f

def close(f):
  f.close()

# Read the time.
def read_time(f):
  t1 = f["/Base/Zone0/Etc/Time/ data"][0]
  return t1

# Read the particle positions.
def read_part_position(f):
  x1 = numpy.array(f["/Base/Zone0/GridCoordinates/CoordinateX/ data"])
  y1 = numpy.array(f["/Base/Zone0/GridCoordinates/CoordinateY/ data"])
  z1 = numpy.array(f["/Base/Zone0/GridCoordinates/CoordinateZ/ data"])
  return (x1,y1,z1)

# Read the particle velocities.
def read_part_velocity(f):
  u1 = numpy.array(f["/Base/Zone0/Solution/VelocityX/ data"])
  v1 = numpy.array(f["/Base/Zone0/Solution/VelocityY/ data"])
  w1 = numpy.array(f["/Base/Zone0/Solution/VelocityZ/ data"])
  return (u1,v1,w1)

# do periodic counting
def periodic_crossings(x1, y1, z1, x0, y0, z0, Lx, Ly, Lz, a):
  np = len(x1)
  dx = x1 - x0
  ret_x = numpy.zeros(np)
  ret_x[dx > (Lx-a)] = 1
  ret_x[dx < -(Lx-a)] = -1
  
  dy = y1 - y0
  ret_y = numpy.zeros(np)
  ret_y[dy > (Ly-a)] = 1
  ret_y[dy < -(Ly-a)] = -1

  dz = z1 - z0
  ret_z = numpy.zeros(np)
  ret_z[dz > (Lz-a)] = 1
  ret_z[dz < -(Lz-a)] = -1

  return (ret_x, ret_y, ret_z)

def cor(u, v, w, u0, v0, w0):
  ret_u = numpy.multiply(u0,u)
  ret_v = numpy.multiply(v0,v)
  ret_w = numpy.multiply(w0,w)
  return (ret_u, ret_v, ret_w)
