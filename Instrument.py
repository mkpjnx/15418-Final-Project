import csv
import subprocess
import random

mpi = False
topText = 0
args = 4
topText = ("graph size, steps, divides, processes, " +
    "time startup, per startup, " +
    "time jacobi, per jacobi, ")

if mpi:
  topText += ("time local com, per local com, " +
    "time global com, per global com, ")
  args += 2

topText += ("time unknown, per unknown, " +
  "time elapsed, per elapsed\n")

f2 = open("Log/run.csv", "w")
f2.write(topText)

nps = [1,2,3,4,5,6,7,8]
runs = 5
for t in range(8):
  print("range ", t)
  g = 500 * t
  s = 100
  d = 1
  np = 1# nps[t]
  times = [0,0,0,0,0,0]
  percents = [0,0,0,0,0,0]
  for i in range(runs):
    command = ["./grayscott-seq" , "-g", str(g), "-r", "1", "-s", str(s), "-d", str(d), "-I"]
    if mpi:
      command = ["mpirun", "-np", str(np), "./grayscott-mpi", "-g", str(g), "-r", "1", "-s", str(s), "-d", str(d), "-I"]
    print("iteration", i)

    ou = open("Log/log1.txt", "wb")
    subprocess.run(command, stdout=ou)

    f1 = open("Log/log1.txt", "r")
    for a in range(args):
      line = f1.readline()
      spl = line.split()
      times[a] += int(spl[0])
      percents[a] += float(spl[2])

  description = str(g) + ", " + str(s) + ", " + str(d) + ", " + str(np)
  f2.write(description)
  for a in range(args):
    time = times[a]/runs
    percent = percents[a]/runs
    st = ", " + str(time) + ", " + str(percent)
    f2.write(st)
    
  f2.write("\n")