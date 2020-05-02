#! /usr/bin/python3
import csv
import subprocess
import random

mpi = True
perf = True
topText = 0
args = 4
runs = 1
average = 1
topText = ("graph size, steps, divides, processes, " +
    "time startup, per startup, " +
    "time jacobi, per jacobi, "
    "time unknown, per unknown, ")

if mpi:
  topText += ("time local com, per local com, " +
    "time global com, per global com, ")
  args += 2

topText += ("time elapsed, per elapsed, " +
    "cache references, cache misses, cpi\n")


f2 = open("Log/run.csv", "w")
f2.write(topText)

nps = [1,2,3,4,5,6,7,8]
for r in range(runs):
  print("range ", r)
  g = 256
  s = 1000
  np = nps[r]
  d = 1
  times = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
  percents = [0,0,0,0,0,0,0,0,0,0,0]
  for i in range(average):
    if perf:
      command = ["perf",  "stat", "-e", "cache-references,cache-misses,cycles"]
      perfcount = 3
    else:
      command = []
      perfcount = 0
    if mpi:
      command += ["mpirun", "-np", str(np), "./grayscott-mpi", "-g", str(g), "-r", "1", "-s", str(s), "-d", str(d), "-I"]
    else:
      command += ["./grayscott-seq" , "-g", str(g), "-r", "1", "-s", str(s), "-I"]
    print("iteration", i)

    ou = open("Log/log1.txt", "wb")
    er = open("Log/log1.txt", "a")
    subprocess.run(command, stdout=ou, stderr=er)

    f1 = open("Log/log1.txt", "r")
    for a in range(args):
      line = f1.readline()
      spl = line.split()
      times[a] += float(spl[0])
      percents[a] += float(spl[2])

    for f in range(3):
      f1.readline()
    for f in range(perfcount):
      line = f1.readline()
      spl  = line.split()
      val = float(spl[0].replace(",",""))
      times[args+f] += val

  description = str(g) + ", " + str(s) + ", "  + str(d) + ", " + str(np)
  f2.write(description)
  for a in range(args):
    time = times[a]/average
    percent = percents[a]/average
    st = ", " + str(time) + ", " + str(percent)
    f2.write(st)
  for p in range(perfcount):
    ref = times[args + p]/average
    st = ", " + str(ref)
    f2.write (st)
    
  f2.write("\n")

topText = 0
args = 4
topText = ("graph size, steps, divides, processes, " +
    "time startup, per startup, " +
    "time jacobi, per jacobi, ")