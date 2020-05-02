#! /usr/bin/python3
import csv
import subprocess
import random

mpi = True
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

nps = [12,24,36,48]
runs = 1
for t in range(3):
  print("range ", t)
  g = 1000
  s = 100
  d = 4
  np = nps[t]
  times = [0,0,0,0,0,0]
  percents = [0,0,0,0,0,0]
  for i in range(runs):
    command = ["./grayscott-seq" , "-g", str(g), "-r", "1", "-s", str(s), "-d", str(d), "-I"]
    if mpi:
      command = ["mpirun", "-np", str(np), "./grayscott-mpi", "-g", str(g), "-r", "1", "-s", str(s), "-d", str(d), "-I"]
    print("iteration", i)

    ou = open("Log/log1.txt", "wb")
    subprocess.call(command, stdout=ou)
'''
    f1 = open("Log/log1.txt", "r")
    for a in range(args):
      line = f1.readline()
      spl = line.split()
      times[a] += int(spl[0])
      percents[a] += float(spl[2])


  description = "%s,%s,%s,%s"%(g,s,d,np)
  f2.write(description)
  for a in range(args):
    time = times[a]/runs
    percent = percents[a]/runs
    st = ",%s,%s"%(time,percent)
    f2.write(st)
    
  f2.write("\n")
  '''