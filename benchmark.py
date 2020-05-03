#! /usr/bin/python3
import csv
import subprocess
import random

runs = 5
repeats = 5
mpi = True
omp = False
perf = False
doAverage = False
if (omp and mpi): print("DO NOT DO THIS")

#TEXT THAT WILL BE DESPLAYED AT THE TOP
#AND PARAMATERS THAT ARE BEING TESTED
paramaters = 4 #0: startup 1: jacobi 2: unknown 3: elapsed
topText = "graph size, steps, "
if mpi: topText += "divides, processes, "
if omp: topText += "threads, "
topText += ("time startup, per startup, " +
    "time step, per step, "
    "time unknown, per unknown, ")
if mpi: #comm if  MPI
  topText += ("time local com, per local com, time global com, per global com, ")
  paramaters += 2
topText += ("time elapsed, per elapsed, " +
    "cache references, cache misses, cpi\n")

#INCLUDE THE PERF RESULTS
perfstr = []
perfcount = 0
if perf:
  perfParamaters = ["cache-references", "cache-misses", "cycles"]
  perfstr = ["perf",  "stat", "-e", ",".join(perfParamaters)]
  perfcount = len(perfParamaters) #make sure to equal the number of things in perf

#PREPARE THE CSV FILE
f2 = open("Log/benchmark.csv", "w")
f2.write(topText)

#PRESET VALUES TO CHANGE OVER RUNS
nps = [(i)%8+1  for i in range(runs)]
ts  = [(i)%8+1 for i in range(runs)]
gs  = [500 for i in range(runs)]
ss  = [100 for i in range(runs)]

#OVERALL RUNS
for r in range(runs):
  print("range ", r)
  g = gs[r]   #grid size
  s = ss[r]    #steps
  np = nps[r] #number of processes MPI ONLY
  d = 1       #divides MPI ONLY
  t = ts[r]   #threads OMP ONLY
  if doAverage:
    times = [0 for x in range(paramaters)]
    perfRes = [0 for x in range(perfcount)]
    percents = [0 for x in range(paramaters)]
  else:
    times = [[] for x in range(paramaters)]
    perfRes = [[] for x in range(perfcount)]
    percents = [[] for x in range(paramaters)]

  #DONE FOR AVERAGING
  for i in range(repeats):

    #DIFFERENT COMMANDS
    if mpi:
      command = ["mpirun", "-np", str(np), "./grayscott-mpi", "-g", str(g), "-r", "1", "-s", str(s), "-d", str(d), "-I"]
    elif omp:
      command = ["./grayscott-omp", "-g", str(g), "-r", "1", "-s", str(s), "-t", str(t), "-I"]
    else:
      command = ["./grayscott-seq" , "-g", str(g), "-r", "1", "-s", str(s), "-I"]
    print("iteration", i, " ".join(perfstr + command))

    #REDIRECT STDOUT AND STDIN
    logfile = "Log/log" + str(r) + ".txt"
    ou = open(logfile, "wb") #for instrumentation
    er = open(logfile, "a")  #for perf
    subprocess.run((perfstr + command), stdout=ou, stderr=er)

    #BEGIN READING FROM TXT
    f1 = open(logfile, "r")
    for p in range(paramaters):
      split = f1.readline().split()
      if doAverage:
        times[p] += float(split[0])
        percents[p] += float(split[2])
      else:
        times[p].append(float(split[0]))
        percents[p].append(float(split[2]))

      

    if perf:
      [f1.readline() for IGNORE_LINES in range(3)] #ignore black perf response
      for p in range(perfcount):
        split = f1.readline().split()
        if doAverage:
          perfRes[p] += float(split[0].replace(",","")) #replace 999,999 to 999999
        else:
          perfRes[p].append(float(split[0].replace(",","")))
        


  #WRITE THE SIMULATION DESCRIPTION
  description = str(g) + ", " + str(s) 
  if mpi: description += ", " + str(d) + ", " + str(np)
  if omp: description += ", " + str(t)

  if doAverage:
    f2.write(description)

    #WRITE OUT THE AVERAGE VALUES
    for p in range(paramaters):
      time = times[p]/repeats
      percent = percents[p]/repeats
      st = ", " + str(time) + ", " + str(percent)
      f2.write(st)
    if perf:
      for p in range(perfcount):
        ref = perfRes[p]/repeats
        st = ", " + str(ref)
        f2.write (st)
    f2.write("\n") #
  else:
    for t in range(repeats):
      f2.write(description)

      for p in range(paramaters):
        time = times[p][t]
        percent = percents[p][t]
        st = ", " + str(time) + ", " + str(percent)
        f2.write(st)
      if perf:
        for p in range(perfcount):
          ref = perfRes[p][t]
          st = ", " + str(ref)
          f2.write (st)
      f2.write("\n")
          
    