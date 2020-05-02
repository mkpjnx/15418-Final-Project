import csv
import subprocess
import random

omp = True
perf = True
topText = 0
args = 4
runs = 8
average = 5
topText = ("graph size, steps, threads, " +
    "time startup, per startup, " +
    "time jacobi, per jacobi, "
    "time unknown, per unknown, " +
    "time elapsed, per elapsed," +
    "cache refrences, cache misses, cpi\n")

if perf:
  command = ["perf",  "stat", "-e", "cache-references,cache-misses,cycles"]
  perfcount = 3
else:
  command = []
  perfcount = 0

f2 = open("Log/run.csv", "w")
f2.write(topText)

ts = [1,2,3,4,5,6,7,8]
for r in range(runs):
  print("range ", r)
  g = 256
  s = 1000
  t = ts[r]
  times = [0.0,0.0,0.0,0.0,0.0,0.0,0.0]
  percents = [0,0,0,0]
  for i in range(average):
    if omp:
      command += ["./grayscott-omp", "-g", str(g), "-r", "1", "-s", str(s), "-t", str(t), "-I"]
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
      print(line)
      spl  = line.split()
      val = float(spl[0].replace(",",""))
      times[args+f] += val

  description = str(g) + ", " + str(s) + ", "  + str(t)
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