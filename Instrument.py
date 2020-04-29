import csv
import subprocess
import random
f2 = open("Log/run.csv", "w")
f2.write("run, graph size, steps, divides, processes, ")
f2.write("time startup, per startup, ")
f2.write("time jacobi, per jacobi, ")
f2.write("time local com, per local com, ")
f2.write("time unknown, per unknown, ")
f2.write("time elapsed, per elapsed\n")

nps = [1,2,3,4,5,6,7,8]
runs = 5
for t in range(8):
  print("range ", t)
  g = 500
  s = 1000
  d = 1
  np =  nps[t]
  times = [0,0,0,0,0]
  percents = [0,0,0,0,0]
  for i in range(runs):
    print("iteration", i)


    ou = open("Log/log1.txt", "wb")
    subprocess.run(["mpirun", "-np", str(np), "./grayscott-mpi", "-g", str(g), "-r", "1", "-s", str(s), "-d", str(d), "-I"],
                  stdout=ou)

    f1 = open("Log/log1.txt", "r")
    for a in range(5):
      line = f1.readline()
      spl = line.split()
      times[a] += int(spl[0])
      percents[a] += float(spl[2])

  description = str(i) + ", " + str(g) + ", " + str(s) + ", " + str(d) + ", " + str(np)
  f2.write(description)
  for a in range(5):
    time = times[a]/runs
    percent = percents[a]/runs
    st = ", " + str(time) + ", " + str(percent)
    f2.write(st)
    
  f2.write("\n")