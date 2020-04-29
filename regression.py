import subprocess
import random
for i in (range(100)):
  print("iteration", i)
  g = random.randint(0,1000)
  print("g", g)
  s = random.randint(5,1000)
  print("s", s)
  d = random.randint(1,4)
  print("d", d)
  np = random.randint(1,2) * d
  print("np", np)
  subprocess.call(["mpirun", "-np", str(np), "./grayscott-mpi", "-g", str(g), "-r", "1", "-s", str(s), "-d", str(d)])
  subprocess.call(["mv", "out/out0.txt", "out/RUN.txt"])
  subprocess.call(["mpirun", "-np", str(np), "./grayscott-mpi", "-g", str(g), "-r", "1", "-s", str(s), "-d", str(d)])
  subprocess.call(["mv", "out/out0.txt", "out/RUN1.txt"])
  subprocess.call(["diff", "out/RUN.txt", "out/RUN1.txt"])