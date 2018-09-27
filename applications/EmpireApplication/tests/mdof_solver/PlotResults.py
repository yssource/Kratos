import matplotlib.pyplot as plt

X, Y = [], []
for line in open('results_mdof_sdof.dat', 'r'):
  values = [float(s) for s in line.split()]
  X.append(values[0])
  Y.append(values[1])

plt.plot(X, Y)
plt.show()