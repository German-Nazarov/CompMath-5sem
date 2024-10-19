import numpy as np
import matplotlib.pyplot as plt

f = open(f'lab-3-spline_data.txt', 'r')

Ns = []
errs = []

for line in f:
    Ns.append(float((line.split("\t"))[0]))
    errs.append(float((line.split("\t"))[-1].replace('\n', '')))

Ns = np.log(np.array(Ns))
errs = np.log(np.array(errs))

plt.figure(figsize=(8, 6), dpi=100)

plt.scatter(Ns, errs, color = 'royalblue')

plt.plot(np.array([0, 6]), -2.85 * np.array([0, 6]) + 14, color = 'darkorange', ls='--', label = '$y = -2.85 \cdot x + 14$')

plt.legend()
plt.xlim([0, 6])
plt.ylim([-2, 10])
plt.xlabel("$\ln(N)$")
plt.ylabel("$\ln(Error)$")
plt.grid(True, linestyle="--")
plt.savefig("lab-3-spline-graph.png")
plt.show()