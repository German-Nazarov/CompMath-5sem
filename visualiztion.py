import numpy as np
import matplotlib.pyplot as plt

plt.figure(figsize=(8, 6), dpi=100)
for k in range(3, 6, 1):
    f = open(f'points_{k}.txt', 'r')
    lines = [list(map(float, line.strip(" ").split(" "))) for line in f]

    steps = np.zeros(len(lines))
    errors = np.zeros(len(lines))
    for j in range(0, len(lines)):
        steps[j] = np.log(10**lines[j][0])
        errors[j] = np.log(abs(lines[j][2]))
        print(abs(lines[j][2]))
    
    plt.scatter(steps, errors, s = 25, alpha=0.6, label=f'N = {k}')

plt.legend()
plt.ylabel("$\ln(error)$")
plt.xlabel("$\ln(step)$")
plt.grid(True, linestyle="--")
plt.savefig("graph.png")
plt.show()
