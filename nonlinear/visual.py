import matplotlib.pyplot as plt
import numpy as np


nums = [1, 2, 5, 8]
plt.figure(figsize=(8, 6), dpi=100)
for i in range(4):
    f = open(f'data0{nums[i]}.txt', 'r')

    steps = []
    errs = []
    for line in f:
        steps.append(float((line.split("\t"))[0]))
        errs.append(float((line.split("\t"))[1]))
    steps = np.array(steps)
    errs = np.array(errs)
    plt.plot((steps), np.log(errs), label = f"data ecc = 0.{nums[i]}")

plt.xlabel("step (i)")
plt.ylabel("ln(error)")
plt.grid(True, linestyle="--")
plt.legend()
plt.savefig("graph.png")
plt.show()


