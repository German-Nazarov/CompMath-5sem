import numpy as np
import matplotlib.pyplot as plt

plt.figure(figsize=(8, 6), dpi=100)
for n in range(3, 6, 1):
    f = open(f'lab-2-interp_data-n{n}-u.txt', 'r')

    lengths = []
    errs = []

    for line in f:
        print(line.split("  "))
        lengths.append(float((line.split(" "))[0]))
        errs.append(float((line.split(" "))[-1].replace("\n", "")))
    
    plt.scatter(lengths, errs, s = 40, alpha=1, label=f'N = {n}')
    print(errs)

plt.legend()
plt.ylabel("$error$")
plt.xlabel("$length$")
plt.grid(True, linestyle="--")
plt.savefig("lab-2-interp-graph-u.png")
plt.show()


plt.figure(figsize=(8, 6), dpi=100)
for n in range(3, 6, 1):
    f = open(f'lab-2-interp_data-n{n}-ch.txt', 'r')

    lengths = []
    errs = []

    for line in f:
        print(line.split("  "))
        lengths.append(float((line.split(" "))[0]))
        errs.append(float((line.split(" "))[-1].replace("\n", "")))
    
    plt.scatter(lengths, errs, s = 40, alpha=1, label=f'N = {n}')
    print(errs)

plt.legend()
plt.ylabel("$error$")
plt.xlabel("$length$")
plt.grid(True, linestyle="--")
plt.savefig("lab-2-interp-graph-ch.png")
plt.show()