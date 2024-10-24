import matplotlib.pyplot as plt
import numpy as np

f_4 = open(f'lab-4-data-4.txt', 'r')

hs_4 = []
errs_4 = []
for line in f_4:
    hs_4.append(float((line.split("\t"))[0]))
    errs_4.append(float((line.split("\t"))[-1].replace("\n", "")))

f_6 = open(f'lab-4-data-6.txt', 'r')

hs_6 = []
errs_6 = []
for line in f_6:
    hs_6.append(float((line.split("\t"))[0]))
    errs_6.append(float((line.split("\t"))[-1].replace("\n", "")))

hs_4 = np.log(np.array(hs_4))
errs_4 = np.log(np.array(errs_4))

coefs_4 = np.polyfit(hs_4[1:], errs_4[1:], 1)
xs_4 = np.array([-20, 5])

coefs_6 = np.polyfit(hs_6, errs_6, 1)
xs_6 = np.array([-20, 5])

plt.figure(figsize=(8, 6), dpi=100)

plt.scatter(hs_4, errs_4, color = "green", label = "N = 4")
plt.plot(xs_4, coefs_4[0]*xs_4 + coefs_4[1], color = "darkgreen", alpha = 0.6, ls = "--", label = f"$y = {round(coefs_4[0],1)}x + ({round(coefs_4[1], 1)})$")

plt.scatter(hs_6, errs_6, color = "purple", label = "N = 6")
plt.plot(xs_6, coefs_6[0]*xs_6 + coefs_6[1], color = "orchid", alpha = 0.6, ls = "--", label = f"$y = {round(coefs_6[0],1)}x + ({round(coefs_6[1], 1)})$")

plt.xlabel("$ln step$")
plt.ylabel("$ln error$")
plt.xlim([-12, 6])
plt.ylim([-17, 3])
plt.legend()
plt.grid(True, linestyle="--")
plt.savefig("lab-4-graph.png")
plt.show()