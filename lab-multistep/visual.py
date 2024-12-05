import matplotlib.pyplot as plt
import numpy as np

f = open(f'multi-data.txt', 'r')
#.replace("\n", "")
steps = []
oscill = []
linear = []
for line in f:
    steps.append(float((line.split("\t"))[0]))
    oscill.append(float((line.split("\t"))[1]))
    linear.append(float((line.split("\t"))[2]))
steps = np.array(steps)
oscill = np.array(oscill)
linear = np.array(linear)
plt.plot(np.log(steps), np.log(np.abs(oscill)), color = "royalblue", label = f"$y'' = -y$")
plt.plot(np.log(steps), np.log(np.abs(linear)), color = "orange", label = f"$y' = t^3$")
x = np.array([np.log(steps)[0], np.log(steps)[-1]])
# plt.plot(x, -1.25*x-32.5, color = "red", ls = "--", alpha = 0.6, label = r"$y = -1.25x - 32.5$")
# plt.plot(x, 4*x-4, color = "blue", ls = "--", alpha = 0.6, label = r"$y = 4x - 4$")
# plt.plot(x, 4*x-31.5, color = "red", ls = "--", alpha = 0.2, label = r"$y = 4x - 31.5$")
plt.xlabel(r"$ln(\text{step})$")
plt.ylabel(r"$ln(\text{error})$")
plt.legend()
plt.savefig("rk-graph-2.png")
plt.show()
