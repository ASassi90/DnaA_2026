
import numpy as np
import matplotlib.pyplot as plt

# Parameters
K = 1
z = 10
a = 1000
n_star = 975

# Define V*(alpha)
def V_star(alpha):
    return (alpha * n_star * z) / ((alpha * a - z)  * (alpha * K + z))

# Alpha domain that shows behavior and skips poles
pole1 = z / a      # 0.01
pole2 = -z / K     # -10
alpha_min, alpha_max = 1.1*z/a, 1.


alpha = np.linspace(alpha_min, alpha_max, 5000)

V_plot = V_star(alpha)

def alpha_of_V(volume):
    return 1.-0.9/volume

# Plot

volumes=np.linspace(0.1, 5., 1000)
alphas=alpha_of_V(volumes)
plt.figure(figsize=(8, 5))
plt.plot(alpha, V_plot, color='navy', lw=2)
plt.plot(alphas, volumes, color='red', lw=2)


plt.xlim(0., 1.)
plt.xlabel(r"$\alpha$")
plt.ylabel(r"$V^*(\alpha)$")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
