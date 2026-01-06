
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
from src.utils.helpers import create_figure

# Define the function and solver
def h(t2, t1, C, lamb, gamma):
    rho = min(1, t1 / C)
    LHS = (1 + rho+gamma) * np.exp(lamb * t2)
    new_rho = min(1, t2 / C)
    RHS = 2 * (1 + new_rho+gamma)

    return LHS - RHS

def solve_t2(t1, C, lamb, gamma, t3_bounds=(1e-6, 100)):
    try:
        result = root_scalar(h, args=(t1, C, lamb, gamma), bracket=t3_bounds, method='brentq')
        return result.root if result.converged else 0
    except ValueError:
        return 0

# Parameters
t1_i=35.
t1 = t1_i
C = 40.
tau=25.
lamb = np.log(2.) / tau
n_iter = 15
gamma=0.6

# Collect cobweb points
points = []
for _ in range(n_iter):
    t2 = solve_t2(t1, C, lamb, gamma)
    points.append((t1, t2))  # (x_n, x_{n+1})
    t1=t2

# Plotting
fig, ax = create_figure(layout='single', figsize=(6,6), xlabel=r'$t_a$', ylabel=r'$t_b$', labelsize=30)
plt.subplots_adjust(bottom=0.15, top=0.85, left=0.15, right=0.85)
# Identity line
#x_vals = np.linspace(min(p[0] for p in points) - 1, max(p[0] for p in points) + 1, 100)
x_vals=np.linspace(0., 80., 100)
ax.plot(x_vals, x_vals, 'k--', label='$x_{n+1} = x_n$')

# Helper to draw oriented segment with small triangle
def draw_oriented_segment(x0, y0, x1, y1, color):
    # Draw the line segment
    ax.plot([x0, x1], [y0, y1], color=color, linewidth=1.8)

    # Compute midpoint
    xm, ym = (x0 + x1) / 2, (y0 + y1) / 2
    dx, dy = x1 - x0, y1 - y0

    # Normalize direction
    norm = np.hypot(dx, dy)
    if norm == 0:
        return
    dx /= norm
    dy /= norm

    # Small offset to make the triangle visible
    scale = 0.5  # triangle size
    ax.plot(xm, ym, marker=(3, 0, np.degrees(np.arctan2(-dx, dy))), 
            markersize=8, color=color)

print(points)
xx=np.linspace(0., 80., 1000)
yy=[solve_t2(x, C, lamb, gamma) for x in xx]
ax.plot(xx, yy, color='b', lw=2)

gamma_threshold=tau/C*(2./np.log(2.)-1.)-1.

yy_th=[solve_t2(x, C, lamb, gamma=gamma_threshold) for x in xx]
ax.plot(xx, yy_th, '-', color='black', lw=2)

yy0=[solve_t2(x, C, lamb, gamma=0) for x in xx]
ax.plot(xx, yy0, '--', color='b', lw=2)



ax.set_xlim(5., 45.)
ax.set_ylim(5., 45.)

# Draw cobweb as oriented segments
x_prev = points[0][0]
for (x_n, x_np1) in points:
    # Vertical: (x_prev, x_prev) -> (x_prev, x_np1)
    draw_oriented_segment(x_prev, x_prev, x_prev, x_np1, color='r')
    # Horizontal: (x_prev, x_np1) -> (x_np1, x_np1)
    draw_oriented_segment(x_prev, x_np1, x_np1, x_np1, color='r')
    x_prev = x_np1

# Final styling
#ax.set_xlabel(rf'$t_a$')
#ax.set_ylabel(rf'$t_b$')
#ax.set_title('Cobweb Plot with Oriented Segments')
ratio=C/tau
#ax.set_title(r"$\frac{{C}}{{\tau}} = {:.2f}$".format(ratio))
ax.set_aspect('equal', adjustable='box')
ax.grid(True)
#ax.legend()
fig.savefig(rf'C:\Users\Albi\Desktop\Pigolotti\DnaA\figures\cobweb_{tau}.jpeg')
plt.tight_layout()
plt.show()
