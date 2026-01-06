
import numpy as np
import matplotlib.pyplot as plt

# --- Parameters (from user) ---
C = 40.0
tau = 25.0
Lambda = np.log(2.0) / tau
n_genomes = 2
n_sites = 300

# --- V* coefficients (from user) ---
a = 1000.0
z = 10.0

# --- chi(t) = chi0 * n_forks(t); set chi0 to your desired value ---
#chi0 = 0.22  # placeholder; adjust if needed


COLORS = {
    "V_star": "#0072B2",      # blue
    "V": "#D55E00",           # vermillion
    "fill": "#FFBB78",        # light orange (use alpha)
    "init_line": "#009E73",   # green
    "term_line": "#CC79A7",   # purple
    "text": "#6E6E6E"         # slate gray
}



for chi0 in (0.,  0.22):
    # --- Time grid requested: [-10, +30] ---
    T_min = -10.0
    T_max = 30.0
    N = 1601
    t = np.linspace(T_min, T_max, N)

    # --- Fork schedule with added change at t = tau ---
    C_minus_tau = C - tau  # 15.0
    # Right-continuous steps:
    # 4 for t<0, 12 for 0<=t<15, 8 for 15<=t<tau, 24 for t>=tau
    n_forks = np.where(t < 0.0, 4,
                np.where(t < C_minus_tau, 12,
                np.where(t < tau, 8, 24)))

    # --- n_star(t) ---
    # n_star(t) = n_genomes * n_sites * (1 + min(1, (tau + t)/C) + 2*max(0, t/C))
    term1 = np.minimum(1.0, (tau + t) / C)
    term2 = np.minimum(1.0, np.maximum(0.0, t / C))
    term3 = np.minimum(1.0, np.maximum(0.0, (t-tau) / C))

    n_star = n_genomes * n_sites * (1.0 + term1 + 2.0 * (term2+2.0*term3))

    # --- chi(t) ---
    chi = chi0 * n_forks

    # --- V*(t) ---
    # V* = ((n_star + a * chi) * (1 + sqrt(z/a))) / (2 * (a - z))
    if a <= 0 or a <= z:
        raise ValueError("Parameters must satisfy a > 0 and a > z for V* to be finite.")
    V_star = ((n_star + a * chi) * (1.0 + np.sqrt(z / a))) / (2.0 * (a - z))

    # --- Set V(0)=V*(0^-) BEFORE fork change ---
    n_star_0 = n_genomes * n_sites * (1.0 + min(1.0, (tau + 0.0) / C) + 2.0 * max(0.0, 0.0 / C))
    chi_pre0 = chi0 * 4  # forks before initiation
    V_star_pre0 = ((n_star_0 + a * chi_pre0) * (1.0 + np.sqrt(z / a))) / (2.0 * (a - z))

    V0 = V_star_pre0
    V = V0 * np.exp(Lambda * t)
    V/=V0
    V_star/=V0
    # --- Plot ---
    fig, ax = plt.subplots(figsize=(9.5, 4.5))
    ax.plot(t, V_star, label=r"$V^*(t)/V(0)$", color=COLORS["V_star"], linewidth=2.6)
    ax.plot(t, V, label=r"$V(t)/V(0)$", color=COLORS["V"], linewidth=2.0, linestyle="--")

    # Shade area between curves from t=0 to t=tau
    mask = (t >= 0.0) & (t <= tau)
    ax.fill_between(t[mask], V_star[mask], V[mask], color=COLORS["fill"], alpha=0.35)
    
    # Mark fork change events within the plotting window
    for change_time, color in [(0.0, COLORS["init_line"]), (C_minus_tau, COLORS["term_line"]), (tau, COLORS["init_line"])]:
        if T_min <= change_time <= T_max:
            ax.axvline(change_time, color=color, linestyle=":", linewidth=2.)
            ymax = plt.ylim()[1]
            #ax.text(change_time, ymax*0.96, label_text, rotation=90,
            #        va='top', ha='right', fontsize=9, color="#7f7f7f")

    
    # Annotate the t<0 region
    #plt.text(T_min + 0.5, plt.ylim()[1]*0.92, "t < 0: forks = 4", fontsize=9, color="#7f7f7f")
    ax.set_xlabel("time (min)", fontsize=12)
    ax.set_ylabel("normalized volume", fontsize=12)
    ax.set_xlim(T_min, T_max)
    ax.set_ylim(0.5, 2.5)
    ax.legend(loc="upper left", frameon=False)
    ax.grid(True, alpha=0.25)
    #fig.tight_layout()
    fig.savefig(rf"C:\Users\Albi\Desktop\Pigolotti\DnaA\figures\chi_%g.png"%chi0)
    plt.clf()
    plt.close('all')

plt.show()
