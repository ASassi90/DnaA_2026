from src.simulation.cycle_updates import make_step
from src.utils.setup import initialize_n_nforks, get_n_star_n_forks, initial_tree
from matplotlib import pyplot as plt
import numpy as np

def log_state(history, **kwargs):
    for key, value in kwargs.items():
        history[key].append(value)

def run_simulation(cfg):
    """
        These simulations returns the values of the main quantities of interest (such as volume, no of sites, 
        no of DnaA-ATP proteins, no of origins etc.) as a function of time. 
        It can be used both in the case the firing rate is given by k=k_max*P_open and in the case we assume 
        perfect step-wise response. 
    """
    simulation_data={
        "time" : [],
        "volume" : [],
        "n_tot" : [], 
        "a_atp" : [],
        "a_adp" : [],
        "c_atp" : [],
        "c_adp" : [],
        "n_forks" : [],
        "origins" : [],
        "fpr" : []
    }
    """
    n_tot, n_forks, chi_0, y, tree_manager=initialize(dnaa=cfg.model.DNAA_CONCENTRATION, 
                                            n_sites=cfg.model.SITES, 
                                            kori=cfg.model.K_OPEN, 
                                            epsilon_cost=cfg.model.E_COST, 
                                            origin_sites=cfg.model.ORIGIN_SITES, 
                                            K=cfg.model.K,
                                            cfg=cfg)
    """
    n_tot, n_forks = initialize_n_nforks(cfg)
    tree_manager=initial_tree(cfg)
    tree_manager.n_forks=n_forks
    _, n_forks_init=get_n_star_n_forks(cfg)
    chi=cfg.model.CHI0/n_forks_init
    y=cfg.model.COOP
    t_max=cfg.simulation.T_MAX
    dt=cfg.simulation.DT
    volume, alpha, =1., 0.999
    a_atp, a_adp = alpha*cfg.model.DNAA_CONCENTRATION, (1.-alpha)*cfg.model.DNAA_CONCENTRATION
    time=0.
    count=1
    while time<t_max:
        if count%20000==0:
            print(f"{time/t_max:.3g}")
        time, a_atp, a_adp, c_atp, c_adp, volume, n_tot, f_rate = make_step(n_forks, n_tot, volume, a_atp, a_adp, time, 
                                                                          dt, y, chi, step=False, cfg=cfg)
        tree_manager.update(time, f_rate, volume, n_tot)
        tree_manager.simulate_step(cfg)
        n_tot = tree_manager.n_tot
        volume = tree_manager.volume
        n_forks = tree_manager.n_forks
        log_state(simulation_data,
                  time=time, 
                  volume=volume, 
                  n_tot=n_tot, 
                  a_atp=a_atp,
                  a_adp=a_adp, 
                  c_atp=c_atp, 
                  c_adp=c_adp, 
                  fpr=f_rate,
                  origins=len(tree_manager.origins.keys()),
                  n_forks=n_forks
        )
        count+=1
    return simulation_data

def initiation_and_division(initiation_times, division_times, ax, cfg):
    """
        This is just needed for the plots. 
    """
    for in_time in initiation_times:
        ax.axvspan(in_time - cfg.model.LICENSING,  # left boundary
                    in_time,                          # right boundary
                    color='k',                        # fill color (black)
                    alpha=0.2,                        # transparency
                    linewidth=0)                      # no edge line

    for div_time in division_times:
        ax.axvline(div_time,linestyle='--', color='k')

def main(cfg):
    simulation_data=run_simulation(t_max=200000., dt=0.1)
    time=simulation_data["time"]
    origins=np.array(simulation_data["origins"])
    sites=np.array(simulation_data["n_tot"])
    volume=np.array(simulation_data["volume"])
    a_ATP=simulation_data["a_atp"]
    a_ADP=simulation_data["a_adp"]
    f_rate=simulation_data["fpr"]
    division_times=[]
    initiation_times=[]
    for volume_index in range(1, len(volume)):
        if volume[volume_index]<volume[volume_index-1]:
            division_times.append(time[volume_index-1])
            initiation_times.append(time[volume_index-1]-60.)
 
    start=int(len(volume)/50)
    end=start+4000
    starting_time=time[start]
    time =np.array(time)-starting_time
    initiation_times =np.array(initiation_times)-starting_time
    division_times =np.array(division_times)-starting_time

    colors = [
        "#1f77b4",  # Blue
        "#ff7f0e",  # Orange
        "#2ca02c",  # Green
        "#d62728",  # Red
        "#9467bd"   # Purple
        ]

    #fig, axs = create_figure(layout='stacked', figsize=(6,8), n_stacked=5, sharex=True, sharey=False, 
    #              xlim=(time[start], time[end]), xticks=None, yticks=None, labelsize=14)
    

    a_tot=np.array(a_ATP)+np.array(a_ADP)    
    alphas=np.array(a_ATP)/a_tot
    """
    axs[0].plot(time, sites, '--', label="sites", color=colors[0])
    axs[0].plot(time, volume*a_tot, label="DnaA proteins", color=colors[0])
    axs[0].set_ylim(300., 1500.)
    axs[0].set_ylabel("copy numbers")
    initiation_and_division(initiation_times, division_times, axs[0])
    axs[0].legend(frameon=True, loc='upper right', framealpha=1, facecolor='white', edgecolor='black')

    conc_ratio=(sites/volume)/cfg.model.DNAA_CONCENTRATION
    axs[1].plot(time, conc_ratio, color=colors[1])
    axs[1].axhline(1., color='k')
    initiation_and_division(initiation_times, division_times, axs[1])
    axs[1].set_ylim(0.9, 1.2)
    axs[1].set_ylabel(r"$c_{tot}/a$")
    
    axs[2].plot(time, origins, color=colors[2])   
    axs[2].set_ylim(-0.1, 10.)
    axs[2].set_ylabel("origins")
    initiation_and_division(initiation_times, division_times, axs[2])

    axs[3].plot(time, alphas, color=colors[3])
    axs[3].set_ylim(-0.05, 0.3)
    axs[3].set_ylabel(rf"$\alpha$")
    initiation_and_division(initiation_times, division_times, axs[3])

    axs[4].plot(time, f_rate, color=colors[4])   
    axs[4].set_ylim(-5., 101.)
    initiation_and_division(initiation_times, division_times, axs[4])
    axs[4].set_ylabel(r"$P_{open}k_{max}$")

    """
    fig_volume = plt.figure("volume")
    plt.plot(time, volume) 
    a_tot=np.array(a_ATP)+np.array(a_ADP)
    alphas=np.array(a_ATP)/a_tot
    plt.xlim(time[start], time[end])
    #plt.ylim(-0.05, 1.5)
    #fig_volume.savefig(rf"C:\Users\Albi\Desktop\DnaA_manuscript\figures\volume.png")

    fig_alpha = plt.figure("alpha")
    plt.plot(time, alphas)    
    plt.xlim(time[start], time[end])
    plt.ylim(-0.02, 0.4)
    #fig_alpha.savefig(rf"C:\Users\Albi\Desktop\DnaA_manuscript\figures\alpha.png")


    plt.figure("origins")
    plt.plot(time, origins)   
    plt.xlim(time[start], time[end])
    plt.ylim(-0.1, 10.)

    plt.figure("sites and initiators")
    plt.plot(time, sites)
    plt.plot(time, volume*a_tot)
    #initiation_and_division(initiation_times, division_times)
    plt.xlim(time[start], time[end])
    plt.ylim(-5., 2000)


    plt.figure("firing rate")
    plt.plot(time, f_rate)
    #initiation_and_division(initiation_times, division_times)
    plt.xlim(time[start], time[end])


    plt.figure("c_tot")
    conc_ratio=(sites/volume)/cfg.model.DNAA_CONCENTRATION
    plt.plot(time, conc_ratio)
    plt.axhline(1., color='k')
    #initiation_and_division(initiation_times, division_times)
    plt.xlim(time[start], time[end])
    plt.ylim(0.5, 1.5)
    
    #axs[4].set_xlabel("time (min)")
    #fig.savefig(rf"C:\Users\Albi\Desktop\DnaA_manuscript\figures\time_traces.png")
    plt.show()

if __name__=="__main__":
    main()